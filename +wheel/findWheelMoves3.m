function [moveOnsets, moveOffsets, moveAmps, peakVelTimes] = findWheelMoves3(pos, t, Fs, varargin)
% [moveOnsets, moveOffsets] = findWheelMoves3(pos, t, Fs, params)
%
% algorithm is: for each point, is there > posThresh max movement in the
% next tThresh seconds. If there is, then that tThresh window is part of a
% movement. Merge small gaps. Now for every time you go from not-moving to
% moving, jump ahead by tThresh and look backwards in time until you find a
% point that's very close to the starting point (different by <
% posThreshOnset). Finally, drop movements that are too brief.
% 
% params may be structure or name-value pairs and can include:
% - posThresh = 8; % if position changes by less than this
% - tThresh = 0.2; % over at least this much time, then it is a quiescent period
% - minGap = 0.1; % any movements that have this little time between the end of one and
% % the start of the next, we'll join them
% - posThreshOnset = 1.5; % a lower threshold, used when finding exact onset times.     
% - minDur = 0.05; % seconds, movements shorter than this are dropped
% - makePlots = false; % plot position and velocity showing detected movements
% - batchSize = 10000; % compute in batches of this size.  The lager the matrix the higher
% % the memory use, but not by much.

%% Validate input
if numel(varargin) == 1 && isempty(varargin{1})
  varargin(1) = []; % Backward compatility for old param input
end
p = inputParser;
validator = @(x) isnumeric(x) && isscalar(x) && (x > 0);
addRequired(p, 'pos', @(x) isnumeric(x) && isvector(x));
addRequired(p, 't', @(x) isnumeric(x) && isvector(x));
addRequired(p, 'Fs', validator);
% If position changes by less than this...
addParameter(p, 'posThresh', 8, validator);
% over at least this much time, then it is a quiescent period.
addParameter(p, 'tThresh', 0.2, validator);
% Any movements that have this little time between the end of one and the
% start of the next, we'll join them.
addParameter(p, 'minGap', 0.1, validator);
% A lower threshold, used when finding exact onset times.
addParameter(p, 'posThreshOnset', 1.5, validator);
% Minimum duration in second.  Movements shorter than this are dropped.
addParameter(p, 'minDur', 0.05, validator);
% Compute matrix in batches of this size.  The lager the matrix the higher
% the memory use.
addParameter(p, 'batchSize', 10000, validator);
% Plot position and velocity showing detected movements.
addParameter(p, 'makePlots', false, @(x) islogical(x) && isscalar(x));
parse(p, pos, t, Fs, varargin{:});

p = p.Results; % Final parameters

% First compute an evenly-sampled position (in case that's not what was
% provided, if it is, that's ok this is fast)
rawT = t; rawPos = pos;
t = rawT(1):1/Fs:rawT(end);
pos = interp1(rawT, rawPos, t);

% Convert the time threshold into number of samples given the sampling
% frequency
tThreshSamps = round(p.tThresh*Fs);

%% Compute the approximate movement onset and offset samples
% Computing the full Toeplitz/Hankel matrix blows up memory for long
% timeseries so we will do it in batches.  This is quite ugly, but memory
% controlled and reasonably fast. This is stricly equivalent to the 3 lines
% commented below:
%   wheelToep = toeplitz(fliplr(pos), nan(1,tThreshSamps));
%   wheelToep = flipud(wheelToep);
%   totalDev = max(wheelToep,[],2)-min(wheelToep,[],2);
warning('off', 'MATLAB:hankel:AntiDiagonalConflict')
totalDev = zeros(length(t), 1); % Initialize vector of total deviations
c = 0;
while true
    i2proc = (1:p.batchSize) + c;
    i2proc(i2proc > length(t)) = [];
    w2e = hankel(pos(i2proc), nan(1, tThreshSamps));
    totalDev(i2proc) = max(w2e, [], 2) - min(w2e, [], 2);
    c = c + p.batchSize - tThreshSamps;
    if i2proc(end) == length(t), break, end
end

isMoving = totalDev > p.posThresh;
isMoving = [false; isMoving]; 
isMoving(end) = false; % make sure we end on an offset

% fill in small gaps - this is like a dilation/contraction
moveOnsetSamps = find(~isMoving(1:end-1) & isMoving(2:end));
moveOffsetSamps = find(isMoving(1:end-1) & ~isMoving(2:end));
tooShort = find((moveOnsetSamps(2:end) - moveOffsetSamps(1:end-1)) / Fs < p.minGap);
for q = 1:length(tooShort)
    isMoving(moveOffsetSamps(tooShort(q)):moveOnsetSamps(tooShort(q)+1)) = true;
end

%% Compute precise movement onset samples
% definition of move onset: 
% - starting from the end of the tThresh window, look back until you find
% one that's not different from the moveOnset, by some smaller threshold.
moveOnsetSamps = find(~isMoving(1:end-1) & isMoving(2:end));

% Same as above, this is to replace the two lines below in a memory
% controlled way without looping on every sample
%  wheelToepOnsets = wheelToep(moveOnsetSamps,:);
%  wheelToepOnsetsDev = abs(bsxfun(@minus, wheelToepOnsets, wheelToepOnsets(:,1)));
wheelToepOnsetsDev = zeros(length(moveOnsetSamps), tThreshSamps);
c = 0; cwt = 0;
while ~isempty(moveOnsetSamps)
    i2proc = (1:p.batchSize) + c;
    [icomm] = intersect( i2proc(1:end-tThreshSamps-1), moveOnsetSamps );
    [~, itpltz] = intersect( i2proc(1:end-tThreshSamps-1), moveOnsetSamps );
    i2proc(i2proc > length(t)) = [];
    if ~isempty(icomm)        
        w2e = hankel(pos(i2proc), nan(1,tThreshSamps));
        w2e = abs(bsxfun(@minus, w2e, w2e(:,1)));
        wheelToepOnsetsDev(cwt + (1:length(icomm)), :) =  w2e(itpltz, :);
        cwt = cwt + length(icomm);
    end
    c = c + p.batchSize - tThreshSamps;
    if i2proc(end) >= moveOnsetSamps(end), break, end
end
warning('on', 'MATLAB:hankel:AntiDiagonalConflict')

hasOnset = wheelToepOnsetsDev > p.posThreshOnset;
[a, b] = find(~fliplr(hasOnset));
onsetLags = tThreshSamps-accumarray(a(:), b(:), [], @min);
moveOnsetSamps = moveOnsetSamps + onsetLags;
moveOnsets = t(moveOnsetSamps);

% we won't do the same thing for offsets, instead just take the actual end
% of the isMoving. This is because we're just not so concerned about being
% temporally precise with these. 
moveOffsetSamps = find(isMoving(1:end-1) & ~isMoving(2:end));
moveOffsets = t(moveOffsetSamps);

moveDurs = moveOffsets - moveOnsets;
tooShort = moveDurs < p.minDur;
moveOnsetSamps = moveOnsetSamps(~tooShort);
moveOnsets = moveOnsets(~tooShort); 
moveOffsetSamps = moveOffsetSamps(~tooShort);
moveOffsets = moveOffsets(~tooShort); 

moveGaps = moveOnsets(2:end) - moveOffsets(1:end-1);
gapTooSmall = moveGaps < p.minGap;
% for these, drop the offending offset and onset, which effectively joins
% the two
if ~isempty(moveOnsets)
    moveOnsets = moveOnsets([true ~gapTooSmall]); % always keep first onset
    moveOnsetSamps = moveOnsetSamps([true ~gapTooSmall]);
    moveOffsets = moveOffsets([~gapTooSmall true]); % always keep last offset
    moveOffsetSamps = moveOffsetSamps([~gapTooSmall true]); % always keep last offset
end
moveOnsets = moveOnsets(:); % return a column
moveOffsets = moveOffsets(:); % return a column

moveAmps = pos(moveOffsetSamps) - pos(moveOnsetSamps);
moveAmps = moveAmps(:); % return a column
vel = conv(diff([0 pos]), wheel.gausswin(10), 'same');
peakVelTimes = nan(size(moveOnsets));
for m = 1:numel(moveOnsets)
    thisV = abs(vel(moveOnsetSamps(m):moveOffsetSamps(m)));
    peakVelTimes(m) = moveOnsets(m) + find(thisV == max(thisV), 1) / Fs;
end
%% see how it looks
if p.makePlots
    figure('Name', 'Wheel movements'); 
    % Plot the wheel position
    ax1 = subplot(2,1,1);
    hold on; 
    on = plot(moveOnsets, pos(moveOnsetSamps), 'go', 'DisplayName', 'onset');
    off = plot(moveOffsets, pos(moveOffsetSamps), 'bo', 'DisplayName', 'offset');
    hold on; 
    inMove = logical(WithinRanges(t, [moveOnsets moveOffsets]));
    in = plot(t(inMove), pos(inMove), 'r.', 'DisplayName', 'in movement');
    plot(t(~inMove), pos(~inMove), 'k.');
    ylabel('position');
    legend([on off in], 'Location', 'SouthEast')
    
    % Plot the velocity trace
    ax2 = subplot(2,1,2);
    vel = wheel.computeVelocity2(pos, 0.015, Fs);
    hold on; 
    plot(moveOnsets, vel(moveOnsetSamps), 'go');
    plot(moveOffsets, vel(moveOffsetSamps), 'bo');
    plot(t(inMove), vel(inMove), 'r.');
    plot(t(~inMove), vel(~inMove), 'k.');
    ylabel('velocity');
    xlabel('time (sec)');
    
    linkaxes([ax1 ax2], 'x');
end