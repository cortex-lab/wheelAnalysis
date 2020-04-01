function [onsets, offsets, displacement, peakVelTimes, peakAmps] = ...
  findWheelMoves3(pos, t, Fs, varargin)
% [onsets, offsets, s, peakVelTimes, peakAmps] = findWheelMoves3(pos, t, Fs, params)
%
% Algorithm: for each point, is there > posThresh max movement in the
% next tThresh seconds. If there is, then that tThresh window is part of a
% movement. Merge small gaps. Now for every time you go from not-moving to
% moving, jump ahead by tThresh and look backwards in time until you find a
% point that's very close to the starting point (different by <
% posThreshOnset). Finally, drop movements that are too brief.
% 
% Required Inputs:
%   pos : an array of wheel positions
%   t : an array of wheel sample timestamps 
%   Fs : the sampling frequency used for linear interpolation
%
% Optional Parameters (may be struct or name-value pairs): 
%   posThresh = 8 : if position changes by less than this
%   tThresh = 0.2 : over at least this much time, then it is a quiescent period
%   minGap = 0.1 : any movements that have this little time between the end 
%     of one and the start of the next, we'll join them
%   posThreshOnset = 1.5 : a lower threshold, used when finding exact onset times.     
%   minDur = 0.05 : seconds, movements shorter than this are dropped.
%   makePlots = false : plot position and velocity showing detected movements.
%   batchSize = 10000 : compute in batches of this size.  The larger the 
%     matrix the higher the memory use, but not by much.  Must be >= length(pos).
%
% Outputs:
%   onsets : an array of detected movement onset times
%   offsets : an array of detected movement offset times
%   displacement : the total displacement of each movement
%   peakVelTimes : the time of peak velocity for each detected movement
%   peakAmplitude : the absolute maximum amplitude of each detected 
%     movement, relative to onset position.

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
onsetSamps = find(~isMoving(1:end-1) & isMoving(2:end));
offsetSamps = find(isMoving(1:end-1) & ~isMoving(2:end));
tooShort = find((onsetSamps(2:end) - offsetSamps(1:end-1)) / Fs < p.minGap);
for q = 1:length(tooShort)
    isMoving(offsetSamps(tooShort(q)):onsetSamps(tooShort(q)+1)) = true;
end

%% Compute precise movement onset samples
% definition of move onset: 
% - starting from the end of the tThresh window, look back until you find
% one that's not different from the moveOnset, by some smaller threshold.
onsetSamps = find(~isMoving(1:end-1) & isMoving(2:end));

% Same as above, this is to replace the two lines below in a memory
% controlled way without looping on every sample
%  wheelToepOnsets = wheelToep(moveOnsetSamps,:);
%  wheelToepOnsetsDev = abs(bsxfun(@minus, wheelToepOnsets, wheelToepOnsets(:,1)));
wheelToepOnsetsDev = zeros(length(onsetSamps), tThreshSamps);
c = 0; cwt = 0;
while ~isempty(onsetSamps)
    i2proc = (1:p.batchSize) + c;
    [icomm] = intersect( i2proc(1:end-tThreshSamps-1), onsetSamps );
    [~, itpltz] = intersect( i2proc(1:end-tThreshSamps-1), onsetSamps );
    i2proc(i2proc > length(t)) = [];
    if ~isempty(icomm)        
        w2e = hankel(pos(i2proc), nan(1,tThreshSamps));
        w2e = abs(bsxfun(@minus, w2e, w2e(:,1)));
        wheelToepOnsetsDev(cwt + (1:length(icomm)), :) =  w2e(itpltz, :);
        cwt = cwt + length(icomm);
    end
    c = c + p.batchSize - tThreshSamps;
    if i2proc(end) >= onsetSamps(end), break, end
end
warning('on', 'MATLAB:hankel:AntiDiagonalConflict')

hasOnset = wheelToepOnsetsDev > p.posThreshOnset;
[a, b] = find(~fliplr(hasOnset));
onsetLags = tThreshSamps-accumarray(a(:), b(:), [], @min);
onsetSamps = onsetSamps + onsetLags;
onsets = t(onsetSamps);

% we won't do the same thing for offsets, instead just take the actual end
% of the isMoving. This is because we're just not so concerned about being
% temporally precise with these. 
offsetSamps = find(isMoving(1:end-1) & ~isMoving(2:end));
offsets = t(offsetSamps);

moveDurs = offsets - onsets;
tooShort = moveDurs < p.minDur;
onsetSamps = onsetSamps(~tooShort);
onsets = onsets(~tooShort); 
offsetSamps = offsetSamps(~tooShort);
offsets = offsets(~tooShort); 

moveGaps = onsets(2:end) - offsets(1:end-1);
gapTooSmall = moveGaps < p.minGap;
% for these, drop the offending offset and onset, which effectively joins
% the two
if ~isempty(onsets)
    onsets = onsets([true ~gapTooSmall]); % always keep first onset
    onsetSamps = onsetSamps([true ~gapTooSmall]);
    offsets = offsets([~gapTooSmall true]); % always keep last offset
    offsetSamps = offsetSamps([~gapTooSmall true]); % always keep last offset
end
onsets = onsets(:); % return a column
offsets = offsets(:); % return a column
% Calculate displacement
displacement = pos(offsetSamps) - pos(onsetSamps);
displacement = displacement(:); % return a column
% Calculate peak velocity times and peak amplitudes
vel = conv(diff([0 pos]), wheel.gausswin(10), 'same');
peakVelTimes = nan(size(onsets));
peakAmps = nan(size(onsets));
for m = 1:numel(onsets)
    thisV = abs(vel(onsetSamps(m):offsetSamps(m)));
    peakVelTimes(m) = onsets(m) + find(thisV == max(thisV), 1) / Fs;
    % Get index of maximum absolute position relative to move onset
    [~,I] = max(abs(pos(onsetSamps(m):offsetSamps(m)) - pos(onsetSamps(m))));
    peakAmps(m) = pos(onsetSamps(m) + I) - pos(onsetSamps(m));
end

%% see how it looks
if p.makePlots
    figure('Name', 'Wheel movements'); 
    % Plot the wheel position
    ax1 = subplot(2,1,1);
    hold on; 
    on = plot(onsets, pos(onsetSamps), 'go', 'DisplayName', 'onset');
    off = plot(offsets, pos(offsetSamps), 'bo', 'DisplayName', 'offset');
    hold on; 
    inMove = logical(WithinRanges(t, [onsets offsets]));
    in = plot(t(inMove), pos(inMove), 'r.', 'DisplayName', 'in movement');
    plot(t(~inMove), pos(~inMove), 'k.');
    ylabel('position');
    legend([on off in], 'Location', 'SouthEast')
    
    % Plot the velocity trace
    ax2 = subplot(2,1,2);
    vel = wheel.computeVelocity2(pos, 0.015, Fs);
    hold on; 
    plot(onsets, vel(onsetSamps), 'go');
    plot(offsets, vel(offsetSamps), 'bo');
    plot(t(inMove), vel(inMove), 'r.');
    plot(t(~inMove), vel(~inMove), 'k.');
    ylabel('velocity');
    xlabel('time (sec)');
    
    linkaxes([ax1 ax2], 'x');
end