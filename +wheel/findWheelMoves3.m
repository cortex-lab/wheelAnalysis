

function [moveOnsets moveOffsets] = findWheelMoves3(pos, t, Fs, params)
% function [moveOnsets moveOffsets] = findWheelMoves3(pos, t, Fs, params)
%
% algorithm is: for each point, is there >posThresh max movement in the next 
% tThresh seconds. If there is, then that tThresh window is part of a
% movement. Merge small gaps. Now for every time you go from not-moving to
% moving, jump ahead by tThresh and look backwards in time until you find a
% point that's very close to the starting point (different by
% <posThreshOnset). Finally, drop movements that are too brief. 

posThresh = 8; % if position changes by less than this
tThresh = 0.2; % over at least this much time, then it is a quiescent period
minGap = 0.1; % any movements that have this little time between the end of one and
    % the start of the next, we'll join them
posThreshOnset = 1.5; % a lower threshold, used when finding exact onset times.     
minDur = 0.05; % seconds, movements shorter than this are dropped
makePlots = false;

if ~isempty(params)
    fn = fieldnames(params);
    for f = 1:length(fn)
        eval(sprintf('%s = params.%s;', fn{f}, fn{f}));
    end
end

% first compute an evenly-sampled position (in case that's not what was
% provided, if it is, that's ok this is fast) and velocity
rawT = t; rawPos = pos; 
t = rawT(1):1/Fs:rawT(end);
pos = interp1(rawT, rawPos, t);


tThreshSamps = round(tThresh*Fs);
wheelToep = toeplitz(fliplr(pos), nan(1,tThreshSamps));
wheelToep = flipud(wheelToep);
totalDev = max(wheelToep,[],2)-min(wheelToep,[],2);
isMoving = totalDev>posThresh;
isMoving = [false; isMoving]; 
isMoving(end) = false; % make sure we end on an offset

% fill in small gaps - this is like a dilation/contraction
moveOnsetSamps = find(~isMoving(1:end-1) & isMoving(2:end));
moveOffsetSamps = find(isMoving(1:end-1) & ~isMoving(2:end));
tooShort = find((moveOnsetSamps(2:end)-moveOffsetSamps(1:end-1))/Fs<minGap);
for q = 1:length(tooShort)
    isMoving(moveOffsetSamps(tooShort(q)):moveOnsetSamps(tooShort(q)+1)) = true;
end

% definition of move onset: 
% - starting from the end of the tThresh window, look back until you find
% one that's not different from the moveOnset, by some smaller threshold.
moveOnsetSamps = find(~isMoving(1:end-1) & isMoving(2:end));
wheelToepOnsets = wheelToep(moveOnsetSamps,:);
wheelToepOnsetsDev = abs(bsxfun(@minus, wheelToepOnsets, wheelToepOnsets(:,1)));
hasOnset = wheelToepOnsetsDev>posThreshOnset;
[a, b] = find(~fliplr(hasOnset));
onsetLags = tThreshSamps-accumarray(a, b, [], @min);
moveOnsetSamps = moveOnsetSamps+onsetLags;
moveOnsets = t(moveOnsetSamps);

% we won't do the same thing for offsets, instead just take the actual end
% of the isMoving. This is because we're just not so concerned about being
% temporally precise with these. 
moveOffsetSamps = find(isMoving(1:end-1) & ~isMoving(2:end));
moveOffsets = t(moveOffsetSamps);

moveDurs = moveOffsets-moveOnsets;
tooShort = moveDurs<minDur;
moveOnsetSamps = moveOnsetSamps(~tooShort);
moveOnsets = moveOnsets(~tooShort); 
moveOffsetSamps = moveOffsetSamps(~tooShort);
moveOffsets = moveOffsets(~tooShort); 

moveGaps = moveOnsets(2:end)-moveOffsets(1:end-1);
gapTooSmall = moveGaps<minGap;
% for these, drop the offending offset and onset, which effectively joins
% the two
moveOnsets = moveOnsets([true ~gapTooSmall]); % always keep first onset
moveOnsetSamps = moveOnsetSamps([true ~gapTooSmall]); 
moveOffsets = moveOffsets([~gapTooSmall true]); % always keep last offset
moveOffsetSamps = moveOffsetSamps([~gapTooSmall true]); % always keep last offset

moveOnsets = moveOnsets(:); % return a column
moveOffsets = moveOffsets(:); % return a column

% see how it looks
if makePlots
    figure; 
    
    ax1 = subplot(2,1,1);
    %plot(t, pos); 
    hold on; 
    plot(moveOnsets, pos(moveOnsetSamps), 'go');
    plot(moveOffsets, pos(moveOffsetSamps), 'bo');
    hold on; 
    inMove = logical(WithinRanges(t, [moveOnsets; moveOffsets]'));
    plot(t(inMove), pos(inMove), 'r.');
    plot(t(~inMove), pos(~inMove), 'k.');
    ylabel('position');
    
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