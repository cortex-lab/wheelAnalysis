

function moveType = classifyWheelMoves(t, pos, moveOnsets, moveOffsets, intStartTime, respTime, resp)
% function moveType = classifyWheelMoves(t, pos, moveOnsets, moveOffsets, intStartTime, respTime, resp)
%
% Label movements as "left" (1), "right" (2), "flinch" (0), or unclassified
% (3). 
%
% intStartTime is the time that interactive period started
% respTime is the choiceworld-defined response time
% resp is the response type (1, 2, or 3)

crossWindow = 0.5; % time you have to cross the threshold in order to count

moveOnsets = moveOnsets(:); 
moveOffsets = moveOffsets(:); 

% determine the threshold (empirically)
trialsWithMove = resp==1|resp==2;% only left or right
respTime = respTime(trialsWithMove); 
intStartTime = intStartTime(trialsWithMove); 
resp = resp(trialsWithMove);

posAtIntStart = interp1(t, pos, intStartTime);
posAtResp = interp1(t, pos, respTime);

leftThresh = mean(posAtResp(resp==1)-posAtIntStart(resp==1));
rightThresh = mean(posAtResp(resp==2)-posAtIntStart(resp==2));

% now for each movement, determine whether it resulted in crossing that
% threshold within the first crossWindow (or within the duration of the move, if
% it was shorter than a crossWindow). 

Fs = 1000; rawT = t; rawPos = pos;
t = rawT(1):1/Fs:rawT(end);
pos = interp1(rawT, rawPos, t);

win = [0 crossWindow];
winSamps = win(1):1/Fs:win(2);
nM = numel(moveOnsets);
wheelSamps = bsxfun(@plus, moveOnsets(:), winSamps);
wheelMoves = interp1(t, pos, wheelSamps); % position trace triggered on move onset
wheelMoves = bsxfun(@minus, wheelMoves, wheelMoves(:,1)); % start of movement is zero

moveDurs = moveOffsets-moveOnsets;
moveType = zeros(size(moveOnsets)); % any that stay 0 are "flinches"
for m = 1:nM
    makesLeft = max(wheelMoves(m,winSamps<moveDurs(m)))>leftThresh;
    makesRight = min(wheelMoves(m,winSamps<moveDurs(m)))<rightThresh;
    if makesLeft && makesRight
        moveType(m) = 3; % some kind of equivocal thing, going both ways quickly
        % 3 is "unclassified" basically
    elseif makesLeft
        moveType(m) = 1;
    elseif makesRight
        moveType(m) = 2;
    end
end