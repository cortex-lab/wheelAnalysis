


function [moveTimes, moveAmp, movePeakVel] = findAllMoves(times, vel, pos, thresh, minBetweenMoves)
% function [moveTimes, moveAmp, movePeakVel] = findAllMoves(times, vel, pos, thresh, minBetweenMoves)
% Here, we look for a threshold crossing, then find the subsequent time
% that the velocity is less than the threshold for some amount of time
% (minBetweenMoves). That's the end of the movement, and the total
% amplitude is computed as the difference in position between the end and
% the start. 
%
% - moveTimes will be 2 x N for N detected movements, a start and end time
% - moveAmp and movePeakVel are 1xN, the total signed amplitude of the
% movement and the peak velocity attained (these may have opposite sign!)

timesFs = mean(diff(times)); % assuming times/vel/pos are evenly sampled in time! Use interp1 to do this if it isn't true.

velAboveThresh = vel>thresh | vel<-thresh;

possibleMoveStarts = find( (~velAboveThresh(1:end-1) & velAboveThresh(2:end)) );


betweenMovesSamps = ceil(minBetweenMoves/timesFs);

% this will be length(vel), each entry true if it is not near any
% movement according to minBetweenMoves
isEndCondition = conv(double(velAboveThresh), ones(1, betweenMovesSamps), 'same')==0;
                
moveTimes = zeros(2, length(possibleMoveStarts));                
moveAmp = zeros(1, length(possibleMoveStarts));    
movePeakVel = zeros(1, length(possibleMoveStarts));  
                
pMoveStartInd = 1;
moveInd = 1;
lastEndTime = 0;
% try
while pMoveStartInd<length(possibleMoveStarts)
    if times(possibleMoveStarts(pMoveStartInd))<lastEndTime
        % this possibleMoveStart happened during the last movement; just
        % move on.
        pMoveStartInd = pMoveStartInd+1;
    else
        
        startSamp = possibleMoveStarts(pMoveStartInd);
        endSamp = find(isEndCondition(possibleMoveStarts(pMoveStartInd)+1:end),1) + startSamp;
        
        if isempty(endSamp) % this can happen when the movement is going on at the time the recording ends
            endSamp = length(isEndCondition);
        end
        
        moveTimes(1, moveInd) = times(startSamp);
        moveTimes(2, moveInd) = times(endSamp);
        
        moveAmp(moveInd) = pos(endSamp)-pos(startSamp);
        thisVel = vel(startSamp:endSamp);
        movePeakVel(moveInd) = thisVel(find(abs(thisVel)==max(abs(thisVel)),1));
        
        lastEndTime = moveTimes(2, moveInd);
        moveInd = moveInd+1;
        pMoveStartInd = pMoveStartInd+1;
        
    end
    
end
% catch me
%     keyboard
% end
    
% pick out just the entries of moveTimes we had actually filled
moveTimes = moveTimes(:,1:moveInd-1);
moveAmp = moveAmp(:,1:moveInd-1);
movePeakVel = movePeakVel(:,1:moveInd-1);