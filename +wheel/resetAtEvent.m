

function posRel = resetAtEvent(t, pos, eventTimes)
% function resetAtEvent(t, pos, eventTimes)
% - t is the time points at which wheel position was measured
% - pos are the position measurements
% - eventTimes are the times at which pos should reset to zero

if ~isrow(t)
    t = t';
end
if ~isrow(eventTimes)
    eventTimes = eventTimes';
end 

% find the next t after each eventTimes
[~,ii] = sort([t eventTimes]);
eventTs = ii(length(t)+1:end)+1;

posAtEvents = zeros(size(pos));
posAtEvents(eventTs) = pos(eventTs);

posRel = pos-cumsum(posAtEvents);
