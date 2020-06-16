

function posOut = correctCounterDiscont(pos)
% rotary encoder sometimes gets weird values. This fixes them. Algorithm
% from CB. 

posDiff = diff(pos(:)); % compute velocity
% correct for counter flow discontinuities
posDiff(posDiff > 2^31) = posDiff(posDiff > 2^31) - 2^32;
posDiff(posDiff < -2^31) = posDiff(posDiff < -2^31) + 2^32;
posOut = cumsum([0; posDiff]); % back to position
posOut = reshape(posOut, size(pos))+pos(1);