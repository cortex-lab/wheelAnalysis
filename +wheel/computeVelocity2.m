function [vel, acc] = computeVelocity2(pos, smoothSize, Fs)
% function [vel, acc] = computeVelocity(pos, smoothSize, Fs)
%
% assumes pos is uniformly sampled in time
%
% smooth size is in units of seconds

% area of this smoothing window is 1 so total values are unchanged - units
% don't change
% smoothWin = wheel.gausswin(smoothSize)./sum(gausswin(smoothSize));
smoothWin = myGaussWin(smoothSize, Fs); 
pos = pos(:);
vel = [0; conv(diff(pos), smoothWin, 'same')]*Fs; % multiply by Fs to get cm/sec

if nargout>1
    % here we choose to apply the smoothing again - it's sort of
    % empirically necessary since derivatives amplify noise. 
    acc = [0; conv(diff(vel), smoothWin, 'same')]*Fs; %cm/sec^2
end