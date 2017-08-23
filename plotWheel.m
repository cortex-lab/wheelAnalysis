

function h = plotWheel(t, pos, detectedMoves)
% function h = plotWheel(t, pos, detectedMoves)
%
% function to plot wheel traces (position and velocity) in a useful way. 
% That includes both position and velocity, on linked axes
% 
% t is timestamps, pos is positions, both are vectors of the same length. 
% detectedMoves is a struct which can have fields:
% - moveOnsets
% - moveOffsets
% - moveType (if classified)

% first compute an evenly-sampled position (in case that's not what was
% provided, if it is, that's ok this is fast) and velocity
rawT = t; rawPos = pos; 
Fs = 1000; 
t = rawT(1):1/Fs:rawT(end);
pos = interp1(rawT, rawPos, t);
vel = wheel.computeVelocity2(pos, 0.015, Fs);

if ~isempty(detectedMoves)
    if isfield(detectedMoves, 'moveOnsets'); moveOnsets = detectedMoves.moveOnsets(:); else; moveOnsets = []; end;
    if isfield(detectedMoves, 'moveOffsets'); moveOffsets = detectedMoves.moveOffsets(:); else; moveOffsets = []; end;
    if isfield(detectedMoves, 'moveType'); moveType = detectedMoves.moveType(:); else; moveType = []; end;
end

h = figure; 

ax1 = subplot(2,1,1); ax2 = subplot(2,1,2);
hold(ax1, 'on'); hold(ax2, 'on');

if ~isempty(moveOnsets)
    plot(ax1, moveOnsets, interp1(t, pos, moveOnsets), 'o', 'Color', [0 0.7 0]);
    plot(ax2, moveOnsets, interp1(t, vel, moveOnsets), 'o', 'Color', [0 0.7 0]);    
end
if ~isempty(moveOffsets)
    plot(ax1, moveOffsets, interp1(t, pos, moveOffsets), 'ro');
    plot(ax2, moveOffsets, interp1(t, vel, moveOffsets), 'ro');
end
if isempty(moveType) 
    % in this case just plot in and out of move as black and gray
    inMove = logical(WithinRanges(t, [moveOnsets moveOffsets]));
    plot(ax1,t(inMove), pos(inMove), '.', 'Color', [0 0.5 1]);
    plot(ax1,t(~inMove), pos(~inMove), '.', 'Color', 0.8*[1 1 1]);
    plot(ax2, t(inMove), vel(inMove), '.', 'Color', [0 0.5 1]);
    plot(ax2, t(~inMove), vel(~inMove), '.', 'Color', 0.8*[1 1 1]);
else
    colors(1,:) = [0 0.5 1];
    colors(2,:) = [1 0.5 0];
    colors(3,:) = 0.8*[1 1 1]; % unclassified, light gray
    colors(4,:) = [0 0 0]; % flinches, black
    moveType(moveType==0) = 4;
    
    % not in any move? gray
    inMove = logical(WithinRanges(t, [moveOnsets moveOffsets]));
    plot(ax1, t(~inMove), pos(~inMove), '.', 'Color', 0.8*[1 1 1]);
    plot(ax2, t(~inMove), vel(~inMove), '.', 'Color', 0.8*[1 1 1]);
    
    for r = 1:4
        inMove = logical(WithinRanges(t, [moveOnsets(moveType==r) moveOffsets(moveType==r)]));
        plot(ax1, t(inMove), pos(inMove), '.', 'Color', colors(r,:));
        plot(ax2, t(inMove), vel(inMove), '.', 'Color', colors(r,:));
    end    
end

linkaxes([ax1 ax2], 'x');
xlabel(ax2, 'time (sec)');
ylabel(ax2, 'velocity');
ylabel(ax1, 'position');


