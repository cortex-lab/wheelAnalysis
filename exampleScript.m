

% Script to demonstrate basic usage of wheel analysis package

%% pull the wheel data from a block

rawPos = block.inputSensorPositions;
rawTimes = block.inputSensorPositionTimes;

%% interpolate it to be regularly sampled

Fs = 1000;
t = rawTimes(1):1/Fs:rawTimes(end);
pos = interp1(times, rawPos, t);

wheelRadius = 5; % cm (burgess wheel)
wheelRadius = 15; % cm (running wheel)

rotaryEncoderResolution = 360; % number of ticks for one revolution

pos = pos/rotaryEncoderResolution*2*pi*wheelRadius; % convert to cm

%% find velocity and acceleration

[vel, acc] = wheel.computeVelocity(t, pos);

%% make the position signal be relative to a particular time(s)

posRel = wheel.resetAtEvent(t, pos, eventTimes);

%% compute movement onsets and offsets

thresh = 2; % cm/sec, a velocity threshold 
minBetweenMoves = 0.1; % sec
[moveTimes, moveAmplitudes, movePeakVelocities] = wheel.findAllMoves(...
    t, vel, pos, thresh, minBetweenMoves);


%% compute event triggered average and plot some traces


window = [-0.2 0.2];

[thisTr, thisStd, allT] = eventTrigAvgAllTraces(posRel, theseEventTimes-t(1), window, Fs);

timepnts = (1:numel(thisTr))/Fs+window(1);

plot(timepnts, allT, 'k');
hold on; 
plotWithErr(timepnts, thisTr, thisStd./sqrt(size(allT,1)), 'r');