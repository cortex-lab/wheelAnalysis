

% Script to demonstrate basic usage of wheel analysis package

%% pull the wheel data from a block

rawPos = block.inputSensorPositions;
rawPos = wheel.correctCounterDiscont(rawPos); % correction because sometimes negative wheel positions wrap around
rawTimes = block.inputSensorPositionTimes;

%% interpolate it to be regularly sampled

Fs = 1000;
t = rawTimes(1):1/Fs:rawTimes(end);
pos = interp1(rawTimes, rawPos, t, 'linear');

wheelRadius = 31; % mm (burgess wheel) (measured by Chris)
wheelRadius = 150; % mm (running wheel) (guess!)

rotaryEncoderResolution = 360*4; % number of ticks for one revolution (factor of 4 is according to CB)

pos = pos./(rotaryEncoderResolution)*2*pi*wheelRadius; % convert to mm

nSampToPlot = min(50000, length(pos));

figure;
subplot(3,1,1);
plot(t(1:nSampToPlot), pos(1:nSampToPlot));
xlim([t(1) t(nSampToPlot)]);
xlabel('time (sec)'); 
ylabel('wheel position (mm)');

%% find velocity and acceleration

[vel, acc] = wheel.computeVelocity(pos, 50, Fs);

subplot(3,1,2);
plot(t(1:nSampToPlot), vel(1:nSampToPlot));
xlim([t(1) t(nSampToPlot)]);
xlabel('time (sec)'); 
ylabel('wheel velocity (mm/sec)');

subplot(3,1,3);
plot(t(1:nSampToPlot), acc(1:nSampToPlot));
xlim([t(1) t(nSampToPlot)]);
xlabel('time (sec)'); 
ylabel('wheel acceleration (mm/sec/sec)');

%% make the position signal be relative to a particular time(s)
tr = [block.trial];
tr = tr(1:block.numCompletedTrials);
eventTimes = [tr.interactiveStartedTime];

posRel = wheel.resetAtEvent(t, pos, eventTimes);

%% compute movement onsets and offsets

thresh = 30; % cm/sec, a velocity threshold 
minBetweenMoves = 0.2; % sec
[moveTimes, moveAmplitudes, movePeakVelocities] = wheel.findAllMoves(...
    t, vel, pos, thresh, minBetweenMoves);


figure; 
plot(t(1:nSampToPlot), pos(1:nSampToPlot));
hold on;
xlim([t(1) t(nSampToPlot)]);
xlabel('time (sec)'); 
ylabel('wheel position (mm)');
yl = ylim();
[xx,yy] = rasterize(moveTimes(1,moveTimes(1,:)<t(nSampToPlot))); 
plot(xx,yy*diff(yl)+yl(1));

%% compute event triggered average and plot some traces

interactiveStart = [tr.interactiveStartedTime];
responseMade = [tr.responseMadeTime];
responseMadeStartTime = nan(size(responseMade));
for r = 1:length(responseMade)
    theseMoves = find(moveTimes(1,:)>interactiveStart(r) & moveTimes(1,:)<responseMade(r));
    if ~isempty(theseMoves)
        responseMadeStartTime(r) = moveTimes(1,theseMoves(end));
    end
end
window = [0 0.7];

[thisTr, thisStd, allT] = eventTrigAvgAllTraces(posRel, responseMadeStartTime(~isnan(responseMadeStartTime))-t(1), window, Fs);

timepnts = (1:numel(thisTr))/Fs+window(1);

figure;
plot(timepnts, allT, 'k');
hold on; 
plotWithErr(timepnts, thisTr, thisStd./sqrt(size(allT,1)), 'r');