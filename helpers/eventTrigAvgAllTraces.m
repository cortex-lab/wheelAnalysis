



function [trace, stdev, allTraces] = eventTrigAvgAllTraces(data, events, window, Fs)
% function [trace, stdev, allTraces] = eventTrigAvgAllTraces(data, events, window, Fs)
% 
% Returns the mean and standard deviation of "data" around each of the timestamps in "events". 
%
% This function differs from eventTrigAvg in that it will return you all
% snippets of the data around the events, rather than just the average and
% standard deviation. This makes it slower, and it can only work with one
% data trace at a time. 
%
% Inputs:
%   - data [1xN] Continuous data trace(s) (e.g. LFP, eye position), N
%   samples
%   - events [1xE] List of E timestamps, in seconds
%   - window [1x2] times, in seconds, relative to the events that you're
%   interested in. For example, to average data from 50ms prior to 200ms
%   after each event, use window = [-0.05 0.2]
%   - Fs - sampling frequency of data
%


if size(data,1)>size(data,2) 
    % make sure data is oriented correctly; has mild assumption that there
    % are more data points than channels.
    data = data';
end


numSamples = ceil((window(2)-window(1))*Fs)+1;

allDat = zeros(1,numSamples);

allTraces = zeros(length(events), numSamples);

% the algorithm here is to go through each event and just pull out that
% segment of the data, adding it to the running total and dividing, at the
% end, by the number of events. This algorithm turns out to actually be
% significantly faster, for most data sizes, than methods using bsxfun to
% avoid the for-loop, or the alternate for-loop method in which each
% element of allDat is computed sequentially. The pre-for-loop 
% parts take care of events whose window exeeds the bounds
% of the data. The "counts" variable has to have an entry to each different
% index of allDat since some might have more or less constituent data
% points due to these edge events. 
events = sort(events);
counts = zeros(1,numSamples);
frontEdgers = 0;
backEdgers = 0;

i = 1;
firstSample = ceil((events(i)+window(1))*Fs);
while firstSample<=0
    allDat((2-firstSample):end) = allDat((2-firstSample):end) + data(1:(ceil((events(i)+window(1))*Fs)+numSamples-1));
    allTraces(i,(2-firstSample):end) = data(1:(ceil((events(i)+window(1))*Fs)+numSamples-1));
    frontEdgers = frontEdgers+1;
    i = i+1;
    counts((2-firstSample):end) = counts((2-firstSample):end)+1;
    firstSample = ceil((events(i)+window(1))*Fs);
end

% now start at the end and work backwards to find all those events whose
% window exceeds the data in that direction
i = length(events); 
lastSample = ceil((events(i)+window(1))*Fs)+numSamples-1;
while lastSample>length(data)
    allDat(1:(numSamples-(lastSample-length(data)))) = allDat(1:(numSamples-(lastSample-length(data)))) + data(ceil((events(i)+window(1))*Fs):end);
    allTraces(i,1:(numSamples-(lastSample-length(data)))) = data(ceil((events(i)+window(1))*Fs):end);
    backEdgers = backEdgers+1;
    i = i-1;
    counts( 1:(numSamples-(lastSample-length(data)))) = counts( 1:(numSamples-(lastSample-length(data))))+1;
    lastSample = ceil((events(i)+window(1))*Fs)+numSamples-1;
end

for i = (1+frontEdgers):(length(events)-backEdgers)
    allDat = allDat + data(ceil((events(i)+window(1))*Fs):ceil((events(i)+window(1))*Fs)+numSamples-1);
    allTraces(i,:) = data(ceil((events(i)+window(1))*Fs):ceil((events(i)+window(1))*Fs)+numSamples-1);
end
counts = counts + length(events)-frontEdgers-backEdgers;

trace = allDat ./ counts;

if nargout>1 % requested standard deviation - go back through now that you have the mean and compute it
    
    % keep a running sum of the squared difference between observations
    % and the mean
    var = zeros(1, numSamples);
    counts = zeros(1,numSamples);
    frontEdgers = 0;
    backEdgers = 0;
    
    i = 1;
    firstSample = ceil((events(i)+window(1))*Fs);
    while firstSample<=0
        var((2-firstSample):end) = var((2-firstSample):end) + (trace((2-firstSample):end) - data(1:(ceil((events(i)+window(1))*Fs)+numSamples-1))).^2;
        frontEdgers = frontEdgers+1;
        i = i+1;
        counts((2-firstSample):end) = counts((2-firstSample):end)+1;
        firstSample = ceil((events(i)+window(1))*Fs);
    end
    
    % now start at the end and work backwards to find all those events whose
    % window exceeds the data in that direction
    i = length(events);
    lastSample = ceil((events(i)+window(1))*Fs)+numSamples-1;
    while lastSample>length(data)
        var(1:(numSamples-(lastSample-length(data)))) = var(1:(numSamples-(lastSample-length(data)))) + (trace(1:(numSamples-(lastSample-length(data)))) - data(ceil((events(i)+window(1))*Fs):end)).^2;
        backEdgers = backEdgers+1;
        i = i-1;
        counts( 1:(numSamples-(lastSample-length(data)))) = counts( 1:(numSamples-(lastSample-length(data))))+1;
        lastSample = ceil((events(i)+window(1))*Fs)+numSamples-1;
    end
    
    for i = (1+frontEdgers):(length(events)-backEdgers)
        firstSample = ceil((events(i)+window(1))*Fs);
        lastSample = ceil((events(i)+window(1))*Fs)+numSamples-1;
        var = var+(trace - data(firstSample:lastSample)).^2;
    end
    
    counts = counts + length(events)-frontEdgers-backEdgers;
    for c = 1:size(var,1)
        var(c,:) = var(c,:) ./ (counts-1);
    end
    stdev = sqrt(var);    
    
else
    stdev= [];
end

