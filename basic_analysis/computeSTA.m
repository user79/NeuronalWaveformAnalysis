function [STA, randomizedSTA, peakToTroughDelta, peakToTroughRandomizedDelta, PeakToTroughNormalToRandomizedRatio] = computeSTA(spikeTimestamps, rawWaveData, Fs, hFilter)
% computes the spike-triggered-average of the LFP or other waveform signal, as in Avila 2010
% TG 5/24/2010
% input parameters in computeSettings:

  STAwindowWidth = 3; % seconds, including before and after spike

% %pseudocode:
% - filter the waveform (e.g. from 12 to 25 Hz or 25 to 40 Hz, cf. Avila et al 2010)
% - compute the STA
% - compute the randomized STA
%   - convert the timestamps to ISIs
%   - randomize the ISIs
%   - convert ISIs back to timestamps
%   - compute the STA

if exist('hFilter','var')
    waveData = filter(hFilter, rawWaveData);% bandpass filter the LFP data before taking the STA
else
    waveData = rawWaveData;
end

STAwindowWidthIndex = Fs*STAwindowWidth + 1;  % this needs to be an odd number
tempArray = zeros(size(spikeTimestamps,1),STAwindowWidthIndex);
for n=1:size(spikeTimestamps,1)
    if (floor(spikeTimestamps(n)*Fs-(STAwindowWidthIndex-1)/2) >= 1) && (floor(spikeTimestamps(n)*Fs+(STAwindowWidthIndex-1)/2) <= size(waveData,1)) % check for out-of-bounds
        tempArray(n,:) = waveData(floor(spikeTimestamps(n)*Fs)-(STAwindowWidthIndex-1)/2:floor(spikeTimestamps(n)*Fs)+(STAwindowWidthIndex-1)/2); % add the LFP data surrounding this spike to the temp array
    end
end
STA = mean(tempArray); % take the average to get the STA
peakToTroughDelta = max(STA)-min(STA);

isi = [diff(spikeTimestamps) rand(length(spikeTimestamps)-1,1)]; % shuffle the ISI's, then convert back to timestamps
randISI = sortrows(isi,2); 
randISI = randISI(:,1);
randomizedSpikeTimestamps = cumsum(randISI);
tempArray = zeros(size(spikeTimestamps,1),STAwindowWidthIndex);
for n=1:size(randomizedSpikeTimestamps,1)
    if (floor(randomizedSpikeTimestamps(n)*Fs-(STAwindowWidthIndex-1)/2) >= 1) && (floor(randomizedSpikeTimestamps(n)*Fs+(STAwindowWidthIndex-1)/2) <= size(waveData,1))
        tempArray(n,:) = waveData(floor(randomizedSpikeTimestamps(n)*Fs)-(STAwindowWidthIndex-1)/2:floor(randomizedSpikeTimestamps(n)*Fs)+(STAwindowWidthIndex-1)/2);
    end
end
randomizedSTA = mean(tempArray);
peakToTroughRandomizedDelta = max(randomizedSTA)-min(randomizedSTA);
PeakToTroughNormalToRandomizedRatio = peakToTroughDelta / peakToTroughRandomizedDelta;

STA = STA(:); % put into column vector form
randomizedSTA = randomizedSTA(:);
