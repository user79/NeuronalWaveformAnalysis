function [filteredData, dFs, clipping] = preprocessWaveforms14v4(data, Fs, settings) 
clipping = 0;

if (max(data) == 5) || (min(data) == -5)
    disp('ADC saturation detected!');
    clipping = 1;
    % todo - eliminate this clipping somehow?
end

if (round(Fs) == 510) % Andrew's recordings
    savedDF = settings.DecimationFactor;
    settings.DecimationFactor = 2.55102040816326530000;  % to reduce the dFs to ~200 samples/sec
end
    
% decimatedData = decimate(data,settings.DecimationFactor); %decimate() also does its own automatic antialiasing filtering
[p,q] = rat(settings.DecimationFactor);
decimatedData = resample(data, q, p); %resample() does its own automatic antialiasing filtering

filteredData = filter(settings.HdpreprocessWaveforms, decimatedData);% Lowpass filter the LFP data
% filteredData = data;

if settings.NormalizeLFP 
    % filteredData = filteredData / max(abs(filteredData));   % divide by max of absolute value for normalization -- after this step, the new absolute maximum of every lfp recording is 1
    filteredData = filteredData / sqrt(sum(filteredData.*conj(filteredData))/size(filteredData,1));   % divide by the RMS value for normalization
end

filteredData = detrend(filteredData); % remove mean before FFT

dFs = Fs/settings.DecimationFactor;

if (round(Fs) == 510) % Andrew's recordings
    settings.DecimationFactor = savedDF;  % restore the decimation factor to what it was before...
end

