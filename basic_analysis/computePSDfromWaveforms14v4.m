function [percentOfBandPowers, totalPower, temp2, tempf2, normPSD] = computePSDfromWaveforms14v4(decimatedData, dFs, settings) 
% TPG 10-24-2007
% Fs is the sampling frequency of the preprocessed/predecimated data in Hz.
% The frequency bands to compute are in the "settings" struct.  For example,
% settings.LFPFrequencyBands might look like this to compute the percent power in the
% 5-10Hz and 15-20Hz bands:  settings.LFPFrequencyBands = [5 10; 15 20;]';

[temp2,tempf2]=pwelch(decimatedData,settings.LFPNfft,settings.LFPNoverlap,[],dFs); 

numBands = size(settings.LFPFrequencyBands,2);
powers = -1*ones(1,numBands);
for n = 1:numBands
    powers(n) = sum(temp2(floor(settings.LFPFrequencyBands(1,n)/(dFs/2)*((settings.LFPNfft/2)-1))+2: ...
                            floor(settings.LFPFrequencyBands(2,n)/(dFs/2)*((settings.LFPNfft/2)-1))+2));  % spectral power in this frequency band
end
totalPower = sum(temp2);
if (settings.F_BaseNormalizingLimit > 61) || (settings.F_LowerNormalizingLimit < 59) % if the requested normalization band does NOT include the 60 Hz line noise spike
    normalizingPower = sum(temp2(floor(settings.F_BaseNormalizingLimit/(dFs/2)*((settings.LFPNfft/2)-1))+2: ...
                                floor(settings.F_LowerNormalizingLimit/(dFs/2)*((settings.LFPNfft/2)-1))+2)); % normalize only by the power in (for example) the 0-50Hz band, to avoid harmonic peaks of 60Hz (like 60Hz or 180Hz) from throwing off the normalization
else  % if the requested normalization band DOES include the 60 Hz line noise spike, then subtract that power out of the normalization power
    normalizingPower = sum(temp2(floor(settings.F_BaseNormalizingLimit/(dFs/2)*((settings.LFPNfft/2)-1))+2: ...
                                floor(settings.F_LowerNormalizingLimit/(dFs/2)*((settings.LFPNfft/2)-1))+2)) ...
                     - sum(temp2(floor(59/(dFs/2)*((settings.LFPNfft/2)-1))+2: ...
                                floor(61/(dFs/2)*((settings.LFPNfft/2)-1))+2)); 
end
normPSD = temp2 / normalizingPower;
%percentPowers = powers * 100 / totalPower;
percentOfBandPowers = powers * 100 / normalizingPower;

if settings.showFreqGraphs % for debugging
    figure; plot(1/dFs:1/dFs:1/dFs*length(decimatedData),decimatedData); title(sprintf('Waveform, DecimatedFs=%g',dFs));
    figure; plot(tempf2, temp2); title(sprintf('PSD, Nfft=%g,Noverlap=%g',settings.LFPNfft,settings.LFPNoverlap));
    keyboard;
end
