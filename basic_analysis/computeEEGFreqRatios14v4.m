function [EEGFreqPeakRatios, EEGFreqSumRatios, EEGFreqCentroid] = computeEEGFreqRatios14v4(EEGpsd, dFs, settings) 
% TPG 10-24-2007
% Fs is the sampling frequency of the preprocessed/predecimated data in Hz.
% The frequency bands to compute are in the "settings" struct.  For example,
% settings.EEG_Ratio_Bands might look like this to compute ratio of the peak in the 0-2Hz
% band and the 2-5Hz band: settings.EEG_Ratio_Bands = [0 2; 2 5]';

EEGLowPeak = max(EEGpsd(floor(settings.EEG_Ratio_Bands(1,1)/(dFs/2)*((settings.LFPNfft/2)-1))+2 : floor(settings.EEG_Ratio_Bands(2,1)/(dFs/2)*((settings.LFPNfft/2)-1))+2));
EEGHighPeak = max(EEGpsd(floor(settings.EEG_Ratio_Bands(1,2)/(dFs/2)*((settings.LFPNfft/2)-1))+2 : floor(settings.EEG_Ratio_Bands(2,2)/(dFs/2)*((settings.LFPNfft/2)-1))+2));
EEGFreqPeakRatios = EEGHighPeak / EEGLowPeak;

EEGLowSum = sum(EEGpsd(floor(settings.EEG_Ratio_Bands(1,1)/(dFs/2)*((settings.LFPNfft/2)-1))+2 : floor(settings.EEG_Ratio_Bands(2,1)/(dFs/2)*((settings.LFPNfft/2)-1))+2));
EEGHighSum = sum(EEGpsd(floor(settings.EEG_Ratio_Bands(1,2)/(dFs/2)*((settings.LFPNfft/2)-1))+2 : floor(settings.EEG_Ratio_Bands(2,2)/(dFs/2)*((settings.LFPNfft/2)-1))+2));
EEGFreqSumRatios = EEGLowSum / EEGHighSum;

EEGFreqCentroid = centroid(EEGpsd(floor(settings.EEG_Ratio_Bands(1,1)/(dFs/2)*((settings.LFPNfft/2)-1))+2 : floor(settings.EEG_Ratio_Bands(2,2)/(dFs/2)*((settings.LFPNfft/2)-1))+2)');
