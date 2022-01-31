function [computeSettings, graphSettings] = defaultMultipleRunSettings14v4()
% TPG 10-03-07
% TODO: Move ALL settings into this structure

computeSettings = struct(...
    ...
    ... % ***** General Settings ********************************************************
    ...
    'version', {14.4}, ... % June 7 2010
    'NumberOfMeasures', {31}, ... % this number is only used to form the "measures" array... other measures may be returned in structs
    'UnixPathDependencies', {''}, ... % the directory in which the Matlab main codebase resides, for the worker nodes
    'filelistPathDependency', {''}, ... % the directory containing the lists of files to be processed, for the initial file (before distribution to the workers)
    'baseDataDirectoryPrefix', {''}, ...  % the directory in which the data resides
    'outputDirectoryPrefix', {''}, ... % the directory in which to write or look for the output
    'processBasic', {1}, ... %(FR+CV)
    'processSE', {1}, ...  % Sample Entropy
    'processWaveformPSD', {1}, ...   %spectra on raw LFP waveform
    'processSpikePSD', {1}, ...   %Spike Oscillations
    'processBursts', {1}, ...  %Burst analysis
    'processISIHistsAndACs', {1}, ... %ISIhistograms and Autocorrelograms
    'processAdaptiveFilters', {0}, ... % LMS, Kalman, etc
    ...
    ... % ***** Debug settings **********************************************************
    ...
    'verbose',{0}, ... % set this to be 1 or higher for more debugging info
    'showFreqGraphs',{0}, ... % produce graphs of EVERY waveform and periodogram! for debugging
    'showPoissonGraphs',{0}, ... % produce graphs of EVERY poisson comparison
    'interactiveACmode',{0}, ...  % set to 1 if manually checking all AC and class 1,2, or 3 classifications,
    ...    % otherwise set to 0 to have these done "automatically" (lower quality though) and saved
    'saveACinfoFiles',{0}, ... % save into the "output" directory for later checking
    'debugForcedEndblock', {-1}, ... % if "-1", no effect, otherwise, block (eg. block 10) to stop analyzing each record, for debugging purposes (to speed things up) - cf. the SON toolbox documentation
    ...
    ... % ***** Basic Spike processing and Burst processing *****************************
    ...
    'maxISI',{5}, ... %5 seconds, but TODO: fix this somehow
    'processOnlyFirstSortedSpike',{1}, ...  % if 1, this only uses Spike2 spike code "01", i.e. the
    ...    % first sorted neuron, regardless of how many sorted neurons there are
    ...    % in the file.  If 0, this code will eventually process each spike separately.     
    'CVpmRefractoryPeriod', {0.001}, ... % assumed refractory period (in seconds) for computing the CVpm, the "CV percent of maximum" defined by Chelvanayagam and Vidyasagar JNeurosci 2006
    'spiketrainLengthLimit', {0},... {2000}, ... % if this is zero, will not impose any limit.  If > 0, will crop the spiketrain to this number of spikes.  If the spiketrain is shorter than this number, it will be used as is and a warning will be printed.
    'spiketrainTimeLimit', {0}, ... % % if this is zero, will not impose any limit.  If > 0, will crop the spiketrain to this time limit. If the the last spike in the spiketrain is sooner than this duration (in seconds), the spiketrain will be used as is and a warning will be printed.
    'maxSpikesExpectedPerBin', {9}, ... % for calculating the discharge density histogram - % no spikes expected closer together than this rate (per bin, of width 1/meanFiringRate)
    'PoissonComparisonEMIN',{5}, ... % if zero, uses all values from maxSpikesExpectedPerBin.  If positive, combines bins that have counts less than this number.
    'PoissonComparisonBlockSize', {0}, ... % or {100}... if zero, uses whole spiketrain, otherwise, splits into blocks and computes chi-square comparison separately on each block
    'PoissonMean', {1}, ...
    'ClassLabels', {{'''Regular''', '''Irregular''', '''Bursty''', '''Unclassified'''}}, ... % Brigitte's three types, type 1, 2, and 3
    'OscillatoryLabels', {{'Oscillatory','Non-Oscillatory' }}, ... % % Soares et al. JNeurosci. 2004
    'PoissonComparisonLabels', {{'''Regular''', '''Poisson/random''', '''Bursty'''}}, ... % re Levy paper- three types, either poisson-random, or more regular or more bursty
    'BurstDetectionMethod', {'WichmannPoissonSurprise'}, ... % or, 'HanesPoissonSurprise' or 'RankSurprise'
    'burstLimit', {0.05}, ...   %0.020;  %TODO: verify this number (chosen rather arbitrarily!)
    ...    % From the paper:
    ...    % A nonparametric approach for detection of bursts in spike trains
    ...    % Boris Gourévitch, and Jos J. Eggermont
    ...    % Journal of Neuroscience Methods
    ...    % Volume 160, Issue 2, 15 March 2007, Pages 349-358
    'RSalpha', {0.5}, ...  %-log(0.01); % cf. also the paper above by Gourevitch/Eggermont
    'burstFactor', {2}, ... % see Tom Wichmann's "legendy_new3.m" for details on the burst detection parameters
    'burstSamplingRate', {1}, ... % previously I tried 25000, but the way Tom's code is written may be slightly misleading on the name of this param
    'minBurstLength', {4}, ...
    'burstLocalLength', {0}, ...
    'burstSurpriseCutoff', {5}, ...
    'arOrderArray', {[250 250 1]}, ... % see ARX() documentation - first # is # of AR ('a') coeffs, 2nd # is the input ('b') coeffs, 3rd is delay
    ...
    ... % ***** Frequency Domain / Power Spectral Density parameters ********************
    ...
    'DecimationFactor', {5}, ... reduces 1000 samples/sec to, for example, 200 samples/sec if DecimationFactor = 5
    ...% %     'Fpass', {90}, ...     % preprocessing-filter Passband Frequency
    ...% %     'Fstop', {100}, ...    % preprocessing-filter Stopband Frequency
    ...% %     'Apass', {1},   ...    % preprocessing-filter Passband Ripple (dB)
    ...% %     'Astop', {80}, ...     % Stopband Attenuation (dB)
    'Fstop1', {0.1}, ...     % preprocessing-filter stop Frequency1
    'Fpass1', {1}, ...     % preprocessing-filter Passband Frequency1
    'Fpass2', {99}, ...    % preprocessing-filter Passband Frequency2
    'Fstop2', {100}, ...    % preprocessing-filter Stopband Frequency2
    'Apass', {1},   ...    % preprocessing-filter Passband Ripple (dB)
    'Astop', {80}, ...     % Stopband Attenuation (dB)
    'sta12t25Fstop1', {10}, ...     % 12-25Hz filter for STA processing, stop Frequency1
    'sta12t25Fpass1', {12}, ...     % 12-25Hz filter for STA processing, preprocessing-filter Passband Frequency1
    'sta12t25Fpass2', {25}, ...    % 12-25Hz filter for STA processing, preprocessing-filter Passband Frequency2
    'sta12t25Fstop2', {27}, ...    % 12-25Hz filter for STA processing, preprocessing-filter Stopband Frequency2
    'sta25t40Fstop1', {23}, ...     % 25-40Hz filter for STA processing, stop Frequency1
    'sta25t40Fpass1', {25}, ...     % 25-40Hz filter for STA processing, preprocessing-filter Passband Frequency1
    'sta25t40Fpass2', {40}, ...    % 25-40Hz filter for STA processing, preprocessing-filter Passband Frequency2
    'sta25t40Fstop2', {42}, ...    % 25-40Hz filter for STA processing, preprocessing-filter Stopband Frequency2
    'NormalizeLFP', {1}, ... normalize the LFP amplitude before filtering, downsampling, and Fourier transforming (i.e. by the max, or by the rms value), default=1 (yes)
    'LFPNfft', {256}, ... length of FFT to use in Welch's periodogram PSD 
    'LFPNoverlap', {64}, ... number of points to overlap in Welch's PSD
    'LFPFrequencyBands',{[5 10; 10 15; 15 20; 20 25; 25 30; 30 35; 35 40; 40 45; 45 50; 50 55; 55 60; 60 65; 65 70; 70 75; 75 80; 80 85; 85 90; 90 95; 95 100;]'}, ... % the first row contains the starting frequencies, the second row contains the ending frequencies
    'LFPFrequencyBandLabels',{{'5-10Hz','10-15Hz','15-20Hz','20-25Hz','25-30Hz','30-35Hz','35-40Hz','40-45Hz','45-50Hz','50-55Hz','55-60Hz','60-65Hz','65-70Hz','70-75Hz','75-80Hz','80-85Hz','85-90Hz','90-95Hz','95-100Hz'}},...
    'saveFreqGraphsDuringComputation',{1},... % this is a kludge... normally graphs would be formed and saved as a second step
    'F_upperGraphLimit',{61}, ... % (in Hz) top of spectral graphs
    'F_lowerGraphLimit',{0}, ... % (in Hz) bottom of spectral graphs
    'lineLocation',{58}, ... % y-height (in Hz) at which to draw the black and white bar showing which rat is which 
...%     'F_BaseNormalizingLimit', {5}, ... from this frequency (in Hz) to the F_LowerNormalizingLimit is the first band to take the normalizing power from, to avoid the 60Hz peak
...%     'F_LowerNormalizingLimit', {55}, ...
    'F_BaseNormalizingLimit', {5}, ... from this frequency (in Hz) to the F_LowerNormalizingLimit is the first band to take the normalizing power from, to avoid the 60Hz peak ...  if it surrounds the 60 Hz band, the computePSDsFromWaveforms script should still eliminate the 60 Hz peak. 
    'F_LowerNormalizingLimit', {95}, ... % try 65-80Hz to address LD paper reviewers' comments about normalizing band
    'EEG_Ratio_Bands', {[0 2; 2 5]'}, ...
    ...
    ... % ***** Sample Entropy  *********************************************************
    ...
    'SEm',{2}, ... % SampEntropy "m" value, the embedding dimension or maximum run length to check...
    ...    % if [], default is 5-1 = 4 returned values... see sampen.m
    ...    % cf. Comparison of the Use of Approximate Entropy and Sample Entropy:
    ...    % Applications to Neural Respiratory Signal, 1Xinnian Chen, 1,2Irene C.
    ...    % Solomon and 1,2Ki H. Chon
    ...    % http://ieeexplore.ieee.org/iel5/10755/33900/01615393.pdf?arnumber=1615393
    ...    % and Richman JS, Moorman JR.  Physiological time-series analysis using
    ...    % approximate entropy and sample entropy. Am J Physiol Heart Circ Physiol.
    ...    % 2000 Jun;278(6):H2039-49.
    ...    % http://ajpheart.physiology.org/cgi/content/full/278/6/H2039
    'SampleEntropyLabels',{{'SEm=2','SEm=3','SEm=4','SEm=5'}}, ...
    'SEr', {0.2}, ... % tolerance, normalized to standard deviation
    'SampEnStandardLength',{-1},... %400; % minimum number of ISI's to use file, and number to
    ...    % truncate ISI sequence to. Recommended (see above) min length somewhere
    ...    % between 10^m and 20^m, chose 20^2=400.  If this is "-1", uses all spikes
    ...    % and does not check for length.    
    ...
    ... % ***** ISI and Autocorrelation parameters **************************************
    ...
    'isiHistogramNumOfBins', {1000}, ...
    'isiHistogramTime', {0.01}, ...
    'countACpeaksTime', {0.7}, ...
    'forceSpikewavePSD',{1} ...  % use this flag to force computations of the low-frequency spectral
    ...    % power on the fast-sampled, highpass-filtered spike waveform data, for
    ...    % example if a comparison is needed with Brigitte's or Milind's old data,
    ...    % which did NOT record the frequency bands of LFP (i.e. those frequencies
    ...    % were filtered out by the highpass filter).
);


% String parameters to load / important for graphing:
strparams = {'GroupName'};

% Single-digit individual-record measures to compute and graph: (these names become the graph labels)
sdparams = {
    'FiringRate'
    'CoefficientOfVariation'
    'CVpm'
    'BurstIndex'
%     'BurstRate'
%     'SpkPSD_3to8Hz'
%     'SpkPSD_8to15Hz'
%     'SpkPSD_gt15Hz'
%     'PoissonComparison'
%     'PoissonComparisonMean'
%     'MeanIntraburstFrequency'
    'ProportionSpikesInBursts'
%     'ProportionTimeInBursts'
    'rangeCounts'
%     'L'
    'SampleEntropy'    
%     'LFPCentroids'
%     'PeakToTroughNormalToRandomizedRatio12t25'
%     'PeakToTroughNormalToRandomizedRatio25t40'
%
%     'PtTD'
%     'PtTDRRD'
%     'PtTD12t25'
%     'PtTDRRD12t25'
%     'PtTD25t40'
%     'PtTDRRD25t40'
%     'PtTDEEG1'
%     'PtTRRDEEG1'
%     'PtTDEEG2'
%     'PtTRRDEEG2'
%     'PtTDEEGdiff'
%     'PtTRRDEEGdiff'
% 
%     'MeanSurprise'          % all burst-related statistics currently will give an error in createGraphs14v3 if none of a particular group's neurons have any bursts... TODO:  fix this...
%     'NumBursts'
%     'MeanSpikesPerBurst'
%     'MedianSpikesPerBurst'
%
%     'EEGdiffLFPcorrZL'
% %     'EEG1LFPcorrZL'
% %     'EEG2LFPcorrZL'
% 
%       'arBestMSE'
%       'arUndecimatedMSE'
%       'bestScale'
% %       'sA'
% %       'sB'
%       'mA'
%       'mB'
% %       'rsAB'
%       'rmAB'
% %       'ARstabilitySum'
%       'lengthA'
%       'lengthB'
%
%     'lmsEndMSEeeg1'
%     'lmsEndToBeginningMSEeeg1'
%     'lmsNormMSEeeg1'
%     'lmsEndMSEeeg2'
%     'lmsEndToBeginningMSEeeg2'
%     'lmsNormMSEeeg2'
%     'lmsEndMSEeegdiff'
%     'lmsEndToBeginningMSEeegdiff'
%     'lmsNormMSEeegdiff'
%     'kalmanEndMSEeeg1'
%     'kalmanEndToBeginningMSEeeg1'
%     'kalmanNormMSEeeg1'
%     'kalmanEndMSEineeg1'
%     'kalmanEndToBeginningMSEineeg1'
%     'kalmanNormMSEineeg1'
%     'kalmanEndMSEeeg2'
%     'kalmanEndToBeginningMSEeeg2'
%     'kalmanNormMSEeeg2'
%     'kalmanEndMSEineeg2'
%     'kalmanEndToBeginningMSEineeg2'
%     'kalmanNormMSEineeg2'
%     'kalmanEndMSEeegdiff'
%     'kalmanEndToBeginningMSEeegdiff'
%     'kalmanNormMSEeegdiff'
%     'kalmanEndMSEineegdiff'
%     'kalmanEndToBeginningMSEineegdiff'
%     'kalmanNormMSEineegdiff'
%     
%     'EEG1FreqPeakRatios'
%     'EEG1FreqSumRatios'
%     'EEG1FreqCentroid'
%     
    'NumSpikes'
    'Duration'
%     'MedianSurprise'
%     'MedianIntraburstFrequency'
    };

% Array parameters to load and graph:
arrparams = {
    'LFPPowers'
    'NormalizedLFPPowers'
    'EEGDiffLFPCOHs'
    'EEG1LFPCOHs'
    'EEG2LFPCOHs'
% %     'SampleEntropy'
    };

% Array parameter labels to load and graph, in the same order as arrparams -
%    These values are the names of a field in the computeSettings struct that contains
%    the actual string labels!
arrparamlabels = {
    'LFPFrequencyBandLabels'
    'LFPFrequencyBandLabels'
    'LFPFrequencyBandLabels'
    'LFPFrequencyBandLabels'
    'LFPFrequencyBandLabels'
%     'SampleEntropyLabels'
    };

% waveform parameters to load and graph:
waveformparams = {
    'LFPPSD'
    'normalizedLFPPSD'
    'EEG1PSD'
    'EEG2PSD'
    'EEGdiffPSD'
    'EEG1LFPCOH'
    'EEG2LFPCOH'
    'EEGdiffLFPCOH'
    'STA'   
    'EEG1STA'   
    'EEG2STA'   
    'EEGdiffSTA'   
    'randomizedSTA'
    'STA12t25'   % todo: change the plot co-argument instead of using the PSD frequencies
    'randomizedSTA12t25'
    'STA25t40'
    'randomizedSTA25t40'
    };

% Parameters that are graphed as percentages of different types rather than as individual
% mean/SE values:  (e.g. "10% of the contralateral normal STN group of recordings had Class=bursty, 60% had
% Class=regular, and 30% had Class=irregular")
groupparams = {
%     'Class'
    'PoissonComparison'
    'Oscillatory'
};
% Array parameter labels to load and graph, in the same order as arrparams -
%    These values are the names of a field in the computeSettings struct that contains
%    the actual string labels!
groupparamlabels = {
%     'ClassLabels'
    'PoissonComparisonLabels'
    'OscillatoryLabels'
    };

% Other parameters that exist, but don't need to be graphed currently:
otherparams = {
    'ErrorString'
    'Date'
    'Filename'
    'CountAC'
    'PoissonComparisonH'
    };


graphSettings = struct( ...
    ...
    ... % ***** Graphing parameters *****************************************************
    ...
    'recordSingleDigitParams', {sdparams}, ...
    'recordArrayParams', {arrparams}, ...
    'recordArrayParamLabels', {arrparamlabels}, ...
    'recordWaveformParams', {waveformparams}, ...
    'recordStringParams', {strparams}, ...
    'recordGroupParams', {groupparams}, ...
    'recordGroupParamLabels', {groupparamlabels}, ...
    'otherParams', {otherparams}, ...
    ...
    'errBarType', {'StandardError'}, ... % options: 'StandardDeviation', 'StandardError', '95%ConfidenceInterval' ... TODO: check to make sure these formulas are right
    'rotation', {15}, ...
    'useLogPlot', {0}, ... % uses log (1) or not (0) in the spectral plots
    'doNormalize', {0}, ... % normalizes the spectral plots column by column, if {1}, in the createGraphs*.m file
    'doAnova',{0}, ...
    'groupDivider', {0}, ... %if zero, plots all singledigitparams for N filename groups in a single bar graph,ungrouped.  If 2, splits the input data set into N/2 groups for the bar plot, if 3, N/3 groups, etc.
    'createGraphs', {1}, ...
    'colorMapSingle', {'gray'}, ... % or 'gray' or 'jet' or [0 0 0] or [1 1 1], etc
    'colorMapPercentage', {'gray'}, ... % or 'gray' or 'jet' or [0 0 0] or [1 1 1], etc
    'saveGraphs', {0}, ...
    'listData', {0} ...
    );
graphSettings.F_upperGraphLimit = computeSettings.F_upperGraphLimit;
graphSettings.F_lowerGraphLimit = computeSettings.F_lowerGraphLimit;
graphSettings.F_BaseNormalizingLimit = computeSettings.F_BaseNormalizingLimit;
graphSettings.F_LowerNormalizingLimit = computeSettings.F_LowerNormalizingLimit;
graphSettings.lineLocation = computeSettings.lineLocation;
graphSettings.processWaveformPSD = computeSettings.processWaveformPSD;
graphSettings.LFPFrequencyBands = computeSettings.LFPFrequencyBands;
graphSettings.LFPFrequencyBandLabels = computeSettings.LFPFrequencyBandLabels;
graphSettings.LFPNfft = computeSettings.LFPNfft;
