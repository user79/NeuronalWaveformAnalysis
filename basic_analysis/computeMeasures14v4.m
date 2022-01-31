function [measuresStruct] = computeMeasures14v4(fileName, channels, manualData, settings)
% TPG created 10-10-07
% fileName is a char array containing a string with the fully qualified filename of
%       the .smr file being processed, (or the .txt file of spiketimes).
% channels is a numerical array of channel numbers of specific types to process.  "-1" in the "channels" array
%       means that this particular type of input is not expected to be found in this
%       file, so don't bother looking for it.
% manualData is a variable with hand-classified information, such as 
%       "bursty"/"irregular"/"regular" (3,2,1) classification scheme (4 is
%       unclassified).
% settings is a struct that contains lots of settings.  See defaultMultipleRunSettings.m
%       for more info.
%
% measuresStruct is a structure with all the output metrics, easy for subsequent graphing
%
%
%
% % Debug/example: uncomment these parameters and comment out the function definition at the top to use this file in script mode
% prefix = 'C:\Documents and Settings\tgilmour\My Documents\thesis\matlab code\'; 
% list = { ...
%     {[prefix 'tract12 dv8p6.smr']}, {[7,8,-1,2,1,4,3,-1,-1,-1,31]},{4};
%     };
% fileName = char(list{1});
% channels = list{2}{1,1};
% manualData = list{1,3};
% settings = defaultMultipleRunSettings14;


% *************************************************************************
% Initialize variables
% *************************************************************************
measuresStruct.Filename = fileName;
manualDataArray = double(manualData{1});
measuresStruct.Class = manualDataArray(1);  % Using Brigitte Piallat's classification scheme, either class 1 ("regular tonic"), class 2 ("irregular"), or class 3 ("irregular and bursty") or 4 for "doesn't fit the classification scheme"
measuresStruct.Rat = -1;
measuresStruct.EstimatedDepth = -1;
measuresStruct.SimpleListedDepth = -1;
if length(manualDataArray) > 1
    measuresStruct.Rat = manualDataArray(2);
    measuresStruct.EstimatedDepth = manualDataArray(3);
    measuresStruct.SimpleListedDepth = manualDataArray(4);
end
measuresStruct.SampleEntropy = -1*ones(1, settings.SEm - 1);
measuresStruct.SampleEntropyLabels = settings.SampleEntropyLabels;
measuresStruct.CoefficientOfVariation = -1;
measuresStruct.CVpm = -1;
measuresStruct.FiringRate = -1;
measuresStruct.BurstIndex = -1;
% measuresStruct.LFPPercentPowers = -1*ones(1,size(settings.LFPFrequencyBands,2));
measuresStruct.LFPPSD = -1;
measuresStruct.EEG1PSD = -1;
measuresStruct.EEG2PSD = -1;
measuresStruct.EEGdiffLFPCOH = -1;
measuresStruct.EEG1LFPCOH = -1;
measuresStruct.EEG2LFPCOH = -1;
measuresStruct.EEGdiffPSD = -1;
measuresStruct.normalizedLFPPSD = -1;
measuresStruct.LFPPSDf = -1;
measuresStruct.LFPFrequencyBands = settings.LFPFrequencyBands;
measuresStruct.LFPFrequencyBandLabels = settings.LFPFrequencyBandLabels;
% measuresStruct.LFPPSD_3to8Hz = -1;
% measuresStruct.LFPPSD_8to15Hz = -1;
% measuresStruct.LFPPSD_15to30Hz = -1;
% measuresStruct.LFPPSD_65to80Hz = -1;
measuresStruct.SpkPSD_3to8Hz = -1;
measuresStruct.SpkPSD_8to15Hz = -1;
measuresStruct.SpkPSD_gt15Hz = -1;
measuresStruct.EEG1FreqPeakRatios = -1;
measuresStruct.EEG1FreqSumRatios  = -1;
measuresStruct.EEG1FreqCentroid   = -1;
measuresStruct.NumBursts = -1;
measuresStruct.BurstRate = -1;
measuresStruct.MeanSpikesPerBurst = -1;
measuresStruct.MedianSpikesPerBurst = -1;
measuresStruct.MeanIntraburstFrequency = -1;
measuresStruct.MedianIntraburstFrequency = -1;
measuresStruct.ProportionTimeInBursts = -1;
measuresStruct.ProportionSpikesInBursts = -1;
measuresStruct.MeanSurprise = -1;
measuresStruct.MedianSurprise = -1;
measuresStruct.L = -1;
measuresStruct.rangeCounts = -1;
measuresStruct.DDH = -1*ones(1,settings.maxSpikesExpectedPerBin);
measuresStruct.DDHmean = -1;
measuresStruct.DDHvar = -1;
measuresStruct.STA = -1;
measuresStruct.STA12t25 = -1;
measuresStruct.STA25t40 = -1;
measuresStruct.EEG1STA = -1;
measuresStruct.EEG2STA = -1;
measuresStruct.EEGdiffSTA = -1;
measuresStruct.PtTDEEGdiff = -1;
measuresStruct.PtTRRDEEGdiff = -1;
measuresStruct.PtTDEEG1 = -1;
measuresStruct.PtTRRDEEG1 = -1;
measuresStruct.PtTDEEG2 = -1;
measuresStruct.PtTRRDEEG2 = -1;
measuresStruct.PtTD = -1;
measuresStruct.PtTRRD = -1;
measuresStruct.PtTD12t25 = -1;
measuresStruct.PtTRRD12t25 = -1;
measuresStruct.PtTD25t40 = -1;
measuresStruct.PtTRRD25t40 = -1;
measuresStruct.randomizedSTA12t25 = -1;
measuresStruct.randomizedSTA25t40 = -1;
measuresStruct.randomizedSTA = -1;
measuresStruct.randomizedEEG1STA = -1;
measuresStruct.randomizedEEG2STA = -1;
measuresStruct.randomizedEEGdiffSTA = -1;
measuresStruct.PeakToTroughNormalToRandomizedRatio = -1;
measuresStruct.PeakToTroughNormalToRandomizedRatio25t40 = -1;
measuresStruct.PeakToTroughNormalToRandomizedRatio12t25 = -1;
measuresStruct.PeakToTroughNormalToRandomizedRatioEEG1 = -1;
measuresStruct.PeakToTroughNormalToRandomizedRatioEEG2 = -1;
measuresStruct.PeakToTroughNormalToRandomizedRatioEEGdiff = -1;
measuresStruct.PeakEEGFreq = -1;
measuresStruct.EEGdiffLFPcorr = -1;
measuresStruct.EEGdiffLFPcorrZL = -1;
measuresStruct.EEG1LFPcorr = -1;
measuresStruct.EEG1LFPcorrZL = -1;
measuresStruct.EEG2LFPcorr = -1;
measuresStruct.EEG2LFPcorrZL = -1;
measuresStruct.PoissonComparison = -1;
measuresStruct.PoissonComparisonP = -1;
measuresStruct.PoissonComparisonSt = [];
measuresStruct.PoissonComparisonMean = -1;
measuresStruct.PoissonComparisonH = [];
measuresStruct.NumSpikes = -1;
measuresStruct.Duration = -1;
measuresStruct.CountAC = -1;
measuresStruct.Oscillatory = -1;
measuresStruct.LFPClipping = 0; % no typo: 0 means no amplifier saturation/clipping, 1 means yes clipping
measuresStruct.EEG1Clipping = 0;
measuresStruct.EEG2Clipping = 0;
measuresStruct.Date = datestr(now);
measuresStruct.ErrorString='';
measuresStruct.Fs_spikes = -1;
measuresStruct.Fs_lfp = -1;
measuresStruct.Fs_lfpDfs = -1;
measuresStruct.Fs_spikes_and_lfp = -1;
measuresStruct.starttime = -1;
measuresStruct.endtime = -1;

measuresStruct.bestScale = -1;
measuresStruct.arBestMSE = -1;
measuresStruct.arUndecimatedMSE = -1;
measuresStruct.MSEarray = -1;
measuresStruct.models = -1;
measuresStruct.sA = -1;
measuresStruct.sB = -1;
measuresStruct.mA = -1;
measuresStruct.mB = -1;
measuresStruct.rsAB = -1;
measuresStruct.rmAB = -1;
measuresStruct.lengthA = -1;
measuresStruct.lengthB = -1;



%measures = -1*ones(1,settings.NumberOfMeasures); % "-1" returned means that a particular measure has NOT been computed
if settings.verbose,
    disp('In ComputeMeasures14v4')     % cmSettings = settings
end

% *************************************************************************
% Load file and compute requested measures
% *************************************************************************
[fid, fileOpenMessage] = fopen(fileName, 'r'); % these are fully-qualified filenames, including directories
if fid ~= -1
    numSortedNeurons =0; %initialize
    measuresStruct.Fs_spikes = 0;
    measuresStruct.Fs_lfp = 0;
    measuresStruct.Fs_spikes_and_lfp = 0;
    measuresStruct.starttime=0;
    measuresStruct.endtime=0;
    if strcmpi(fileName((length(fileName)-2):length(fileName)),'smr')
        if (channels(1) ~= -1) % if there is a raw high speed highpass filtered spike waveform
		    % Uses the CED SON libraries, version 2.2, from Malcolm Lidierth
            [spikewavedata, header1]=SONGetChannel(fid, channels(1), 'seconds','scale'); % get the waveform (assumes that channel 1 is a "raw" ADC high-freq. waveform channel)
            measuresStruct.Fs_spikes = 1/(header1.sampleinterval * 1e-6); % should be ~ 25 kHz or 41kHz
            measuresStruct.endtime=header1.npoints/measuresStruct.Fs_spikes;
        end
        if (channels(2) ~= -1) % if there is a raw low speed lowpass filtered LFP waveform
            [lfpwavedata, header2]=SONGetChannel(fid, channels(2), 'seconds','scale');
            measuresStruct.Fs_lfp = 1/(header2.sampleinterval * 1e-6); % should be ~ 500 Hz
            measuresStruct.endtime=header2.npoints/measuresStruct.Fs_lfp;
        end
        if (channels(3) ~= -1) % if there is a raw high speed wideband filtered spikes+LFP waveform
            [spikesandlfpwavedata, header3]=SONGetChannel(fid, channels(3), 'seconds','scale');
            measuresStruct.Fs_spikes_and_lfp = 1/(header3.sampleinterval * 1e-6); % should be ~ 25 kHz or 41kHz
            measuresStruct.endtime=header3.npoints/measuresStruct.Fs_spikes_and_lfp;

            % TODO: add code to filter out the LFP signal separately, here...

        end
        if (channels(4) ~= -1) % if there is a raw low speed lowpass filtered EEG waveform
            [eeg1wavedata, headereeg1]=SONGetChannel(fid, channels(4), 'seconds','scale');
            measuresStruct.endtime=headereeg1.npoints/measuresStruct.Fs_lfp;   % assumes the same Fs for EEG and STN LFP
        end
        if (channels(5) ~= -1) % if there is a raw low speed lowpass filtered EEG waveform
            [eeg2wavedata, headereeg2]=SONGetChannel(fid, channels(5), 'seconds','scale');
            measuresStruct.endtime=headereeg2.npoints/measuresStruct.Fs_lfp; % assumes the same Fs for EEG and STN LFP
        end
        %... etc for the EEG channels, EMG channels, etc
        if (channels(7) ~= -1) % if there is one sorted spike channel
            [neuron1tsdata, header7]=SONGetChannel(fid, channels(7), 'seconds','scale');
            % get all of the spike timestamps on this sorted channel (assumes that
            % channel2 is a "Marker" or a "WaveMark" channel which contains the timestamp
            % info for ALL sorted spikes... but if settings.processOnlyFirstSortedSpike is
            % set to 1, this program will only analyze the spikes of code "01", i.e. the
            % first sorted spike only)
        end    %
        %... etc for more sorted spike channels
    elseif strcmpi(fileName((length(fileName)-2):length(fileName)),'txt')
        neuron1tsdata = load(fileName);
        measuresStruct.endtime = neuron1tsdata(end);
        measuresStruct.starttime = neuron1tsdata(1);
        settings.processWaveformPSD = 0; % if a .txt file, no waveform info - can't do PSD
    elseif strcmpi(fileName((length(fileName)-2):length(fileName)),'mat')
        c1 = load(fileName,['chan' int2str(channels(1))],['head' int2str(channels(1))],['chan' int2str(channels(2))], ['head' int2str(channels(2))], ['chan' int2str(channels(7))],['head' int2str(channels(7))], 'FileInfo');
        if exist('c1','var') && isstruct(c1)
            if isfield(c1,['chan' int2str(channels(1))])
                spikewavedata = double(c1.(['chan' int2str(channels(1))]) / c1.(['head' int2str(channels(1))]).max);
                measuresStruct.Fs_spikes = 1/(c1.(['head' int2str(channels(1))]).sampleinterval * 1e-6); % should be ~ 25 kHz or 41kHz
                measuresStruct.endtime=c1.(['head' int2str(channels(1))]).npoints/measuresStruct.Fs_spikes;
            end
            if isfield(c1,['chan' int2str(channels(2))])
                lfpwavedata = double(c1.(['chan' int2str(channels(2))]) / c1.(['head' int2str(channels(2))]).max);
                measuresStruct.Fs_lfp = 1/(c1.(['head' int2str(channels(2))]).sampleinterval * 1e-6); % should be ~ 25 kHz or 41kHz
                measuresStruct.endtime=c1.(['head' int2str(channels(2))]).npoints/measuresStruct.Fs_lfp;
            end
            if isfield(c1,['chan' int2str(channels(7))])
                neuron1tsdata = c1.(['chan' int2str(channels(7))]);
                if isfield(neuron1tsdata,'timings')
                    measuresStruct.endtime = neuron1tsdata.timings(end);
                    measuresStruct.starttime = neuron1tsdata.timings(1);
                end
            end
        end
    else
        measuresStruct.ErrorString = [measuresStruct.ErrorString '\nUnrecognized file format!'];
    end
    if settings.processOnlyFirstSortedSpike % This code assumes that the desired neuron is the one with sortcode "01"... this code ignores all others
        if (isfield(neuron1tsdata,'timings') && isfield(neuron1tsdata,'markers')) % the channel is a Spike2 "WaveMark" or "Marker" type channel
            temp1 = neuron1tsdata.timings;
            temp2 = neuron1tsdata.markers;
            numSortedNeurons = max(temp2(:,1));
            times01 = zeros(length(temp1),1);
            %times02 = zeros(length(temp1),1); % use these for 2nd, 3rd, etc sorted spikes/neurons
            counter01 = 1;
            %counter02 = 1;
            for tn1=1:length(temp1)
                if (temp2(tn1,1) == 01) % i.e. if this particular spike had sortcode "01"
                    times01(counter01) = temp1(tn1); % then copy its timestamp for further analysis
                    counter01 = counter01 + 1;
                    %elseif (temp2(tn1,1) == 02)
                    %    times02(counter02) = temp1(tn1); % then copy its timestamp for further analysis
                    %    counter02 = counter02 + 1;
                end
            end
            spiketimestamps = times01(1:counter01-1);
        else % the channel is a Spike2 "Event+" type channel or from a .txt file - there is no info about which spike came when, so might as well take all timestamps
            spiketimestamps = neuron1tsdata;
            numSortedNeurons = 1;
        end
    else
        % TODO: add code to process each sorted spike separately and
        % return the separate results
    end
    % Check for any errors before proceeding:  (try to weed out any weird problems with ascertaining the Fs, etc)
    if ((channels(1) ~= -1)&&(measuresStruct.Fs_spikes(1) < 1000)) || ((channels(2) ~= -1)&&(measuresStruct.Fs_lfp(1) < 200)) || ...
            ((channels(3) ~= -1)&&(measuresStruct.Fs_spikes_and_lfp(1) < 1000)) || (measuresStruct.endtime(1) == 0) || ...
            (size(measuresStruct.endtime,2) > 1) || isempty(spiketimestamps)  % the sampling frequency should be at least this, and typically much higher... also, for now, avoid multi-segment files
        measuresStruct.ErrorString = [measuresStruct.ErrorString sprintf('\n >>>>  Problematic file: measuresStruct.Fs_spikes=%g, measuresStruct.Fs_lfp=%g, segments=%g, measuresStruct.endtime=%g, ts=%g. Skipped ''%s''.  <<<<',measuresStruct.Fs_spikes,measuresStruct.Fs_lfp,size(measuresStruct.endtime,2),measuresStruct.endtime,length(spiketimestamps),fileName)];
        disp('! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ')
        disp(measuresStruct.ErrorString)
        disp('! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ')
        fclose(fid);
        return; % go to the next file
    end
    if (settings.spiketrainLengthLimit > 0)
        if (length(spiketimestamps) >= settings.spiketrainLengthLimit)
            spiketimestamps = spiketimestamps(1:settings.spiketrainLengthLimit); % crop to the first spiketrainLengthLimit spikes
            measuresStruct.endtime = spiketimestamps(end);
            measuresStruct.starttime = spiketimestamps(1);
        else
            disp(sprintf('Warning - the number of spikes %g is less than the ''spiketrainLengthLimit'' %g.', length(spiketimestamps), settings.spiketrainLengthLimit))
        end
    elseif (settings.spiketrainTimeLimit > 0) % this will not be checked if the spiketrainLengthLimit is < 0
        if ((spiketimestamps(end)-spiketimestamps(1)) >= settings.spiketrainTimeLimit)
            indexTimeLimit = find((spiketimestamps - spiketimestamps(1)) > settings.spiketrainTimeLimit, 1);
            spiketimestamps = spiketimestamps(1:indexTimeLimit-1); % crop to the first spiketrainLengthLimit spikes
            measuresStruct.endtime = spiketimestamps(end);
            measuresStruct.starttime = spiketimestamps(1);
        else
            disp(sprintf('Warning - the duration %g of the spiketrain is less than the ''spiketrainTimeLimit'' %g.', spiketimestamps(end)-spiketimestamps(1), settings.spiketrainTimeLimit))
        end
    end

%     % % make spike/EEG pics: special ad-hoc code for 'classifying' all neurons by 'hypothermic' eeg, 'SWS' or 'GA' slow wave sleep or global activation...
%     if exist('eeg1wavedata', 'var') && exist('eeg2wavedata', 'var') 
%         h=figure(1); clf(h)
%         extraheight=0.03;extrawidth=0.05;
%         h1=subplot(2,2,1); stem(spiketimestamps,ones(length(spiketimestamps),1)), title(fileName), axis([0 20 -.25 1.25])
%         h2=subplot(2,2,2); stem(spiketimestamps,ones(length(spiketimestamps),1)), axis([spiketimestamps(end)-20 spiketimestamps(end) -.25 1.25]),  title('last 20 sec')
%         if exist('eeg1wavedata', 'var') && exist('eeg2wavedata', 'var')
%             h3=subplot(2,2,3); plot(eeg1wavedata), axis([0 20000 0 1]), axis 'auto y', title('initial 20 sec')
%             h4=subplot(2,2,4); plot(eeg1wavedata), axis([spiketimestamps(end)*1000-20000 spiketimestamps(end)*1000 0 1]), axis 'auto y', title('last 20 sec')
%         end
%         p = get(h1, 'pos'); p(4) = p(4) + extraheight; p(3) = p(3) + 2*extrawidth; p(1) = p(1) - 2*extrawidth;
%         set(h1, 'pos', p);
%         p = get(h2, 'pos'); p(4) = p(4) + extraheight; p(3) = p(3) + 2*extrawidth; p(1) = p(1) - extrawidth;
%         set(h2, 'pos', p);
%         if exist('eeg1wavedata', 'var') && exist('eeg2wavedata', 'var')
%             p = get(h3, 'pos'); p(4) = p(4) + extraheight; p(3) = p(3) + 2*extrawidth; p(1) = p(1) - 2*extrawidth;
%             set(h3, 'pos', p);
%             p = get(h4, 'pos'); p(4) = p(4) + extraheight; p(3) = p(3) + 2*extrawidth; p(1) = p(1) - extrawidth;
%             set(h4, 'pos', p);
%         end
%         saveas(h,[settings.outputDirectoryPrefix 'hyp_' datestr(now, 30) '.png'])
%     end
%     %
%     % % % pause
%     % % % keyboard
% 
%     % make EEG and EEG-PSD pics for each neuron (that has EEG)
%     if exist('eeg1wavedata', 'var') && exist('eeg2wavedata', 'var') % special ad-hoc code for 'classifying' all neurons by 'hypothermic' eeg, 'SWS' or 'GA' slow wave sleep or global activation...
%         h=figure(1); clf(h)
%         extraheight=0.03;extrawidth=0.05;
%         %     [Pxx,F] = PWELCH(X,WINDOW,NOVERLAP,NFFT,Fs)
%         graphTime = 40;
%         NfftSWS = 4096;
%         [PsdTemp1,Ftemp1] = pwelch(eeg1wavedata(1:graphTime*measuresStruct.Fs_lfp),[],128,4*NfftSWS,measuresStruct.Fs_lfp);
%         [PsdTemp2,Ftemp2] = pwelch(eeg1wavedata,[],128,NfftSWS,measuresStruct.Fs_lfp);
%         [maxPeak,peakIndex] = max(PsdTemp2);
%         measuresStruct.PeakEEGFreq = Ftemp2(peakIndex);
%         [PsdTemp3,Ftemp3] = pwelch(eeg1wavedata(length(eeg1wavedata)-graphTime*measuresStruct.Fs_lfp:length(eeg1wavedata)),[],128,4*NfftSWS,measuresStruct.Fs_lfp);
%         h1=subplot(2,3,1); plot(Ftemp1,PsdTemp1), title(fileName), axis([0 10 -.25 1.25]), axis 'auto y'
%         h2=subplot(2,3,2); plot(Ftemp2,PsdTemp2), axis([0 10 -.25 1.25]), axis 'auto y'
%         h3=subplot(2,3,3); plot(Ftemp3,PsdTemp3), axis([0 10 -.25 1.25]), axis 'auto y'
%         h4=subplot(2,3,4); plot([1/(graphTime*measuresStruct.Fs_lfp):1/(graphTime*measuresStruct.Fs_lfp):1]*graphTime,eeg1wavedata(1:graphTime*measuresStruct.Fs_lfp)), axis([0 graphTime 0 1]), axis 'auto y', title('initial 20 sec')
%         h5=subplot(2,3,5); plot((1:length(eeg1wavedata))/measuresStruct.Fs_lfp,eeg1wavedata), axis([0 length(eeg1wavedata)/measuresStruct.Fs_lfp 0 1]), axis 'auto y',
%         h6=subplot(2,3,6); plot([1/(graphTime*measuresStruct.Fs_lfp):1/(graphTime*measuresStruct.Fs_lfp):1]*graphTime,eeg1wavedata(length(eeg1wavedata)-graphTime*measuresStruct.Fs_lfp+1:length(eeg1wavedata))), , axis([0 graphTime 0 1]), axis 'auto y', title('final 20 sec')
%     %keyboard
%         p = get(h1, 'pos'); p(4) = p(4) + extraheight; p(3) = p(3) + 2*extrawidth; p(1) = p(1) - 2*extrawidth;
%         set(h1, 'pos', p);
%         p = get(h2, 'pos'); p(4) = p(4) + extraheight; p(3) = p(3) + extrawidth; p(1) = p(1) - extrawidth;
%         set(h2, 'pos', p);
%         p = get(h3, 'pos'); p(4) = p(4) + extraheight; p(3) = p(3) + 2*extrawidth; p(1) = p(1) - extrawidth;
%         set(h3, 'pos', p);
%         p = get(h4, 'pos'); p(4) = p(4) + extraheight; p(3) = p(3) + 2*extrawidth; p(1) = p(1) - 2*extrawidth;
%         set(h4, 'pos', p);
%         p = get(h5, 'pos'); p(4) = p(4) + extraheight; p(3) = p(3) + extrawidth; p(1) = p(1) - extrawidth;
%         set(h5, 'pos', p);
%         p = get(h6, 'pos'); p(4) = p(4) + extraheight; p(3) = p(3) + 2*extrawidth; p(1) = p(1) - extrawidth;
%         set(h6, 'pos', p);
%         saveas(h,[settings.outputDirectoryPrefix 'sws_' datestr(now, 30) '.png'])
%     end
%     %
%     % % % pause

    
    
    

    % *************************************************************************
    % Spike Timestamp-based processing
    % *************************************************************************
    if settings.processBasic
        if(settings.verbose) disp('processing FiringRate and CoefficientOfVariation...'); end
        isi = diff(spiketimestamps);  % TODO: add checks for outliers, from the gaps between multiple separated recording sessions
        stdISI = std(isi);
        meanISI = mean(isi);
        measuresStruct.CoefficientOfVariation = stdISI/meanISI;
        measuresStruct.NumSpikes = length(spiketimestamps);
        measuresStruct.Duration = (sum(measuresStruct.endtime) - measuresStruct.starttime);
        measuresStruct.CVpm = measuresStruct.CoefficientOfVariation / (sqrt(measuresStruct.NumSpikes - 2) * (1 - (measuresStruct.NumSpikes - 1)*(settings.CVpmRefractoryPeriod / measuresStruct.Duration)));
        measuresStruct.FiringRate = measuresStruct.NumSpikes/measuresStruct.Duration; % Firing Rate for the recorded data segment(s)
        %modeISI = max(isi);
        edges = 0:0.0005:settings.maxISI; % 10ms bins for finding the mode of the ISI histogram for the burst index
        [nhistc,binhistc] = histc(isi,edges);
        [maxisicount,modeISIindex] = max(nhistc);
        if modeISIindex < length(edges)
            measuresStruct.BurstIndex = meanISI / edges(modeISIindex+1); % Burst index
        end
    end
    if settings.processSE
        if(settings.verbose) disp('processing SampleEntropy...'); end
        if settings.SampEnStandardLength > 0 %if settings.SampEnStandardLength == 0, then use ALL spikes, otherwise, check for length
            if length(isi) >= settings.SampEnStandardLength
                SampEntropy = sampen(isi(1:settings.SampEnStandardLength), settings.SEm, settings.SEr);  % m=settings.SEm (default 5, returns 4 values), r=0.2*STD,
                % but this function normalizes the ISI sequence by default
                % to unit standard deviation, so r is effectively 0.2*std(isi)
                % see http://www.physionet.org/physiotools/sampen/ for more details
                measuresStruct.SampleEntropy = SampEntropy';
            else
                measuresStruct.ErrorString = [measuresStruct.ErrorString '\nToo few spikes for SampleEntropy.  Skipped file.'];
            end
        else
            SampEntropy = sampen(isi, settings.SEm, settings.SEr);
            measuresStruct.SampleEntropy = SampEntropy'; %vector: returns multiple values in a row array
        end
    end
    if settings.processSpikePSD  % TODO: pass in various frequency ranges in an array, and return the normalized spectral power in those ranges in another array
        if(settings.verbose) disp('processing spectral power from ISI''s...'); end   %this is TomWichmann's procedure for looking at spike oscillations
        [b2_3_to_8, b2_8_to_15, b2_gt_15, totalpowerISI, failureFlag, oscillatory] = computePSDfromISIs14v4(isi);
        if ~failureFlag
            measuresStruct.SpkPSD_3to8Hz = b2_3_to_8/totalpowerISI*100;
            measuresStruct.SpkPSD_8to15Hz = b2_8_to_15/totalpowerISI*100;
            measuresStruct.SpkPSD_gt15Hz = b2_gt_15/totalpowerISI*100;
            measuresStruct.Oscillatory = oscillatory;
        end
    end
    if settings.processBursts
        if(settings.verbose) disp('processing bursts...'); end
        c1 = [];
        if length(isi) > 2
            % First do Poisson-comparison (Levy 2001 "Effects of Apomorphine on Subthalamic Nucleus and Globus
            % Pallidus Internus Neurons in Patients With Parkinson's Disease" modification of Kaneoke and Vitek
            % 1996 "Burst and oscillation as disparate neuronal properties" method
            [poissontest,mtemp,ptemp,sttemp,ddhtemp,Htemp,ctemp,Ltemp,rangeCountstemp,mu,variance]  = computeDDH(spiketimestamps, settings);
            measuresStruct.PoissonComparison = poissontest;
            measuresStruct.PoissonComparisonP = ptemp;
            measuresStruct.PoissonComparisonMean = mtemp;
            measuresStruct.PoissonComparisonH = Htemp;
            measuresStruct.PoissonComparisonSt = sttemp;
            measuresStruct.L = Ltemp;
            measuresStruct.rangeCounts = rangeCountstemp;
            measuresStruct.DDH = ddhtemp;
            measuresStruct.DDHmean = mu;
            measuresStruct.DDHvar = variance;
            % Next, detect actual bursts and do processing on them
            switch settings.BurstDetectionMethod
                case 'HanesPoissonSurprise'
                    disp('HanesPoissonSurprise method no longer used')
                case 'WichmannPoissonSurprise'
                    Bursts = legendy_new3(isi,settings.burstFactor,settings.burstSamplingRate,settings.minBurstLength,settings.burstLocalLength,settings.burstSurpriseCutoff);
                    if Bursts(1).num_bursts > 0
                        measuresStruct.NumBursts = Bursts(1).num_bursts;
                        measuresStruct.BurstRate = Bursts(1).num_bursts / (sum(measuresStruct.endtime) - measuresStruct.starttime);
                        measuresStruct.MeanSpikesPerBurst = Bursts(1).mean_spikes_per_burst;
                        measuresStruct.MedianSpikesPerBurst = Bursts(1).median_spikes_per_burst;
                        measuresStruct.MeanIntraburstFrequency = Bursts(1).mean_intra_burst_frequency;
                        measuresStruct.MedianIntraburstFrequency = Bursts(1).median_intra_burst_frequency;
                        measuresStruct.ProportionTimeInBursts = Bursts(1).proportion_time_in_bursts;
                        measuresStruct.ProportionSpikesInBursts = Bursts(1).proportion_spikes_in_bursts;
                        measuresStruct.MeanSurprise = mean([Bursts.surprise]);
                        measuresStruct.MedianSurprise = median([Bursts.surprise]);
                    else
                        measuresStruct.NumBursts = 0;
                        measuresStruct.BurstRate = 0;
                        measuresStruct.ProportionTimeInBursts = 0;
                        measuresStruct.ProportionSpikesInBursts = 0;
                    end
                    % figure; hold on
                    % stem(spiketimestamps,ones(length(spiketimestamps),1),'k')
                    % stem(spiketimestamps([Bursts.begin]),[Bursts.surprise],'b')
                    % hold off
                case 'RankSurprise'
                    disp('RankSurprise method no longer used')
            end
        end
    end
    if settings.processISIHistsAndACs && settings.processBasic
        [isihist,isihistcenters] = hist(isi, settings.isiHistogramNumOfBins);
        [countAC, bin_centers, ltmean, threshold, cgram, smoothed_cgram] = countACpeaks01(spiketimestamps, settings.countACpeaksTime, 1/settings.isiHistogramNumOfBins);
        % TODO: add oscillation-checks here
        if settings.interactiveACmode
            figure(1);
            subplot(3,1,1);
            hist(isi, settings.isiHistogramNumOfBins);
            xlim([0 settings.isiHistogramTime])
            title(fileName);
            %shg; %pause; %close(1);
            subplot(3,1,2);
            bar(bin_centers, cgram)
            subplot(3,1,3);
            bar(bin_centers, smoothed_cgram)
            lc = length(smoothed_cgram);
            line(bin_centers(end)/lc:bin_centers(end)/lc:bin_centers(end),threshold*ones(1,lc))
            line(bin_centers(end)/lc:bin_centers(end)/lc:bin_centers(end),ltmean*ones(1,lc),'LineStyle','--')
            %fprintf('\n Count = %g',count)
            xlabel([int2str(countAC) ' peaks?']);
            manually_classified_num_peaks = input(sprintf('\nHow many peaks in autocorrelation? (guess=%g peaks)',countAC));
            if ~isempty(manually_classified_num_peaks)
                countAC = manually_classified_num_peaks;
            end
            %shg; %pause %close(1);
            if (measuresStruct.CoefficientOfVariation < 0.5) && (countAC >= 3)  %this section will override the section above in processBasic
                measuresStruct.Class = 1;
            elseif ((measuresStruct.CoefficientOfVariation >= 0.5) && (measuresStruct.CoefficientOfVariation < 0.8)) && ((countAC == 1) || (countAC == 2))
                measuresStruct.Class = 2;
            elseif (measuresStruct.CoefficientOfVariation >= 0.8) && (countAC == 0)
                measuresStruct.Class = 3;
            else
                measuresStruct.Class = 4; %i.e. doesn't fit the other preset categories
            end
        end
        measuresStruct.CountAC = countAC;
        if settings.saveACinfoFiles
            fn_ac = [settings.outputDirectoryPrefix 'ac' datestr(now, 30)];
            if length(spiketimestamps) > 500
                croppedTS = spiketimestamps(1:500);
            else
                croppedTS = spiketimestamps;
            end
            save(fn_ac, 'countAC', 'measuresStruct', 'settings', 'cgram', 'croppedTS', 'smoothed_cgram', 'bin_centers', 'ltmean', 'threshold','isihist','isihistcenters');
        end
    end

    % *************************************************************************
    % Waveform-based processing
    % *************************************************************************
    if settings.processWaveformPSD
        if exist('lfpwavedata', 'var') 
            if(settings.verbose) disp(sprintf('Processing spectral power in waveforms... Fs=%g.',measuresStruct.Fs_lfp)); end
            % % % %             h = fdesign.lowpass(settings.Fpass, settings.Fstop, settings.Apass, settings.Astop, measuresStruct.Fs_lfp/settings.DecimationFactor);
            % % %             h = fdesign.bandpass(settings.Fstop1, settings.Fpass1, settings.Fpass2, settings.Fstop2, settings.Astop, settings.Apass, settings.Astop, measuresStruct.Fs_lfp/settings.DecimationFactor);
            % % %             settings.HdpreprocessWaveforms = design(h, 'cheby1', 'MatchExactly', 'passband'); % make lowpass filter, given this Fs
            [decimatedlfpwavedata, dFs, clipping] = preprocessWaveforms14v4(lfpwavedata, measuresStruct.Fs_lfp, settings); %filter and decimate (from, e.g. 1000 samples/sec to 200 samples/sec)
            %
            [STA12t25,rSTA12t25,PtTD12t25,PtTRRD12t25,PeakToTroughNormalToRandomizedRatio12t25] = computeSTA(spiketimestamps, decimatedlfpwavedata, dFs, settings.HdSTA12t25);
            measuresStruct.STA12t25 = STA12t25;
            measuresStruct.PtTD12t25 = PtTD12t25;
            measuresStruct.PtTRRD12t25 = PtTRRD12t25;
            measuresStruct.randomizedSTA12t25 = rSTA12t25;
            measuresStruct.PeakToTroughNormalToRandomizedRatio12t25 = PeakToTroughNormalToRandomizedRatio12t25;
            [STA25t40,rSTA25t40,PtTD25t40,PtTRRD25t40,PeakToTroughNormalToRandomizedRatio25t40] = computeSTA(spiketimestamps, decimatedlfpwavedata, dFs, settings.HdSTA25t40);
            measuresStruct.STA25t40 = STA25t40;
            measuresStruct.PtTD25t40 = PtTD25t40;
            measuresStruct.PtTRRD25t40 = PtTRRD25t40;
            measuresStruct.randomizedSTA25t40 = rSTA25t40;
            measuresStruct.PeakToTroughNormalToRandomizedRatio25t40 = PeakToTroughNormalToRandomizedRatio25t40;
            [STA,rSTA,PtTD,PtTRRD,PeakToTroughNormalToRandomizedRatio] = computeSTA(spiketimestamps, decimatedlfpwavedata, dFs); % no filtering
            measuresStruct.STA = STA;
            measuresStruct.PtTD = PtTD;
            measuresStruct.PtTRRD = PtTRRD;
            measuresStruct.randomizedSTA = rSTA;
            measuresStruct.PeakToTroughNormalToRandomizedRatio = PeakToTroughNormalToRandomizedRatio;
            %
            [LFPPercentPowers, totalpowerW, PSDcp, PSDfcp, normPSD] = computePSDfromWaveforms14v4(decimatedlfpwavedata, dFs, settings);
            %measuresStruct.LFPPercentPowers = LFPPercentPowers;
            measuresStruct.LFPPSD = PSDcp;
            measuresStruct.Fs_lfpDfs = dFs;
            measuresStruct.LFPClipping = clipping;
            measuresStruct.normalizedLFPPSD = normPSD;
            measuresStruct.LFPPSDf = PSDfcp; % uses the decimated dFs, i.e. 200, not 1000
        else
            measuresStruct.ErrorString = [measuresStruct.ErrorString '\nNo channel with lfp waveform data found... LFPPSD not computed.'];
        end
        if exist('eeg1wavedata', 'var') && exist('eeg2wavedata', 'var')  % only my recordings
            [decimatedeeg1wavedata, dFs, clipping] = preprocessWaveforms14v4(eeg1wavedata, measuresStruct.Fs_lfp, settings);
            [PercentPowers, totalpowerW, PSDcp, PSDfcp, normPSD] = computePSDfromWaveforms14v4(decimatedeeg1wavedata, dFs, settings);
            measuresStruct.EEG1PSD = PSDcp;
            measuresStruct.EEG1Clipping = clipping;
            [decimatedeeg2wavedata, dFs, clipping] = preprocessWaveforms14v4(eeg2wavedata, measuresStruct.Fs_lfp, settings);
            [PercentPowers, totalpowerW, PSDcp, PSDfcp, normPSD] = computePSDfromWaveforms14v4(decimatedeeg2wavedata, dFs, settings);
            measuresStruct.EEG2PSD = PSDcp;
            measuresStruct.EEG2Clipping = clipping;
            eegdiff = eeg1wavedata(1:min(length(eeg1wavedata),length(eeg2wavedata))) - eeg2wavedata(1:min(length(eeg1wavedata),length(eeg2wavedata)));
            [decimatedeegdiffwavedata, dFs, clipping] = preprocessWaveforms14v4(eegdiff, measuresStruct.Fs_lfp, settings);
            [PercentPowers, totalpowerW, PSDcp, PSDfcp, normPSD] = computePSDfromWaveforms14v4(decimatedeegdiffwavedata, dFs, settings);
            measuresStruct.EEGdiffPSD = PSDcp;
            
            [EEGSTA,rEEGSTA,PtTD,PtTRRD,PeakToTroughNormalToRandomizedRatioEEG] = computeSTA(spiketimestamps, decimatedeeg1wavedata, dFs); % no filtering
            measuresStruct.EEG1STA = EEGSTA;
            measuresStruct.PtTDEEG1 = PtTD;
            measuresStruct.PtTRRDEEG1 = PtTRRD;
            measuresStruct.randomizedEEG1STA = rEEGSTA;
            measuresStruct.PeakToTroughNormalToRandomizedRatioEEG1 = PeakToTroughNormalToRandomizedRatioEEG;
            [EEGSTA,rEEGSTA,PtTD,PtTRRD,PeakToTroughNormalToRandomizedRatioEEG] = computeSTA(spiketimestamps, decimatedeeg2wavedata, dFs); % no filtering
            measuresStruct.EEG2STA = EEGSTA;
            measuresStruct.PtTDEEG2 = PtTD;
            measuresStruct.PtTRRDEEG2 = PtTRRD;
            measuresStruct.randomizedEEG2STA = rEEGSTA;
            measuresStruct.PeakToTroughNormalToRandomizedRatioEEG2 = PeakToTroughNormalToRandomizedRatioEEG;
            [EEGSTA,rEEGSTA,PtTD,PtTRRD,PeakToTroughNormalToRandomizedRatioEEG] = computeSTA(spiketimestamps, decimatedeegdiffwavedata, dFs); % no filtering
            measuresStruct.EEGdiffSTA = EEGSTA;
            measuresStruct.PtTDEEGdiff = PtTD;
            measuresStruct.PtTRRDEEGdiff = PtTRRD;
            measuresStruct.randomizedEEGdiffSTA = rEEGSTA;
            measuresStruct.PeakToTroughNormalToRandomizedRatioEEGdiff = PeakToTroughNormalToRandomizedRatioEEG;

            measuresStruct.EEGdiffLFPCOH = mscohere(decimatedeegdiffwavedata(1:min(length(decimatedeegdiffwavedata),length(decimatedlfpwavedata))),decimatedlfpwavedata(1:min(length(decimatedeegdiffwavedata),length(decimatedlfpwavedata))),settings.LFPNfft,settings.LFPNoverlap,[],dFs);
            measuresStruct.EEG1LFPCOH = mscohere(decimatedeeg1wavedata(1:min(length(decimatedeeg1wavedata),length(decimatedlfpwavedata))),decimatedlfpwavedata(1:min(length(decimatedeeg1wavedata),length(decimatedlfpwavedata))),settings.LFPNfft,settings.LFPNoverlap,[],dFs);
            measuresStruct.EEG2LFPCOH = mscohere(decimatedeeg2wavedata(1:min(length(decimatedeeg2wavedata),length(decimatedlfpwavedata))),decimatedlfpwavedata(1:min(length(decimatedeeg2wavedata),length(decimatedlfpwavedata))),settings.LFPNfft,settings.LFPNoverlap,[],dFs);

            % todo: refine / add coherence, TE, MI, etc more waveform measures and EEG-STN similarity-measures here
            
            [EEGFreqPeakRatios, EEGFreqSumRatios, EEGFreqCentroid] = computeEEGFreqRatios14v4(measuresStruct.EEG1PSD, dFs, settings); % to quantify SWS/GA separation
            measuresStruct.EEG1FreqPeakRatios = EEGFreqPeakRatios;
            measuresStruct.EEG1FreqSumRatios  = EEGFreqSumRatios;
            measuresStruct.EEG1FreqCentroid   = EEGFreqCentroid;
            
            % Correlation (todo: maybe don't need to save the actual arrays, just the zero-lag (max) ?  )  
            % ( TODO: split the file into chunks, take correlations of each chunk, return averages and standard deviations)
            segmentLengthToProcess = 20; % in seconds (todo: bring this out into the settings struct)
            tempLFPdata = decimatedlfpwavedata(1:segmentLengthToProcess*dFs);  % take the beginning of each file...
            tempEEG1data = decimatedeeg1wavedata(1:segmentLengthToProcess*dFs);  % take the beginning of each file...
            tempEEG2data = decimatedeeg2wavedata(1:segmentLengthToProcess*dFs);  % take the beginning of each file...
            tempEEGdiffdata = decimatedeegdiffwavedata(1:segmentLengthToProcess*dFs);  % take the beginning of each file...
            measuresStruct.EEGdiffLFPcorr = xcorr(tempLFPdata, tempEEGdiffdata);
            measuresStruct.EEGdiffLFPcorrZL = measuresStruct.EEGdiffLFPcorr(segmentLengthToProcess*dFs); % save the zero-lag value separately for graphing
            measuresStruct.EEG1LFPcorr = xcorr(tempLFPdata, tempEEG1data);
            measuresStruct.EEG1LFPcorrZL = measuresStruct.EEG1LFPcorr(segmentLengthToProcess*dFs);
            measuresStruct.EEG2LFPcorr = xcorr(tempLFPdata, tempEEG2data);
            measuresStruct.EEG2LFPcorrZL = measuresStruct.EEG2LFPcorr(segmentLengthToProcess*dFs);
                         
            if settings.processAdaptiveFilters
                % % LMS - see how well the adaptive filter can adapt/track
                
            end
        else
            measuresStruct.ErrorString = [measuresStruct.ErrorString '\nNo channel with relevant eeg waveform data found... EEG PSDs not computed.'];
        end
    end

    % *************************************************************************
    % Cleanup
    % *************************************************************************
    fclose('all'); %fclose(fid); % close file, in normal case
    if (settings.verbose)
        measuresStruct
        fprintf('\nFinished %s, %g spikes, %gs, (%g neurons)',fileName,length(spiketimestamps),(measuresStruct.endtime - measuresStruct.starttime),numSortedNeurons);
        % save([settings.outputDirectoryPrefix 'debug_cM_' datestr(now, 30)]);  % warning: VERY LARGE file
    end
else % if file-opening error (all the measures will simply be returned as "-1")
    measuresStruct.ErrorString = fileOpenMessage
end
