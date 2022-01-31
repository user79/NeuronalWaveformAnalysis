function [count, bin_centers, ltmean, threshold, cgram, smoothed_cgram] = countACpeaks01(timestamps, width, bin)
% Pseudocode algorithm:
%   take autocorrelogram of timestamps
%   smooth autocorrelogram
%   count number of times the smoothed autocorrelogram rises above a prespecified multiple
%     of the distance between the maximum peak (usually the first peak) and the 
%     long-time-mean (average of the second half of the autocorrelogram data)  
%     (alternatively, use a multiple of the
%     standard deviation for the threshold...?)

msc_multiple = 0.3;
base_threshold = 0.3;
minMeanPeakDiff = 3;
distanceToSearchForTrough = 60;
smoothingWindowSize = 39;
gaussian_window =  gausswin(smoothingWindowSize); %[0.0439    0.2494    0.7066    1.0000    0.7066    0.2494    0.0439]; % gausswin(7);

[bin_edges,bin_centers,cgram] = cor_gram_Hz(timestamps,timestamps,width,bin,0,1);
lc = length(cgram);
smoothed_cgram = conv(gaussian_window, cgram);
ltmean = mean(smoothed_cgram(ceil(lc/2):lc));
% ltstd = std(smoothed_cgram(ceil(lc/2):lc))
[firstpeak,ix_firstpeak] = max(smoothed_cgram);  % called "first", but it's actually just the largest
if (ix_firstpeak + distanceToSearchForTrough) < lc
    [firsttrough,ix2_firsttrough] = min(smoothed_cgram(ix_firstpeak:ix_firstpeak + distanceToSearchForTrough)); % called "first", but it's actually just the largest
else
    [firsttrough,ix2_firsttrough] = min(smoothed_cgram(ix_firstpeak:end)); 
end
ix_firsttrough = ix2_firsttrough + ix_firstpeak;

threshold = base_threshold + ltmean + (firstpeak - ltmean) * msc_multiple;
count = sum(abs(diff(smoothed_cgram >= threshold)))/2;
if ((firstpeak - ltmean) <= minMeanPeakDiff) || ((ltmean - firsttrough) <= minMeanPeakDiff)
    count = 0;
end
smoothed_cgram = smoothed_cgram((smoothingWindowSize-1)/2:end-(smoothingWindowSize+1)/2); % eliminate the extra convolved tails

% figure(1)
% subplot(2,1,1)
% bar(cgram)
% subplot(2,1,2)
% bar(smoothed_cgram)
% line(1:lc,threshold*ones(1,lc))
% line(1:lc,ltmean*ones(1,lc),'LineStyle','--')


% junk
%
% baseline_std_threshold_multiple = 3;
% extra_std_threshold_multiple = 3/(msc-ltmean)
% if extra_std_threshold_multiple < 0
%     extra_std_threshold_multiple = 0;
% end
% threshold = ltmean + baseline_std_threshold_multiple*ltstd - extra_std_threshold_multiple*ltstd;

% line(1:lc,ltmean - (msc - ltmean)*extra_std_threshold_multiple*ltstd*ones(1,lc))

% figure(2); plot(abs(fft(smoothed_cgram-mean(smoothed_cgram))))
% 

% 
