function [] = stretch_color(img, mask)
% Stretch colors to make full use of the available range. Threshold values
% are computed as a function of x using a sliding window quantile filter.

% Arguments
upper_limit = 0.975;
lower_limit = 0.025;
bin_width = 50; % pixels

%% setup bins
bin_center = 1:bin_width:size(img,2);

bin_start = bin_center-bin_width/2;
bin_start = max(bin_start, 1);
bin_start = floor(bin_start);

bin_stop = bin_start+bin_width;
bin_stop = min(bin_stop, size(img,2));
bin_stop = ceil(bin_stop);

bin_start(end) = bin_stop(end)-bin_width;

bin_center = (bin_start+bin_stop)/2;


%% get bin quantiles at upper_limit and lower_limit 

img = rgb2gray(img);
upper = nan(size(bin_center));
lower = nan(size(bin_center));
for i = 1:numel(bin_center)
    img_sub = img(:, bin_start(i):bin_stop(i));
    mask_sub = mask(:, bin_start(i):bin_stop(i));
    bin_data = img_sub(mask_sub);
    limits = quantile(bin_data(:), [lower_limit, upper_limit]);
    lower(i) = limits(1);
    upper(i) = limits(2);    
end



%% correct range

keyboard





