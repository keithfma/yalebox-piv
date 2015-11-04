% Script. Experimental attempt to account for the location of the data centroid
% for sample points near the wedge boundary. The approach is to interpolate from
% an irregular grid defined by the data centroids to the regular sample grid.

%% Initialize 

% define parameters
input_data = 'fault_ss_01_sidef_250_251.mat';
sample_len = 30;
sample_spc = sample_len;
interogation_len = 60;

% read in data
load(input_data, 'ini', 'fin');

% generate regular sample grid
[rr, cc] = yalebox_piv_sample_grid(sample_len, sample_spc, size(ini));


%% Get irregular grid

nc = length(cc);
nr = length(rr);

rr_cntr = nan(nr, nc);
cc_cntr = nan(nr, nc);

% loop over sample grid
for jj = 1:nc
    for ii = 1:nr
        
        % get sample and interrogation windows
        [samp, samp_pos] = yalebox_piv_window(ini, rr(ii), cc(jj), sample_len);
        [intr, intr_pos] = yalebox_piv_window(fin, rr(ii), cc(jj), interogation_len);
        
        if all(samp == 0)
            continue
        end
        
        % get centroid for sample window
        [data_rr, data_cc] = find(samp~=0);
        data_cntr = sum([data_rr, data_cc])/length(data_rr);
        rr_cntr(ii, jj) = data_cntr(1)-1+samp_pos(2);
        cc_cntr(ii, jj) = data_cntr(2)-1+samp_pos(1);
        
    end % ii
end % jj

imagesc(ini);
hold on
plot(cc_cntr, rr_cntr, 'xk');
hold off