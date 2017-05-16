function [] = piv_series_standalone(param_file, version)
% function [] = piv_series_standalone(param_file, version)
%
% Wrapper function to run piv_series() as a compiled standalone executable.
% Includes special configuration needed to run on the BU SCC cluster batch
% system.
%
% Arguments:
%   param_file: String, .MAT file containing all piv_series() input arguments,
%       typically produced using a customized version of piv_get_param_template.m
%   version: String, optional, version hash for yalebox-piv, required as a
%       work-around for path problems on the SCC
% % 

% handle optional args
if nargin < 2 
    version = [];
end

% initialize environment for SCC cluster
if getenv('ENVIRONMENT')
    % SCC: avoid use of remote disk and use only requested number of cores
    fprintf('%s: Init parallel environment for SCC\n', mfilename);
    cluster = parcluster('local');
    cluster.JobStorageLocation = getenv('TMPDIR');
    nslots = str2double(getenv('NSLOTS'));
    parpool(cluster, nslots);
    maxNumCompThreads(nslots);
end

% run PIV analysis
args = load(param_file);
piv_series(args.piv_file, args.image_file, args.step_range, args.gap, ...
    args.samp_len, args.samp_spc, args.intr_len, args.num_pass, ...
    args.valid_radius, args.valid_max, args.valid_eps, ...
    args.min_frac_data, args.min_frac_overlap, true, version);
