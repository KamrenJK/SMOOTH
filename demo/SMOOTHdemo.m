% -------------------------------------------------------------------------
% demo script for SMOOTH group-level stats
%
% demonstrates a minimal end-to-end workflow:
%   1) set up paths (toolbox + FieldTrip)
%   2) load an example FieldTrip "source" dataset
%   3) run SMOOTHstat (group-level inference)
%   4) visualize the resulting t-map and significance mask
%   5) evaluate surrogate data
%
% requirements:
%   - FieldTrip
%   - FreeSurfer installation (for fsaverage surfaces)
%   - example dataset: demo/source.mat (distributed with the toolbox)
%
% Authors: Kamren Khan & Arjen Stolk, v0.1
% -------------------------------------------------------------------------

clear; clc; close all;

% set paths -- edit smooth_setup.m once (repository root); nothing to change here.
addpath(fullfile(fileparts(mfilename('fullpath')),'..'));
paths = smooth_setup();

% cortical mesh (only for visualization)
cortex      = load_fsaverage_mesh(paths.freesurfer);

% read in data
load(fullfile(paths.demo, 'source.mat'))  % example dataset
nsubs = numel(source);
nchans = sum(cellfun(@(s) numel(s.label), source));
fprintf('\nLoaded demo dataset: %d subjects, %d total electrodes.\n', nsubs, nchans);

%% run SMOOTH
tic
cfg                     = [];
cfg.fshome              = paths.freesurfer;
cfg.numrandomization    = 1000;
cfg.randomseed          = 42;      % fixed seed -> reproducible p-values
cfg.kernelwidth         = 12;      % mm FWHM
cfg.graphsigma          = 8;
cfg.smooth              = 'sphere';
cfg.kernel              = 'gaussian';
cfg.rankrescale         = 'exact';
cfg.normalize           = 'no';
cfg.tail                = 'both';  % 'positive' for a directional hypothesis
cfg.keepmaps            = 'yes';
% Surrogates are only needed for the SMOOTHdiag QC step at the end, and they
% dominate memory use. The subject array is nsubjects x nvertices x nsurrsamples,
% and SMOOTHdiag permutes it internally, which transiently doubles it:
%
%   nsurrsamples =  10  ->  ~0.9 GB  (~1.8 GB during the permute)
%   nsurrsamples =  25  ->  ~2.3 GB  (~4.6 GB during the permute)
%   nsurrsamples = 100  ->  ~9.3 GB  (~19  GB during the permute)
%
% plus ~2.4 GB for the group surrogates at 1000 permutations. 10 is ample for a
% smoothness check. Set keepsubsurrogates='no' to skip the QC section entirely.
cfg.keepsubsurrogates   = 'yes';
cfg.keepgroupsurrogates = 'yes';
cfg.subsamplesurr       = 'yes';
cfg.nsurrsamples        = 10;

stat = SMOOTHstat(cfg, source{:});
t = toc;
fprintf('SMOOTHstat completed in %.2f s.\n', t);

% significant clusters
sig = find([stat.posclusters.prob] <= 0.05);
fprintf('%d significant positive cluster(s)\n', numel(sig));
for k = sig
    fprintf('   cluster %d: p = %.4f, mass = %.1f\n', ...
        k, stat.posclusters(k).prob, stat.posclusters(k).clusterstat);
end

%% visualize output

% (1) group-level t-map
figure
ft_plot_mesh(cortex, 'vertexcolor', 'curv');
hold on; ft_plot_mesh(cortex, 'vertexcolor', stat.stat);
clim([-3 3]); colormap('parula')
title('group-level t-map'); view([90 0])

% (2) significant effects
figure
y = stat.stat;  % pos effects
y(~stat.mask | stat.stat<0) = nan;
ft_plot_mesh(cortex, 'vertexcolor', 'curv');
hold on; ft_plot_mesh(cortex, 'vertexcolor', y);
clim([-3 3]); colormap('parula')
title('SMOOTH significant effects'); view([90 0])

%% surrogate robustness check
% verify that surrogate spatial autocorr is equivalent to that of empirical
% data

% surrogate diagnostics
d    = SMOOTHdiag(stat);
emp  = d.subject.spatial.BOTH.moranI.per_subject.emp(:);
null = d.subject.spatial.BOTH.moranI.per_subject.null_mu(:);
[~,p,~,tstat] = ttest(emp, null);  % emp vs null smoothness

% vis
figure; hold on;
scatter(ones(size(emp)), emp, 25, 'filled');
scatter(2*ones(size(null)), null, 25, 'filled');
plot([ones(numel(emp),1) 2*ones(numel(emp),1)]', [emp null]', 'k-');
xlim([0.5 2.5]); xticks([1 2]); xticklabels({'Empirical','Null'});
title("spatial autocorrelation robustness check");
subtitle(sprintf('paired t-test: t = %.2f, p = %.3g', tstat.tstat, p));
box off;
