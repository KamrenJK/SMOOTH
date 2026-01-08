function genFig5
% -------------------------------------------------------------------------
% Figure 5. Group-level results
%   1) Group level empirical, SMOOTH, RSF maps
%   2) SMOOTH & RSF null distributions
%   3) Above threshold clusters
% -------------------------------------------------------------------------

%% config

clear; clc; close all;

% set paths
paths.smooth      = '/LOCAL/PATH/TO/SMOOTH';        % UPDATE
paths.fieldtrip   = '/LOCAL/PATH/TO/FIELDTRIP';     % UPDATE
paths.freesurfer  = '/LOCAL/PATH/TO/FREESURFER';    % UPDATE
paths.demoSource  = fullfile(paths.smooth, 'demo', 'source.mat');
paths.figOut      = fullfile(paths.smooth, 'paper', 'figs', '5');

% add paths
addpath(genpath(paths.smooth));
addpath(paths.fieldtrip);
ft_defaults;

% cortical mesh (for visualization)
cortex      = load_fsaverage_mesh(paths.freesurfer);
cortex_lh   = ft_read_headshape([paths.freesurfer filesep 'subjects' filesep 'fsaverage' filesep 'surf' filesep 'lh.pial']);
cortex_rh   = ft_read_headshape([paths.freesurfer filesep 'subjects' filesep 'fsaverage' filesep 'surf' filesep 'rh.pial']);

% read in data
load(paths.demoSource)  % path to example dataset
nsubs = numel(source);
nchans = sum(cellfun(@(s) numel(s.label), source));
fprintf('\nLoaded demo dataset: %d subjects, %d total electrodes.\n', nsubs, nchans);

%% run stats

% run SMOOTH
cfg                   = [];
cfg.fshome           = paths.freesurfer;
cfg.keepmaps         = 'no';
cfg.keepsurrogates   = 'no';
cfg.normalize        = 'no';
cfg.numrandomization = 1000;
cfg.kernelwidth      = 30;
cfg.rankrescale      = 'exact';
cfg.smooth           = 'sphere';
cfg.kernel           = 'gaussian';
stat = SMOOTHstat(cfg, source{:});
statRSF = SMOOTHdummy(cfg, source{:});  % repeat with RSF

% exemplars
cfg.numrandomization = 1;
cfg.keepmaps         = 'yes';
cfg.keepsurrogates   = 'yes';
exampleSTAT = SMOOTHstat(cfg, source{:});
exampleRSF = SMOOTHdummy(cfg, source{:});

%% visualize maps

% empirical
f = figure('color',[1 1 1]);
ft_plot_mesh(cortex, 'vertexcolor', 'curv');
hold on; ft_plot_mesh(cortex, 'vertexcolor', stat.stat);
clim([-3 3])
colormap(slanCM('Blues'));
view([90 0])
% exportgraphics(f, [figpath filesep 'emp_rh.pdf'], 'contenttype', 'image', 'resolution', 450)
view([-90 0])
% exportgraphics(f, [figpath filesep 'emp_lh.pdf'], 'contenttype', 'image', 'resolution', 450)
% close(f)

% SMOOTH
f = figure('color',[1 1 1]);
ft_plot_mesh(cortex, 'vertexcolor', 'curv');
hold on; ft_plot_mesh(cortex, 'vertexcolor', exampleSTAT.surrogate_group);
clim([-3 3])
colormap(slanCM('Oranges'));
view([90 0])
% exportgraphics(f, [figpath filesep 'smooth_rh.pdf'], 'contenttype', 'image', 'resolution', 450)
view([-90 0])
% exportgraphics(f, [figpath filesep 'smooth_lh.pdf'], 'contenttype', 'image', 'resolution', 450)
% close(f)

% RSF
f = figure('color',[1 1 1]);
ft_plot_mesh(cortex, 'vertexcolor', 'curv');
hold on; ft_plot_mesh(cortex, 'vertexcolor', exampleRSF.surrogate_group);
clim([-3 3])
colormap(slanCM('Greens'));
view([90 0])
% exportgraphics(f, [figpath filesep 'rsf_rh.pdf'], 'contenttype', 'image', 'resolution', 450)
view([-90 0])
% exportgraphics(f, [figpath filesep 'rsf_lh.pdf'], 'contenttype', 'image', 'resolution', 450)
% close(f)

%% visualize null distributions

f = figure('color',[1 1 1]);
histogram(log(stat.posdistribution),'EdgeAlpha',0)
hold on; histogram(log(statRSF.posdistribution),'EdgeAlpha',0)
xline(log(stat.posclusters(1).clusterstat))
xline(log(stat.posclusters(2).clusterstat))
xline(log(stat.posclusters(3).clusterstat))
legend({'SMOOTH','RSF'})
box off; grid on
xlabel('log(cluster energy)')
ylabel('count')
yticks(0:40:160)
% exportgraphics(f, [figpath filesep 'nulldist.pdf'], 'contenttype', 'vector')
% close(f)

%% plot sig clusters

% inflated map
inflated_left     = ft_read_headshape('/Applications/freesurfer/7.3.2/subjects/fsaverage/surf/lh.inflated');
inflated_left.curv(inflated_left.curv<0)=-.3;
inflated_left.curv(inflated_left.curv>0)=-.1;
inflated_right    = ft_read_headshape('/Applications/freesurfer/7.3.2/subjects/fsaverage/surf/rh.inflated');
inflated_right.curv(inflated_right.curv<0)=-.3;
inflated_right.curv(inflated_right.curv>0)=-.1;
inflated.pos      = cat(1,[inflated_left.pos; inflated_right.pos]); % concatenate
inflated.tri      = cat(1,[inflated_left.tri; inflated_right.tri + size(inflated_left.pos,1)]);
inflated.unit     = inflated_left.unit;
inflated.curv     = cat(1,[inflated_left.curv; inflated_right.curv]);


% first sub-threshold cluster
y = stat.posclusterslabelmat;
y(y~=1) = nan;
f = figure('color',[1 1 1]);
ft_plot_mesh(inflated, 'vertexcolor', 'curv');
hold on; ft_plot_mesh(inflated, 'vertexcolor', y);
clim([-3 3])
colormap(slanCM('Blues'))
view([-90 0])
% exportgraphics(f, [figpath filesep 'clus1.pdf'], 'contenttype', 'auto','Resolution',600)
% close(f)

% second sub-threshold cluster
y = stat.posclusterslabelmat;
y(y~=2) = nan;
f = figure('color',[1 1 1]);
ft_plot_mesh(inflated, 'vertexcolor', 'curv');
hold on; ft_plot_mesh(inflated, 'vertexcolor', y);
clim([-3 3])
colormap(slanCM('Blues'))
view([90 0])
% exportgraphics(f, [figpath filesep 'clus2.pdf'], 'contenttype', 'auto','Resolution',600)
% close(f)

% third sub-threshold cluster
y = stat.posclusterslabelmat;
y(y~=3) = nan;
f = figure('color',[1 1 1]);
ft_plot_mesh(inflated, 'vertexcolor', 'curv');
hold on; ft_plot_mesh(inflated, 'vertexcolor', y);
clim([-3 3])
colormap(slanCM('Blues'))
view([-90 0])
% exportgraphics(f, [figpath filesep 'clus3.pdf'], 'contenttype', 'auto','Resolution',600)
% close(f)

%% save colormaps
% f = figure('color',[1 1 1]);
% axis off
% box off
% cb = colorbar;
% cb.Ticks = [];
% colormap(slanCM('Blues'))
% exportgraphics(f, [figpath filesep 'cmap_emp.pdf'], 'contenttype', 'vector')
% colormap(slanCM('Oranges'))
% exportgraphics(f, [figpath filesep 'cmap_smoo.pdf'], 'contenttype', 'vector')
% colormap(slanCM('Greens'))
% exportgraphics(f, [figpath filesep 'cmp_rsf.pdf'], 'contenttype', 'vector')
% close(f)

