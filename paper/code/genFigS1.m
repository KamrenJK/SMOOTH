function genFigS1
% ------------------------------------------------------------
% Supplementary Figure 1. Communication task and channel-level effect estimation.
%   b) Exemplar channel time series (HFB + communicative-action regressor)
%   c) Channel-level statistical maps (FDR q < 0.05)
% Panel a (task illustration) is a hand-drawn schematic, not produced by this script.
% ------------------------------------------------------------

%% config

clear; clc; close all;

% set paths -- edit smooth_setup.m once (repository root); nothing to change here.
% Bootstrap: put the repository root on the path so smooth_setup is findable
% regardless of the current working directory.
addpath(fullfile(fileparts(mfilename('fullpath')),'..','..'));
paths = smooth_setup('quiet');
paths.demoSource = fullfile(paths.demo, 'source.mat');
paths.tsFile     = fullfile(paths.data, 'ts.mat');
paths.figOut     = fullfile(paths.figs, 'S1');
if ~isfolder(paths.figOut), mkdir(paths.figOut); end

% set savefigs = false to render without writing PDFs to paths.figOut
savefigs = true;

% cortical mesh (only for visualization)
cortex = load_fsaverage_mesh(paths.freesurfer);

% read in data
load(paths.demoSource)  % path to example dataset
nsubs = numel(source);
nchans = sum(cellfun(@(s) numel(s.label), source));
fprintf('\nLoaded demo dataset: %d subjects, %d total electrodes.\n', nsubs, nchans);

%% organize channel-level data

CLstats = struct('sub',[],'stat',[],'fsxyz',[],'fssurfxyz',[],'dk',{{}});

for i = 1:numel(source)
    CLstats.sub       = vertcat(CLstats.sub,i.*ones(length(source{i}.label),1));
    CLstats.stat      = vertcat(CLstats.stat,source{i}.stat);
    CLstats.fsxyz     = vertcat(CLstats.fsxyz,source{i}.elec.chanpos);

    % correct jitter
    D = pdist2(cortex.pos,source{i}.elec.chanpos);
    [~, idx] = min(D,[],1);
    fssurfxyz = cortex.pos(idx,:);
    CLstats.fssurfxyz = vertcat(CLstats.fssurfxyz,fssurfxyz);

    % ID channels
    [roi_name, ~, ~, ~] = fsavg_coord_to_dk(fssurfxyz, paths.freesurfer);
    CLstats.dk        = vertcat(CLstats.dk,roi_name);
end

%% visualize channel-level stats (panel c)

% right-tail significance
ppos = 1 - normcdf(CLstats.stat);
[sig, ~, ~, ~]=fdr_bh(ppos);  % FDR correct

% colormap (orange/neutral/blue)
blue     = slanCM('Blues');
orange   = slanCM('Oranges');
berrymap = [orange(180,:); 0 0 0; blue(180,:)];

% plot
fig = figure('units','normalized','outerposition',[0 0 1 1],'visible','on','color',[1 1 1]); hold on;
ft_plot_mesh(cortex, 'vertexcolor', 'curv');
ft_plot_cloud(CLstats.fssurfxyz((sig & CLstats.stat>0),:), 1.*ones(sum(sig & CLstats.stat>0),1), 'cloudtype', 'surf', ...
  'radius', 2.5, 'scalerad', 'no', 'colormap', berrymap, 'clim', [-1 1]);   % sig
ft_plot_cloud(CLstats.fssurfxyz((~sig | CLstats.stat<0),:), zeros(sum(~sig | CLstats.stat<0),1), 'cloudtype', 'surf', ...
  'radius', 1, 'scalerad', 'no', 'colormap', berrymap, 'clim', [-1 1]);   % null

views    = [90 0; -90 0; 0 90; 180 -90; 180 -15; 0 0];
angle    = {'right', 'left', 'dorsal', 'ventral', 'rostral', 'caudal'};
for v = 1:size(views,1)
    view(views(v,:))
    if savefigs, exportgraphics(fig, [paths.figOut filesep 'berries_' angle{v} '.pdf'], 'contenttype', 'image', 'resolution', 450); end
end
close(fig)

%% exemplar channel timeseries (panel b)
load(paths.tsFile)
f = figure('color',[1 1 1],'Position',[100 300 1500 415]);
plot(ts.t,ts.y);
hold on; plot(ts.t,ts.r);
xlim([0 300])
xticks(0:150:300)
ylim([-1 2])
yticks(-1:2)
xlabel('time (s)')
ylabel('zHGB')
legend({'zHGB','mentalizing'})
if savefigs, exportgraphics(f, [paths.figOut filesep 'timeseries.pdf'], 'contenttype', 'vector'); end
close(f)

end

%% subfunctions
% -------------------------------------------------------------------------
function cortex = load_fsaverage_mesh(fshome)
% load cortex mesh given fs path

lh = ft_read_headshape(fullfile(fshome,'subjects','fsaverage','surf','lh.pial'));
rh = ft_read_headshape(fullfile(fshome,'subjects','fsaverage','surf','rh.pial'));

cortex      = [];
cortex.pos  = [lh.pos; rh.pos];
cortex.tri  = [lh.tri; rh.tri + size(lh.pos,1)];
cortex.unit = lh.unit;

if isfield(lh,'curv') && isfield(rh,'curv')
    cortex.curv = [lh.curv; rh.curv];
else
    cortex.curv = zeros(size(cortex.pos,1),1);
end
end
