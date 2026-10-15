function genFigS1
% ------------------------------------------------------------
% Supplementary Figure 1. Communication task and channel-level effect estimation.
%   b) Exemplar channel time series (HFB + communicative-action regressor)
%   c) Channel-level statistical maps (FDR q < 0.05)
%   coverage map: participants contributing to each vertex, lateral and medial views (the same
%      quantity as Figure 5a, read from genFig5's saved results)
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

%% participant coverage (coverage map panel)
% Number of participants whose smoothed map covers each vertex, as in Figure 5a (top), from
% genFig5_results.mat; vertices covered by fewer than three participants are left blank.
R    = load(fullfile(paths.data, 'genFig5_results.mat'), 'primary');
cvg  = double(R.primary.coverage(:)); cvg(cvg < 3) = NaN;
grn  = slanCM('Greens');
nL   = size(cortex.pos,1) / 2;                         % fsaverage: 163,842 vertices per hemisphere
hemi = {'lh', 1:nL, cortex.tri(all(cortex.tri <= nL, 2), :);
        'rh', nL+1:2*nL, cortex.tri(all(cortex.tri > nL, 2), :) - nL};
side = {'lateral', 'medial'}; vw = {[-90 0], [90 0]; [90 0], [-90 0]};   % rows lh, rh
for h = 1:2
    m = struct('pos', cortex.pos(hemi{h,2},:), 'tri', hemi{h,3}, 'curv', cortex.curv(hemi{h,2}));
    c = cvg(hemi{h,2}); kc = isfinite(c);                   % covered vertices, drawn as their own surface
    vm = zeros(size(m.pos,1),1); vm(kc) = 1:nnz(kc);
    sm = struct('pos', m.pos(kc,:), 'tri', vm(m.tri(all(kc(m.tri),2),:)));
    for k = 1:2
        figure('color','w','position',[100 100 560 500]); axes('Position',[0.02 0.04 0.96 0.9]);
        ft_plot_mesh(m, 'vertexcolor', 'curv', 'edgecolor', 'none'); hold on
        ft_plot_mesh(sm, 'vertexcolor', c(kc), 'edgecolor', 'none');
        colormap(grn); clim([0 20]); view(vw{h,k}); axis off
        if savefigs, savepanel(paths.figOut, ['coverage_' hemi{h,1} '_' side{k}]); else, close(gcf); end
    end
end
f = figure('color','w','position',[100 100 430 150]);
ax = axes('Parent', f, 'Position', [0.1 0.6 0.8 0.02]); axis(ax, 'off'); colormap(ax, grn); set(ax, 'CLim', [0 20]);
cb = colorbar(ax, 'southoutside'); cb.Limits = [0 20]; cb.Ticks = [0 20]; cb.Label.String = 'subjects';
if savefigs, exportgraphics(f, [paths.figOut filesep 'cbar_coverage.pdf'], 'contenttype', 'vector'); end
close(f)

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
function savepanel(outdir,name)
% transparent PNG by two-background compositing (exportgraphics drops the alpha channel here)
fw=fullfile(outdir,[name '_w.png']); fk=fullfile(outdir,[name '_k.png']);
set(gcf,'InvertHardcopy','off','PaperPositionMode','auto'); set(findall(gcf,'type','axes'),'Color','none');
set(gcf,'Color','w'); print(gcf,fw,'-dpng','-r300'); set(gcf,'Color','k'); print(gcf,fk,'-dpng','-r300'); close(gcf);
W=double(imread(fw)); K=double(imread(fk)); alpha=min(max(1-mean(W-K,3)/255,0),1); C=min(max(K./max(alpha,1e-3),0),255);
[rr,cc]=find(alpha>0.02); if ~isempty(rr), r=min(rr):max(rr); c2=min(cc):max(cc); C=C(r,c2,:); alpha=alpha(r,c2); end
imwrite(uint8(C),fullfile(outdir,[name '.png']),'Alpha',alpha); delete(fw); delete(fk);
end

% -------------------------------------------------------------------------
function cortex = load_fsaverage_mesh(fshome)
% load cortex mesh given fs path

lh = smooth_read_surf(fullfile(fshome,'subjects','fsaverage','surf','lh.pial'));
rh = smooth_read_surf(fullfile(fshome,'subjects','fsaverage','surf','rh.pial'));

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
