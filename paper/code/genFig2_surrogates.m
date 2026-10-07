function genFig2_surrogates
% -------------------------------------------------------------------------
% Figure 2, panels c-f: surrogate generation procedure.
%   c) Exemplar empirical & surrogate cortical + channel-level maps
%   d) Spectral resampling and reconstruction (mode coefficients before/after sign flip)
%   e) Semivariogram, empirical vs. surrogate (mean +/- SD across 1,000 surrogates)
%   f) Group-level empirical vs. surrogate Moran's I (paired half-violin; the distributions do not
%      differ -- this is the correct plot for that claim, not a scatter -- see genFig2.m for panel b
%      and genFig2_seedcorr.m for panel a)
% -------------------------------------------------------------------------

%% config

clear; clc; close all;

% set paths -- edit smooth_setup.m once (repository root); nothing to change here.
% Bootstrap: put the repository root on the path so smooth_setup is findable
% regardless of the current working directory.
addpath(fullfile(fileparts(mfilename('fullpath')),'..','..'));
paths = smooth_setup('quiet');
paths.demoSource = fullfile(paths.demo, 'source.mat');
paths.figOut     = fullfile(paths.figs, '2_surrogates');
if ~isfolder(paths.figOut), mkdir(paths.figOut); end

% set savefigs = false to render without writing PDFs to paths.figOut
savefigs = true;

% cortical mesh (for visualization)
cortex = load_fsaverage_mesh(paths.freesurfer);
cortex_lh   = ft_read_headshape([paths.freesurfer filesep 'subjects' filesep 'fsaverage' filesep 'surf' filesep 'lh.pial']);
cortex_rh   = ft_read_headshape([paths.freesurfer filesep 'subjects' filesep 'fsaverage' filesep 'surf' filesep 'rh.pial']);

% read in data
load(paths.demoSource)  % path to example dataset
nsubs = numel(source);
nchans = sum(cellfun(@(s) numel(s.label), source));
fprintf('\nLoaded demo dataset: %d subjects, %d total electrodes.\n', nsubs, nchans);

%% generate surrogates

% run alg
cfg                  = [];
cfg.fshome           = paths.freesurfer;
cfg.keepmaps         = 'yes';
cfg.keepsurrogates   = 'yes';
cfg.subsamplesurr    = 'no';
cfg.normalize        = 'no';
cfg.numrandomization = 1000;
cfg.kernelwidth      = 12;
cfg.graphsigma       = 8;
cfg.rankrescale      = 'exact';
cfg.smooth           = 'sphere';
cfg.kernel           = 'gaussian';
stat = SMOOTHsub(cfg, source{10});

%% plot empirical cortical map

f  = figure('color',[1 1 1],'units','centimeters','Position',[0 0 15 15]); hold on

% colormap
blue     = slanCM('Blues');
orange   = slanCM('Oranges');

% plot channels
ft_plot_mesh(cortex, 'vertexcolor', 'curv');
hold on; ft_plot_mesh(cortex, 'vertexcolor', stat.maps);
clim([-0.04 0.04])
colormap(blue)
view([90 0])

if savefigs, exportgraphics(f,[paths.figOut filesep 'empiricalmap.pdf'],'ContentType','auto','Resolution',600); end
close(f)

%% plot surrogate cortical map

f  = figure('color',[1 1 1],'units','centimeters','Position',[0 0 15 15]); hold on

% plot channels
ft_plot_mesh(cortex, 'vertexcolor', 'curv');
hold on; ft_plot_mesh(cortex, 'vertexcolor', stat.surrogate_subject(:,1));
clim([-0.04 0.04])
colormap(orange)
view([90 0])

if savefigs, exportgraphics(f,[paths.figOut filesep 'surrmap.pdf'],'ContentType','auto','Resolution',600); end
close(f)

%% modal power spectrum (RH only)

nL = sum(stat.data.elec.hemi);

% real
f = figure('color',[1 1 1],'units','centimeters','Position',[0 0 20 15]);
plot(stat.data.modalcoeffs(nL+1:end))
xlabel('mode')
ylabel('modal coefficient')
xticks(0:40:160)
ylim([-30 30])
yticks(-30:10:30)
box off
grid on
if savefigs, exportgraphics(f,[paths.figOut filesep 'empcoeffs.pdf'],'ContentType','vector'); end
close(f)

% signflipped 
f = figure('color',[1 1 1],'units','centimeters','Position',[0 0 20 15]);
y = stat.data.modalcoeffs.*stat.data.sflipmat(:,1);  % signflip and reconstruct
plot(y(nL+1:end))
xlabel('mode')
ylabel('modal coefficient')
xticks(0:40:160)
ylim([-30 30])
yticks(-30:10:30)
box off
grid on
if savefigs, exportgraphics(f,[paths.figOut filesep 'surrcoeffs.pdf'],'ContentType','vector'); end
close(f)

%% plot original y

% reconstruct
y = source{10}.stat;

% channel coords (snap to mesh)
chanXYZ = source{10}.elec.chanpos;
D = pdist2(cortex.pos,chanXYZ);
[~, idx] = min(D,[],1);
fssurfxyz = cortex.pos(idx,:);
fssurfxyz(stat.data.elec.hemi,:) = [];

% plot (RH)
f  = figure('color',[1 1 1],'units','centimeters','Position',[0 0 15 15]); hold on
ft_plot_mesh(cortex_rh,'vertexcolor','curv')
ft_plot_cloud(fssurfxyz, y(~stat.data.elec.hemi), 'cloudtype', 'surf', ...
  'radius', 2, 'scalerad', 'no', 'colormap', blue, 'clim', [-10 10]); 
view([90 0])

if savefigs, exportgraphics(f,[paths.figOut filesep 'empberries.pdf'],'ContentType','auto','Resolution',600); end
close(f)

%% plot reconstructed y

% reconstruct
% y_recon = stat.data.modes*(stat.data.modalcoeffs);  % SANITY CHECK
y_recon = stat.data.modes*(stat.data.modalcoeffs.*stat.data.sflipmat(:,1));
y_recon = y_recon(nL+1:end);

% plot
f  = figure('color',[1 1 1],'units','centimeters','Position',[0 0 15 15]); hold on
ft_plot_mesh(cortex_rh,'vertexcolor','curv')
ft_plot_cloud(fssurfxyz, y_recon, 'cloudtype', 'surf', ...
  'radius', 2, 'scalerad', 'no', 'colormap', orange, 'clim', [-10 10]); 
view([90 0])

if savefigs, exportgraphics(f,[paths.figOut filesep 'surrberries.pdf'],'ContentType','auto','Resolution',600); end
close(f)

%% compute spatial stats on surrogates

% empirical
surfmap = stat.maps(end/2:end);
cvg     = isnan(surfmap);
% surfmap(cvg) = [];

% downsample mesh for easy variogram calculation
reg_rh            = ft_read_headshape(fullfile(paths.freesurfer,'subjects','fsaverage','surf','rh.sphere.reg'));
mesh_patch        = struct('vertices', reg_rh.pos, 'faces', reg_rh.tri);  % LH only
reduction_factor  = 0.025;
[tri_ds, pos_ds]  = reducepatch(mesh_patch, reduction_factor);
reg_ds.tri        = tri_ds;
reg_ds.pos        = pos_ds;
D_ds              = pdist2(reg_ds.pos, reg_ds.pos);
Dmap              = pdist2(reg_rh.pos, reg_ds.pos);
[~, dsidx]        = min(Dmap, [], 1);   % project to downsampled mesh
dssurfmap         = surfmap(dsidx);
dscvg             = isnan(dssurfmap);
D_ds(dscvg,:)     = [];
D_ds(:,dscvg)     = [];
surfmap(cvg)      = [];

%% variogram of permuted data
% config
vi      = 5;      
vf      = 30;     
nperms  = 1000;    

% save features 
semivar       = [];
morans        = nan(nperms,1);
mapcorr       = nan(nperms,1);

for i = 1:nperms + 1
    if i <= nperms  % surrogate
        surfmapi = stat.surrogate_subject(end/2:end,i);
    else  % empirical
        surfmapi = stat.maps(end/2:end);
    end
    dssurfmapi = surfmapi(dsidx);  % downsample
    dssurfmapi(dscvg) = [];        % trim nans
    surfmapi(cvg) = [];            % trim nans

    % variogram
    x_rh          = variogram(dssurfmapi, vi, vf, D_ds);
    semivar       = [semivar, x_rh];

    % SA
    morans(i) = spatialautocorr(dssurfmapi, D_ds);

    % corr
    mapcorr(i) = corr(surfmap,surfmapi);
end

%% plot spatial stats

mu    = nanmean(semivar(:,1:end-1),2);
sigma = nanstd(semivar(:,1:end-1),0,2);
lo    = (mu-sigma)';
hi    = (mu+sigma)';
x     = vi:(vf-vi)/20:vf;

% variogram 
f = figure('color',[1 1 1]);
scatter(x,semivar(:,end),'LineWidth',1.5); hold on;
fill([x fliplr(x)], [lo, fliplr(hi)],[0.7 0.7 0.7],'FaceAlpha', 0.25, 'EdgeAlpha',0.2)
legend({'empirical','surrogate'},'Location','northwest')
xlabel('spatial distance (mm)'); ylabel('variance')
box off
grid on
if savefigs, exportgraphics(f,[paths.figOut filesep 'var.pdf'],'ContentType','vector'); end
close(f)

% variogram 
f = figure('color',[1 1 1]);
histogram(mapcorr(1:end-1));
xlabel('corr'); ylabel('count'); title('empirical to surrogate correlations')
xlim([-1 1])
box off
grid on
if savefigs, exportgraphics(f,[paths.figOut filesep 'corr.pdf'],'ContentType','vector'); end
close(f)

%% compute correlations

% permute and save SA
cfg                  = [];
cfg.fshome           = paths.freesurfer;
cfg.keepmaps         = 'yes';
cfg.keepsurrogates   = 'yes';
cfg.subsamplesurr    = 'no';
cfg.normalize        = 'no';
cfg.numrandomization = 500;
cfg.kernelwidth      = 12;
cfg.rankrescale      = 'exact';
cfg.smooth           = 'sphere';
cfg.kernel           = 'gaussian';
stat = SMOOTHstat(cfg, source{:});
d = SMOOTHdiag(stat);

%% visualize group-level emp-surr SA parity

% compare empirical and avg surrogate SA for all subs
emp  = d.subject.spatial.BOTH.moranI.per_subject.emp;
null = d.subject.spatial.BOTH.moranI.per_subject.null_mu;
pair = isfinite(emp) & isfinite(null);
e    = emp(pair);
n    = null(pair);
[~,p,~,STATS] = ttest(e, n); 

% formatting
halfwidth = 0.25;         
colEmp    = [0.2 0.5 0.9];
colNull   = [0.85 0.4 0.2];

% plot
f = figure('color',[1 1 1]); hold on
halfviolin(1, e, 'left', halfwidth, colEmp, 0.25)
scatter(1.25.*ones(size(e)), e, 144, 'v', 'MarkerFaceColor', colEmp, 'MarkerEdgeAlpha', 0)
scatter(1.75.*ones(size(n)), n, 144, '^', 'MarkerFaceColor', colNull, 'MarkerEdgeAlpha', 0)
for i = 1:numel(e)  % connect
    plot([1.25 1.75], [e(i) n(i)], 'Color', [0 0 0 0.3])
end
halfviolin(2, n, 'right', halfwidth, colNull, 0.25)
xlim([0.5 2.5]); xticks([1 2]); xticklabels({'empirical','null'})
ylabel('Moran''s I'); ylim([0.98 1.04]); title('emp vs null Moran''s I')
subtitle(sprintf('t(%d) = %.2f, p = %.3g', STATS.df, STATS.tstat, p))
box off; grid on
axis square
if savefigs, exportgraphics(f,[paths.figOut filesep 'grouplevel.pdf'],'ContentType','vector'); end
close(f)

%% subfunctions
function B = surfacegauss(A,D,Idx,sigma)
% A = original map, D = Distances of nearest points to query points, Idx =
% Indices of nearest points (from rangesearch for mem efficiency)
B = nan(size(A));
for v = 1:length(D)
    w = exp(-(D{v}.^2)/(2*sigma^2));
    w = w./nansum(w);
    y = nansum(w.*A(Idx{v})');
    B(v) = y;
end

function I = spatialautocorr(A, D)
% compute spatial autocorrelation (Moran's I)
% A: n x 1 spatial map
% D: n x n pairwise distance matrix
D(D==0) = Inf;
D(D>100) = Inf; % theshold large distances
Dinv    = 1./D;
Dinv    = Dinv./max(Dinv,[],1);  % weight by inverse distance
A = A - mean(A,'all');  % center on zero
I = (numel(A)/sum(Dinv,'all'))*(A(:)'*Dinv*A(:))/(A(:)'*A(:));

function halfviolin(xc, data, side, halfwidth, faceColor, alphaVal)
% half violin plots
data = data(isfinite(data));
if numel(data) < 5, return; end
[pdf,y] = ksdensity(data,'NumPoints',256);      
pdf = pdf ./ max(pdf); % normalize
w = halfwidth * pdf; % scale
switch lower(side)
    case 'left',  x = xc - w;
    case 'right', x = xc + w;
    otherwise, error('side must be ''left'' or ''right''.');
end
xpoly = [x, repmat(xc,1,numel(x))];
ypoly = [y, fliplr(y)];
patch(xpoly, ypoly, faceColor, 'FaceAlpha', alphaVal, 'EdgeColor', 'none');
