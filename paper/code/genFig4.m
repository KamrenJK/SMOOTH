function genFig4
% -------------------------------------------------------------------------
% Figure 4. Surrogate generation procedure
%   1) Exemplar empirical & surrogate maps
%   2) Exemplar empirical & surrogate channel-level plot
%   3) Spectral resampling and reconstruction
%   4) Group-level surrogate spatial structure
% -------------------------------------------------------------------------

%% config

clear; clc; close all;

% set paths
paths.smooth      = '/Users/kamrenkhan/Desktop/Research/Mu Lab/Ecog/analyses/Kamren/smooth-pub';  % UPDATE
paths.fieldtrip   = '/Users/kamrenkhan/fieldtrip-master';                                         % UPDATE
paths.freesurfer  = '/Applications/freesurfer/7.3.2';                                             % UPDATE
paths.demoSource  = fullfile(paths.smooth, 'demo', 'source.mat');
paths.figOut      = fullfile(paths.smooth, 'paper', 'figs', '4');

% add paths
addpath(genpath(paths.smooth));
addpath(paths.fieldtrip);
ft_defaults;

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
cfg        = [];
cfg.fshome = paths.freesurfer;
cfg.keepmaps = 'yes';
cfg.keepsurrogates = 'yes';
cfg.normalize = 'no';
cfg.numrandomization = 1000;
cfg.kernelwidth = 30;
cfg.rankrescale = 'exact';
cfg.smooth = 'sphere';
cfg.kernel = 'gaussian';
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

% exportgraphics(f,[paths.figOut filesep 'empiricalmap.pdf'],'ContentType','auto','Resolution',600)
% close(f)

%% plot surrogate cortical map

f  = figure('color',[1 1 1],'units','centimeters','Position',[0 0 15 15]); hold on

% plot channels
ft_plot_mesh(cortex, 'vertexcolor', 'curv');
hold on; ft_plot_mesh(cortex, 'vertexcolor', stat.surrogate_subject(:,7));
clim([-0.04 0.04])
colormap(orange)
view([90 0])

% exportgraphics(f,[paths.figOut filesep 'surrmap.pdf'],'ContentType','auto','Resolution',600)
% close(f)


%% modal power spectrum

% real
f = figure('color',[1 1 1],'units','centimeters','Position',[0 0 20 15]);
plot(stat.data.modalcoeffs)
xlabel('mode')
ylabel('modal coefficient')
xticks(0:40:160)
ylim([-30 30])
yticks(-30:10:30)
box off
grid on
% exportgraphics(f,[paths.figOut filesep 'empcoeffs.pdf'],'ContentType','vector')
% close(f)

% signflipped 
f = figure('color',[1 1 1],'units','centimeters','Position',[0 0 20 15]);
plot(stat.data.modalcoeffs.*stat.data.sflipmat(:,7))
xlabel('mode')
ylabel('modal coefficient')
xticks(0:40:160)
ylim([-30 30])
yticks(-30:10:30)
box off
grid on
% exportgraphics(f,[paths.figOut filesep 'surrcoeffs.pdf'],'ContentType','vector')
% close(f)

%% plot original y

% reconstruct
y = source{10}.stat;

% channel coords (snap to mesh)
chanXYZ = source{10}.elec.chanpos;
D = pdist2(cortex.pos,chanXYZ);
[~, idx] = min(D,[],1);
fssurfxyz = cortex.pos(idx,:);

% plot
f  = figure('color',[1 1 1],'units','centimeters','Position',[0 0 15 15]); hold on
ft_plot_mesh(cortex_rh,'vertexcolor','curv')
ft_plot_cloud(fssurfxyz, y, 'cloudtype', 'surf', ...
  'radius', 2, 'scalerad', 'no', 'colormap', blue, 'clim', [-15 15]); 
view([90 0])

% exportgraphics(f,[paths.figOut filesep 'empberries.pdf'],'ContentType','auto','Resolution',600)
% close(f)

%% plot reconstructed y

% reconstruct
y_recon = stat.data.modes*(stat.data.modalcoeffs.*stat.data.sflipmat(:,7));

% plot
f  = figure('color',[1 1 1],'units','centimeters','Position',[0 0 15 15]); hold on
ft_plot_mesh(cortex_rh,'vertexcolor','curv')
ft_plot_cloud(fssurfxyz, y_recon, 'cloudtype', 'surf', ...
  'radius', 2, 'scalerad', 'no', 'colormap', orange, 'clim', [-15 15]); 
view([90 0])

% exportgraphics(f,[paths.figOut filesep 'surrberries.pdf'],'ContentType','auto','Resolution',600)
% close(f)

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
surfmap(cvg) = [];

%%
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
legend({'empirical','surrogate'})
xlabel('spatial distance (mm)'); ylabel('variance')
box off
grid on
% exportgraphics(f,[paths.figOut filesep 'var.pdf'],'ContentType','vector')
% close(f)

% variogram 
f = figure('color',[1 1 1]);
histogram(mapcorr(1:end-1));
xlabel('corr'); ylabel('count'); title('empirical to surrogate correlations')
xlim([-1 1])
box off
grid on
% exportgraphics(f,[paths.figOut filesep 'corr.pdf'],'ContentType','vector')
% close(f)

%% plot group level correlations

% permute and save SA
cfg                  = [];
cfg.fshome           = paths.freesurfer;
cfg.keepmaps         = 'yes';
cfg.keepsurrogates   = 'yes';
cfg.normalize        = 'no';
cfg.numrandomization = 500;
cfg.kernelwidth      = 30;
cfg.rankrescale      = 'exact';
cfg.smooth           = 'sphere';
cfg.kernel           = 'gaussian';
stat = SMOOTHstat(cfg, source{:});
d = SMOOTHdiag(stat);
%%
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
figure; hold on
halfviolin(1, e, 'left', halfwidth, colEmp, 0.25)
scatter(ones(size(e)), e, 36, 'o', 'MarkerFaceColor', colEmp, 'MarkerEdgeColor', 'k')
scatter(2.*ones(size(n)), n, 36, 'o', 'MarkerFaceColor', colNull, 'MarkerEdgeColor', 'k')
for i = 1:numel(e)  % connect
    plot([1 2], [e(i) n(i)], 'Color', [0 0 0 0.3])
end
halfviolin(2, n, 'right', halfwidth, colNull, 0.25)
xlim([0.5 2.5]); xticks([1 2]); xticklabels({'empirical','null'})
ylabel('Moran''s I'); ylim([0.98 1.04]); title('emp vs null Moran''s I')
subtitle(sprintf('t(%d) = %.2f, p = %.3g', STATS.df, STATS.tstat, p))
box on


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
