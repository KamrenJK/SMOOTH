function genFigS2
% -------------------------------------------------------------------------
% Supplementary Figure 2. Participant-specific spatial basis derived from intracranial electrode
% geometry.
%   a) Exemplar subject fsavg-projected coverage (panel a: "channel locations ... subdural grid")
%   b) Pairwise distance matrix (PDM) and normalized graph Laplacian
%   c) Selected eigenmodes plotted on registered spherical coordinates
%   d) Spatial autocorrelation of eigenmodes, modal coefficients, group-average modal power spectrum
%
% Only the fsavg-projected coverage view is rendered for panel a (the rest of the figure operates in
% registered/spherical space, so this is the consistent choice).
% -------------------------------------------------------------------------

%% config

clear; clc; close all;

% set paths -- edit smooth_setup.m once (repository root); nothing to change here.
% Bootstrap: put the repository root on the path so smooth_setup is findable
% regardless of the current working directory.
addpath(fullfile(fileparts(mfilename('fullpath')),'..','..'));
paths = smooth_setup('quiet');
paths.demoSource = fullfile(paths.demo, 'source.mat');
paths.figData    = paths.data;
paths.figOut     = fullfile(paths.figs, 'S2');
if ~isfolder(paths.figOut), mkdir(paths.figOut); end

% set savefigs = false to render without writing PDFs to paths.figOut
savefigs = true;

% cortical mesh (for visualization)
cortex    = load_fsaverage_mesh(paths.freesurfer);
cortex_lh = ft_read_headshape([paths.freesurfer filesep 'subjects' filesep 'fsaverage' filesep 'surf' filesep 'lh.pial']);
cortex_rh = ft_read_headshape([paths.freesurfer filesep 'subjects' filesep 'fsaverage' filesep 'surf' filesep 'rh.pial']);

% read in data
load(paths.demoSource)  % path to example dataset
nsubs = numel(source);
nchans = sum(cellfun(@(s) numel(s.label), source));
fprintf('\nLoaded demo dataset: %d subjects, %d total electrodes.\n', nsubs, nchans);

% demo sub
exemplarID = 10;
%% a) plot coverage (fsavg-projected only; see header)

% plot projected fsavg
c = source{exemplarID}.elec.chanpos;
f = figure('color',[1 1 1]); hold on
ft_plot_mesh(cortex,'vertexcolor','curv','facealpha',0.25);
scatter3(c(:,1),c(:,2),c(:,3),45,'filled','MarkerFaceColor',[0 0 0])
title('fsavg space + cortical projection')
if savefigs, exportgraphics(f, [paths.figOut filesep 'fsavg.pdf'], 'contenttype', 'image', 'resolution', 450); end
close(f)

%% run SMOOTH on exemplar

cfg  = default_smooth_cfg(paths.freesurfer);
cfg.graphsigma = 8;
stat = SMOOTHsub(cfg, source{exemplarID});

%% 3) visualize eigenmodes

% config
nL        = length(stat.data.elec.idxL);
Q         = stat.data.modes(nL+1:end,nL+1:end); % RH eigenmodes
R         = 100;
mrkrsz    = 144;
pos       = stat.data.elec.Rh_sph;  % channel pos on reg sphere
modes2plot = [1 2 3 5 10 15 25 50 75 100 125];
[x, y, z] = sphere(15);  % unit sphere
POV       = [140.5 42.5];
% plot sampled modes
for q = modes2plot
    f = figure('color',[1 1 1],'units','centimeters','Position',[0 0 15 15]);
    surf(R.*x,R.*y,R.*z,'FaceAlpha',0.0,'LineWidth',0.375);  % plot orb
    hold on; scatter3(pos(:,1),pos(:,2),pos(:,3),mrkrsz,Q(:,q),'filled')  % plot channels
    axis equal; grid off; axis off
    clim([min(Q(:,q)) max(Q(:,q))]); xlim([-100 100]); ylim([-100 100]); zlim([-100 100])
    colormap(viridis)
    title(num2str(q))
    view(POV)
    if savefigs, exportgraphics(f, [paths.figOut filesep 'mode_' num2str(q) '.pdf'], 'contenttype', 'vector'); end
    close(f)
end

% plot coefficients
f = figure('color',[1 1 1]);
plot(stat.data.modalcoeffs(nL+1:end),'linewidth',1.5)
xticks([0 50 100 150]); yticks(-30:10:40); ylim([-30 40])
xlabel('mode'); ylabel('coefficient'); grid on
if savefigs, exportgraphics(f, [paths.figOut filesep 'modalcoeffs.pdf'], 'contenttype', 'vector'); end
close(f)

%% 4) modal decomp

% geodesic distance matrix
U     = pos./sqrt(sum(pos.^2, 2));  % unit vectors
C     = U * U.';                    % cos(theta)                 
C     = max(-1, min(1, C));
theta = acos(C);            
D     = 100 * theta;
sigma = 20;
W     = exp(-(D.^2) ./ (2*(sigma^2)));
n     = length(W);
W(1:n+1:end) = 0; % prevent self-loops
d     = sum(W,2); d(d==0) = 1;
Dm    = spdiags(1./sqrt(d), 0, n, n);
L     = speye(n) - (Dm * W * Dm);
L     = (L + L')/2;  % normalized laplacian

% visualize PDM
f = figure('color',[1 1 1]);
imagesc(D); colormap(viridis)
xlabel('chan'); ylabel('chan')
cb = colorbar;
cb.Label.String = 'geodesic distance (mm)';
cb.Label.FontSize = 16;
if savefigs, exportgraphics(f, [paths.figOut filesep 'PDM.pdf'], 'contenttype', 'vector'); end
close(f)

% visualize laplacian
f = figure('color',[1 1 1]);
imagesc(L); colormap(viridis)
xlabel('chan'); ylabel('chan')
cb = colorbar;
cb.Label.String = 'normalized Laplacian (dimensionless)';
cb.Label.FontSize = 16;
clim([-0.5 0.5])
if savefigs, exportgraphics(f, [paths.figOut filesep 'L.pdf'], 'contenttype', 'vector'); end
close(f)

%% 5) compute SA
modeSA = nan(length(Q),1);
for q = 1:length(Q)
    modeSA(q) = spatialautocorr(Q(:,q),D);
end

% visualize
f = figure('color',[1 1 1]);
scatter(1:length(modeSA),modeSA,72,'x','linewidth',1.5)
xticks([0 50 100 150]); yticks(-0.2:0.2:0.8); ylim([-0.2 0.8])
xlabel('mode'); ylabel('SA'); grid on
if savefigs, exportgraphics(f, [paths.figOut filesep 'modalSA.pdf'], 'contenttype', 'vector'); end
close(f)

%% 6) modal power (across subjects)

% run alg
LH = nan(numel(source),1000);
RH = nan(numel(source),1000);

for i = 1:numel(source)

    % run decomp
    cfg  = default_smooth_cfg(paths.freesurfer);
    cfg.graphsigma = 8;
    stat = SMOOTHsub(cfg, source{i});

    % compute power
    nL      = length(stat.data.elec.idxL);
    powLH   = smooth(stat.data.modalcoeffs(1:nL).^2);  % p = weight.^2
    powRH   = smooth(stat.data.modalcoeffs(nL+1:end).^2);
    LH(i,:) = [powLH' nan(1,(1000-length(powLH)))];
    RH(i,:) = [powRH' nan(1,(1000-length(powRH)))];
end

% plot powerspectrum
keepLH = sum(~isnan(LH),1)>3;
keepRH = sum(~isnan(RH),1)>3;
meanLH = nanmean(LH(:,keepLH),1);
meanRH = nanmean(RH(:,keepRH),1);
f = figure('color',[1 1 1]);
plot(1:length(meanLH),meanLH)
hold on; plot(1:length(meanRH),meanRH)
set(gca,'XScale','log','YScale','log')
xlabel('mode'); ylabel('pow'); box off
legend({'LH','RH'})
if savefigs, exportgraphics(f, [paths.figOut filesep 'powerspectrum.pdf'], 'contenttype', 'vector'); end
close(f)

end
%% subfunctions
function cfg = default_smooth_cfg(fshome)
  cfg = [];
  cfg.fshome = fshome;
  cfg.keepmaps = 'yes';
  cfg.keepsurrogates = 'yes';
  cfg.normalize = 'no';
  cfg.numrandomization = 0;
  cfg.kernelwidth = 12;
  cfg.rankrescale = 'exact';
  cfg.smooth = 'sphere';
  cfg.kernel = 'gaussian';
end

