function genFig3
% -------------------------------------------------------------------------
% Figure 3. Spatial Autocorrelation of Intracranial EEG Signal
%   1) Exemplar subject cortical coverage views (native)
%   2) Exemplar LFP / LFB / HFB correlation heatmaps + matching colorbar
%   3) Group FWHM (mm) summary with bootstrap CIs
%   4) Exemplar distance–correlation scatter with fitted Gaussian curves
%% -------------------------------------------------------------------------

%% config

clear; clc; close all;

% set paths
paths.smooth      = '/LOCAL/PATH/TO/SMOOTH';        % UPDATE
paths.fieldtrip   = '/LOCAL/PATH/TO/FIELDTRIP';     % UPDATE
paths.freesurfer  = '/LOCAL/PATH/TO/FREESURFER';    % UPDATE
paths.corrSource  = fullfile(paths.smooth, 'paper', 'figs','3','corrsource.mat'); % this file also not created by this function or demo script. BG
paths.figOut      = fullfile(paths.smooth, 'paper', 'figs', '3');

% add paths
addpath(genpath(paths.smooth));
addpath(paths.fieldtrip);
ft_defaults;

% read in data
load(paths.corrSource) 
load(fullfile(paths.figOut,'natcortex.mat'))
sID = 13;  % exemplar subj
nsubs = numel(corrsource);
nchans = sum(cellfun(@(s) numel(s.label), corrsource));
fprintf('\nLoaded demo dataset: %d subjects, %d total electrodes.\n', nsubs, nchans);

%% visualize coverage

f  = figure('color',[1 1 1],'units','centimeters','Position',[0 0 15 15]); hold on

% plot channels
ft_plot_mesh(natcortex,'vertexcolor','curv','facealpha',0.25)
ft_plot_cloud(corrsource{sID}.elec.nativechanpos, ones(length(corrsource{sID}.elec.nativechanpos),1), 'cloudtype', 'surf', ...
  'radius', 2, 'scalerad', 'no', 'colormap', 'gray', 'clim', [1 2]); 
view([180 0])
% exportgraphics(f,[figpath filesep sID '_anterior.pdf'],'ContentType','auto','Resolution',600)
% view([0 90])
% exportgraphics(f,[figpath filesep sID '_dorsal.pdf'],'ContentType','auto','Resolution',600)
% close(f)

%% exemplar corr data

% visualize
f = figure('color',[1 1 1]);
imagesc(corrsource{sID}.LFPcorr)
clim([-0.25 0.5])
colormap(viridis)
xlabel('channel')
ylabel('channel')
title('LFP')
axis square
% exportgraphics(f,[figpath filesep sID '_LFPcorr.pdf'],'ContentType','auto','Resolution',600)
% close(f)

f = figure('color',[1 1 1]);
imagesc(corrsource{sID}.LFcorr)
clim([-0.25 0.5])
colormap(viridis)
xlabel('channel')
ylabel('channel')
title('LF')
axis square
% exportgraphics(f,[figpath filesep sID '_NBcorr.pdf'],'ContentType','auto','Resolution',600)
% close(f)

f  = figure('color',[1 1 1]);
imagesc(corrsource{sID}.HFBcorr)
clim([-0.25 0.5])
colormap(viridis)
xlabel('channel')
ylabel('channel')
title('HFB')
axis square
% exportgraphics(f,[figpath filesep sID '_BBcorr.pdf'],'ContentType','auto','Resolution',600)
% close(f)

f  = figure('color',[1 1 1]);
cb = colorbar;
cb.Ticks = [];
colormap(viridis)
axis off
% exportgraphics(f,[figpath filesep 'colorbar.pdf'],'ContentType','auto','Resolution',600)
% close(f)

%% fit model

group_fits = nan(numel(corrsource),3);
for i = 1:numel(corrsource)
    if i == sID
        [fits, d, signal] = fitbrainSA(corrsource{i});
        exempfits = fits;
    else
        [fits, ~, ~] = fitbrainSA(corrsource{i});
    end
    group_fits(i,:) = fits;
end

%% visualize FWHM

% plot
f = figure('color',[1 1 1]); hold on;
cmap = get(gca,'ColorOrder');
CI = bootci(1000,@mean,group_fits);
swarmchart(ones(length(group_fits)-1,1),group_fits(1:nsubs~=sID,1),144,'filled','^','MarkerFaceAlpha',0.5,'XJitterWidth',0.25,'MarkerFaceColor',cmap(6,:))
swarmchart(2.*ones(length(group_fits)-1,1),group_fits(1:nsubs~=sID,2),144,'filled','^','MarkerFaceAlpha',0.5,'XJitterWidth',0.25,'MarkerFaceColor',cmap(3,:))
swarmchart(3.*ones(length(group_fits)-1,1),group_fits(1:nsubs~=sID,3),144,'filled','^','MarkerFaceAlpha',0.5,'XJitterWidth',0.25,'MarkerFaceColor',cmap(4,:))
errorbar(1,mean(group_fits(:,1)),mean(group_fits(:,1))-CI(1,1),CI(2,1)-mean(group_fits(:,1)),'color','k')
errorbar(2,mean(group_fits(:,2)),mean(group_fits(:,2))-CI(1,2),CI(2,2)-mean(group_fits(:,2)),'color','k')
errorbar(3,mean(group_fits(:,3)),mean(group_fits(:,3))-CI(1,3),CI(2,3)-mean(group_fits(:,3)),'color','k')
xticks([1 2 3])
xticklabels({'lfp','low freq','high freq band'})
ylabel('FWHM (mm)')
title('Spatial Autocorrelation')
f.Position = [767 360 400 420];
swarmchart(1,group_fits(sID,1),144,'filled','^','MarkerFaceAlpha',0.5,'XJitterWidth',0.25,'MarkerFaceColor','k')
swarmchart(2,group_fits(sID,2),144,'filled','^','MarkerFaceAlpha',0.5,'XJitterWidth',0.25,'MarkerFaceColor','k')
swarmchart(3,group_fits(sID,3),144,'filled','^','MarkerFaceAlpha',0.5,'XJitterWidth',0.25,'MarkerFaceColor','k')
% exportgraphics(f,[figpath filesep '_groupavg.pdf'],'ContentType','auto','Resolution',600)
% close(f)

%% plot curve for exemplar sub

% gaussian 
model   = @(params, x) exp(-x.^2./(2*params(1).^2));
fittedsigma = exempfits./(2*sqrt(2*log(2)));
ticks   = -200:0.5:200;

% plot
mrkrsz = 15;
mrkralpha = 0.5;
f = figure('color',[1 1 1]);
cmap = get(gca,'ColorOrder');
scatter(d,signal(:,1),mrkrsz,'filled','MarkerFaceAlpha',mrkralpha,'MarkerFaceColor',cmap(6,:)); hold on;
scatter(d,signal(:,2),mrkrsz,'filled','MarkerFaceAlpha',mrkralpha,'MarkerFaceColor',cmap(3,:)); hold on;
scatter(d,signal(:,3),mrkrsz,'filled','MarkerFaceAlpha',mrkralpha,'MarkerFaceColor',cmap(4,:))
plot(ticks, model(fittedsigma(1),ticks), 'LineWidth', 2, 'color', cmap(6,:));
plot(ticks, model(fittedsigma(2),ticks), 'LineWidth', 2, 'color', cmap(3,:));
plot(ticks, model(fittedsigma(3),ticks), 'LineWidth', 2, 'color', cmap(4,:));
xlim([0 100])
xlabel('distance')
ylabel('corr')
legend({'LFP','low freq','high freq band'})
% exportgraphics(f,[figpath filesep sID '_curve.pdf'],'ContentType','auto','Resolution',600)
% close(f)
end

function [fits, d, signal] = fitbrainSA(source)

% fitbrainSA fits spatial autocorrelation decay curves for one subject
%
% input
%   source : struct with fields
%       elec.nativechanpos  [N x 3] channel positions (assumed mm)
%       LFPcorr             [N x N] correlation matrix (LFP)
%       LFcorr              [N x N] correlation matrix (LFB)
%       HFBcorr             [N x N] correlation matrix (HFB)
%
% output
%   fits     [1 x 3] FWHM (mm) for each signal (LFP, LFB, HFB)
%   d        [nPairs x 1] Euclidean distances (mm) for included pairs
%   signal   [nPairs x 3] correlations for included pairs

% prune by distance
lim         = 7.5;
keep        = source.elec.dist2surf <= lim | source.elec.ecog;
pos         = source.elec.nativechanpos(keep,:);
nch         = size(pos,1);

% unique intra-hemispheric channel pairs 
hemi        = double(pos(:,1) < 0);
sameHemi    = (hemi + hemi') ~= 1;  % true if NOT cross-hemisphere
upperNoDiag = triu(true(nch), 1);  % unique pairs, exclude diagonal
pairMask    = sameHemi & upperNoDiag;
idx         = find(pairMask);

% config data
D           = pdist2(pos, pos);                
d           = D(idx);
LFP         = source.LFPcorr(keep, keep);
LF          = source.LFcorr(keep, keep);
HFB         = source.HFBcorr(keep, keep);
signal      = [LFP(idx) LF(idx) HFB(idx)];

% fit model
model       = @(sigma, x) exp(-x.^2 ./ (2*sigma.^2));
initial     = 15;                         
fits        = nan(1, size(signal,2));
for k = 1:size(signal,2)
    
    y = signal(:,k);
    good = isfinite(d) & isfinite(y);
    dk = d(good);
    yk = y(good);

    if isempty(dk)
        fits(k) = NaN;
        continue;
    end

    sigfit = lsqcurvefit(model, initial, dk, yk);  % requires Optimization Toolbox
    fits(k) = 2*sqrt(2*log(2)) * sigfit;           % FWHM
end
end
