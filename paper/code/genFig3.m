function genFig3
% -------------------------------------------------------------------------
% Figure 3. SMOOTH controls false positives under spatial autocorrelation.
%
% Consumes the 1/f-GRF cluster FWER results: one `results` .mat per beta in
% paper/data/fwer_sims_1f/, produced by a separate simulation pipeline (not distributed); this
% script reads its cached output.
%   - exemplars selected by beta VALUE = [1 29 50]
%   - linear beta axis, xlim [0 50] (log-spaced grid 1..50, incl. the empirical beta=29)
%   - FWER-vs-beta uses the stored per-beta counts (all betas) with a binomial CI, so only
%     the 3 exemplars need the full per-simulation arrays (distributions + QQ panels)
%
% Panels:
%   a) 1/f power-law spectra          b) example null cortical maps
%   c) null max-cluster-mass distributions + p-value QQ (SMOOTH vs naive)
%   d) FWER vs beta (SMOOTH vs naive) + beta->SA fit
% -------------------------------------------------------------------------

clear; clc; close all;

addpath(fullfile(fileparts(mfilename('fullpath')),'..','..'));
paths = smooth_setup('quiet');
paths.simData = fullfile(paths.data, 'fwer_sims_1f');
paths.figOut  = fullfile(paths.figs, '3');
if ~isfolder(paths.figOut), mkdir(paths.figOut); end
savefigs  = true;
EXEMPLARS = [1 29 50];     % rough / empirical / smooth
XLIM      = [0 50];
DISTXLIM  = [0 1000];      % cluster-mass distribution panels x-range
DISTYLIM  = [0 400];       % ... y-range (count)
DISTBW    = 25;            % ... bin width
rng(0);

% meshes
lh    = ft_read_headshape(fullfile(paths.freesurfer,'subjects','fsaverage5','surf','lh.pial'));
sph_l = ft_read_headshape(fullfile(paths.freesurfer,'subjects','fsaverage5','surf','lh.sphere.reg'));

%% ---- read per-beta results (paper/data/fwer_sims_1f/*beta_<b>_*.mat) ----
d = dir(fullfile(paths.simData, '*.mat'));
assert(~isempty(d), 'No .mat files in %s', paths.simData);
betaFromFile = nan(numel(d),1);
for i = 1:numel(d)
    tok = regexp(d(i).name, 'beta_(\d+)_', 'tokens', 'once');
    if ~isempty(tok), betaFromFile(i) = str2double(tok{1}); end
end
keep = isfinite(betaFromFile); d = d(keep); betaFromFile = betaFromFile(keep);
[betaFromFile, ord] = sort(betaFromFile); d = d(ord);
N = numel(d);
assert(N > 0, 'No beta-tagged results files found.');

beta = nan(N,1); fwerSmooth = nan(N,1); fwerNaive = nan(N,1);
nSim = nan(N,1); kSmooth = nan(N,1); kNaive = nan(N,1); meanSA = nan(N,1);
persim = containers.Map('KeyType','double','ValueType','any');
empSA  = [];
for b = 1:N
    S = load(fullfile(d(b).folder, d(b).name)); r = S.results;
    beta(b)       = r.simcfg.beta;
    fwerSmooth(b) = r.fwerPosSmooth;
    fwerNaive(b)  = r.fwerPosNaive;
    nSim(b)       = r.n;
    kSmooth(b)    = r.kSmooth;
    kNaive(b)     = r.kNaive;
    meanSA(b)     = mean(r.observedSA(:), 'omitnan');
    if isempty(empSA) && isfield(r,'empiricalSA'), empSA = r.empiricalSA(:); end
    if isfield(r,'smoothPosPvals') && numel(r.smoothPosPvals) > 1     % exemplar with per-sim arrays
        persim(r.simcfg.beta) = struct( ...
            'obsMaxMass',        r.observedMaxMass(:), ...
            'smoothPosPvals',    r.smoothPosPvals(:), ...
            'naivePosPvals',     r.naivePosPvals(:), ...
            'smoothNullMassP95', r.smoothNullMassP95(:), ...
            'naiveNullMassP95',  r.naiveNullMassP95(:));
    end
end

%% ---- plotting defaults ----
set(groot,'DefaultAxesLineWidth',1.25); set(groot,'DefaultLineLineWidth',1.5);
set(groot,'DefaultAxesXGrid','on');    set(groot,'DefaultAxesYGrid','on');
set(groot,'DefaultFigureColor','w');   set(groot,'DefaultAxesTickDir','out'); set(groot,'DefaultAxesBox','off');
set(groot,'DefaultAxesFontSize',32);   set(groot,'DefaultTextFontSize',32);  set(groot,'DefaultLegendFontSize',32);
blue = slanCM('Blues'); green = slanCM('Greens'); red = slanCM('Reds');
excmap = containers.Map(EXEMPLARS, {red, green, blue});     % beta value -> colormap

%% ---- 1/f eigenmodes, spectra, example maps ----
[Q, lam] = computeSphereLaplacian(sph_l, 50, 5000);
lam = max(lam(:), 0); c = median(lam);
bg = 1:50; SPEC = nan(numel(bg), size(Q,2)-1);
for k = 1:numel(bg)
    P = 1 ./ ((lam + c).^bg(k)); P(1) = 0;                 % drop DC (matches the sim)
    SPEC(k,:) = P(2:end);
end
Yex = nan(numel(EXEMPLARS), size(Q,1));
for e = 1:numel(EXEMPLARS)
    P = 1 ./ ((lam + c).^EXEMPLARS(e)); P(1) = 0;
    a = sqrt(P) .* randn(size(P)); y = Q * a; Yex(e,:) = y - mean(y);
end

% spectra
f = figure('Position',[0 0 400 900]);
plot(SPEC','Color',[0.7 0.7 0.7]); hold on; xlim([0 750]);
for e = 1:numel(EXEMPLARS), cm = excmap(EXEMPLARS(e)); plot(SPEC(EXEMPLARS(e),:),'Color',cm(192,:),'linewidth',2); end
set(gca,'Yscale','log'); ylabel('log (power)'); xlabel('mode'); xticks(0:250:750); grid off
if savefigs, exportgraphics(f,[paths.figOut filesep 'noisespectra.pdf'],"ContentType","vector"); end
close(f)

% example brainmaps
for e = 1:numel(EXEMPLARS)
    cm = excmap(EXEMPLARS(e));
    f = figure; ft_plot_mesh(lh,'vertexcolor','curv'); hold on;
    ft_plot_mesh(lh,'vertexcolor',Yex(e,:)','colormap',cm); view([-90 0]); colormap(cm);
    if savefigs, exportgraphics(f,[paths.figOut filesep sprintf('brainmap_beta%d.pdf',EXEMPLARS(e))],"ContentType","auto"); end
    close(f)
end

%% ---- distributions + QQ at the exemplars (need per-sim arrays) ----
for e = 1:numel(EXEMPLARS)
    be = EXEMPLARS(e);
    if ~isKey(persim, be), warning('genFig3:noPerSim','no per-sim arrays for beta=%d -- re-run that cell with the per-sim driver', be); continue; end
    ps = persim(be); cm = excmap(be);

    f = figure('position',[0 0 700 280]);
    histogram(ps.naiveNullMassP95, 'FaceColor',cm(128,:),'BinLimits',DISTXLIM,'binwidth',DISTBW,'EdgeAlpha',0,'FaceAlpha',0.75); hold on;
    histogram(ps.smoothNullMassP95,'FaceColor',cm(256,:),'BinLimits',DISTXLIM,'binwidth',DISTBW,'EdgeAlpha',0,'FaceAlpha',0.75);
    histogram(ps.obsMaxMass,       'FaceColor',[0.7 0.7 0.7],'BinLimits',DISTXLIM,'binwidth',DISTBW,'EdgeAlpha',0,'FaceAlpha',0.75);
    legend({'Naive Threshold','SMOOTH Threshold','Observed Max'}); xlim(DISTXLIM); ylim(DISTYLIM); yticks(0:100:400);
    if savefigs, exportgraphics(f,[paths.figOut filesep sprintf('distrib_beta%d.pdf',be)],"ContentType","vector"); end
    close(f)

    f = figure('position',[0 0 300 280]);
    plot(sort(ps.smoothPosPvals),'Color',cm(256,:)); hold on;
    plot(sort(ps.naivePosPvals), 'Color',cm(128,:));
    plot([0 numel(ps.smoothPosPvals)],[0 1],'Color',[0.7 0.7 0.7]);
    xticks(0:250:1000); yticks(0:0.25:1); xlabel('count (sorted)'); ylabel('p-val'); legend({'SMOOTH','Naive'},'Location','northwest');
    if savefigs, exportgraphics(f,[paths.figOut filesep sprintf('QQ_beta%d.pdf',be)],"ContentType","vector"); end
    close(f)
end

%% ---- beta -> SA fit (betaHat = crossing of sim SA and empirical SA) ----
p = polyfit(beta, meanSA, 3); xFit = linspace(min(beta), max(beta), 1000); yFit = polyval(p, xFit);
empSAval = median(empSA, 'omitnan'); [~, idx] = min(abs(yFit - empSAval)); betaHat = xFit(idx);
fprintf('Estimated beta = %.2f  (sim SA crosses empirical SA = %.3f)\n', betaHat, empSAval);
f = figure;
scatter(beta, meanSA, 60, 'filled'); hold on; plot(xFit, yFit, 'LineWidth', 2);
yline(empSAval,'--','LineWidth',1.5); xline(betaHat,'--','LineWidth',1.5);
xlabel('beta'); ylabel('SA'); xlim(XLIM);
legend({'Observed','Polynomial fit','Empirical SA','Estimated beta'},'Location','northwest');
title(sprintf('Estimated beta = %.2f', betaHat)); set(gca,'FontSize',14,'Box','off');
if savefigs, exportgraphics(f,[paths.figOut filesep 'betaSAfit.pdf']); end
close(f)

%% ---- FWER vs beta (from stored counts, all betas) ----
alpha = 0.05; pS = nan(N,1); pN = nan(N,1); ciS = nan(N,2); ciN = nan(N,2);
for b = 1:N
    [pS(b), ciS(b,:)] = binofit(kSmooth(b), nSim(b));
    [pN(b), ciN(b,:)] = binofit(kNaive(b),  nSim(b));
end
% --- manuscript panel: AXES (plot box) exactly 25.868 x 68.5 mm (w x h); figure larger for labels ---
C_SMOOTH = [0.46 0.20 0.60]; C_NAIVE = [0.95 0.58 0.10]; GRY = [0.30 0.30 0.30];
AXW = 2.5868; AXH = 6.85; ML = 1.30; MB = 1.05; MR = 0.35; MT = 0.35;   % axes size + margins (cm)
FW = ML+AXW+MR; FH = MB+AXH+MT;
f = figure('Units','centimeters','Position',[2 2 FW FH],'Color','w');
ax = axes('Parent',f,'Units','centimeters','Position',[ML MB AXW AXH]); hold(ax,'on'); set(ax,'FontSize',6);
hN = errorbar(ax, beta, pN, pN-ciN(:,1), ciN(:,2)-pN, '-o','Color',C_NAIVE, 'MarkerFaceColor','none','MarkerEdgeColor',C_NAIVE, 'MarkerSize',3,'CapSize',3,'LineWidth',0.6);
hS = errorbar(ax, beta, pS, pS-ciS(:,1), ciS(:,2)-pS, '-o','Color',C_SMOOTH,'MarkerFaceColor','none','MarkerEdgeColor',C_SMOOTH,'MarkerSize',3,'CapSize',3,'LineWidth',0.6);
yline(ax, alpha, '--', 'alpha', 'Color',GRY,'LineWidth',0.6,'FontSize',6,'LabelHorizontalAlignment','left','LabelVerticalAlignment','bottom');
xline(ax, 29,      '--', '\beta = 29',  'Color',GRY,'LineWidth',0.6,'FontSize',6,'LabelVerticalAlignment','middle','LabelHorizontalAlignment','left','Interpreter','tex');
xline(ax, betaHat, '--', '\beta = 9.1', 'Color',GRY,'LineWidth',0.6,'FontSize',6,'LabelVerticalAlignment','middle','LabelHorizontalAlignment','right','Interpreter','tex');
xlabel(ax,'beta'); ylabel(ax,'FWER');
set(ax,'XScale','log'); xlim(ax,[1 50]); ylim(ax,[0 1]); xticks(ax,[1 5 10 50]); xticklabels(ax,{'1','5','10','50'}); yticks(ax,0:0.25:1);
set(ax,'TickDir','in','Box','on','XGrid','on','YGrid','on','XMinorGrid','off','YMinorGrid','off','LineWidth',0.5);
lg = legend([hS hN],{'SMOOTH','Naive'},'Location','northwest','FontSize',5.5,'Box','off'); lg.ItemTokenSize = [8 6];
if savefigs
    set(f,'PaperUnits','centimeters','PaperPositionMode','manual','PaperSize',[FW FH],'PaperPosition',[0 0 FW FH]);
    print(f,'-dpdf','-vector',fullfile(paths.figOut,'FDR.pdf'));
end
close(f)
end

%% ========================= subfunction =========================
function [Q, lam, L, W] = computeSphereLaplacian(sph, graphsigma, k_modes)
pos = sph.pos; tri = sph.tri; n = size(pos,1);
if n == 0, Q = zeros(0,0); lam = []; L = sparse(0,0); W = sparse(0,0); return
elseif n == 1, Q = 1; lam = 0; L = sparse(1,1); W = sparse(1,1); return; end
if nargin < 3 || isempty(k_modes), k_modes = min(500, n-1); end
r = sqrt(sum(pos.^2, 2)); U = pos ./ r;
edges = [tri(:,[1 2]); tri(:,[2 3]); tri(:,[3 1])];
edges = unique(sort(edges, 2), 'rows');
i = edges(:,1); j = edges(:,2);
c = max(-1, min(1, sum(U(i,:) .* U(j,:), 2))); D = acos(c);
w = exp(-(D.^2) ./ (2 * graphsigma^2));
W = sparse([i; j], [j; i], [w; w], n, n);
dd = full(sum(W, 2)); dd(dd == 0) = 1;
Dm = spdiags(1 ./ sqrt(dd), 0, n, n); L = speye(n) - Dm * W * Dm; L = (L + L') / 2;
opts.isreal = true; opts.issym = true; opts.tol = 1e-8; opts.maxit = 1000;
try
    [Q, Lam] = eigs(L, k_modes, 'smallestabs', opts); lam = diag(Lam);
catch
    warning('eigs failed; falling back to full eig.'); [Qtmp, Lam] = eig(full(L));
    lam = diag(Lam); [lam, ord] = sort(lam, 'ascend'); Q = Qtmp(:, ord); return
end
[lam, ord] = sort(lam, 'ascend'); Q = Q(:, ord);
end
