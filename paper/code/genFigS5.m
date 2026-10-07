function genFigS5
% ------------------------------------------------------------------------------------------------
% Supplementary Figure 5. Detection sensitivity with false-discovery-rate-corrected atlas analyses.
% COMPANION to genFig4.m: identical panels, but DK/Schaefer-300 significance is assessed by per-parcel
% one-tailed t-test + Benjamini-Hochberg FDR (q<=0.05) instead of the matched sign-flip max-statistic
% permutation FWER used in the main figure. SMOOTH is untouched (reused as-is).
%   a) left: sensitivity KDE (log) | right: sensitivity vs coverage | inset: specificity at sens>50%
%   b) brainmaps (SMOOTH / DK-FDR / Schaefer-300-FDR)
%   c) pairwise 2D histograms with "more sensitive" percentages (as in genFig4.m's hist2pair)
% Data: paper/data/benchmark_fig_data_fdr.mat (same field layout as benchmark_fig_data.mat; built by
% a one-off Python export from a separate simulation pipeline, not distributed: full 14,075-vertex
% FDR sweep, niter=1000, noise=1f beta=29 a=0.53).
% Run WITH a display. Panels -> paper/figs/S5/*.pdf
% ------------------------------------------------------------------------------------------------
clear; clc; close all;
addpath(fullfile(fileparts(mfilename('fullpath')),'..','..'));
paths  = smooth_setup('quiet');
D      = load(fullfile(paths.data,'benchmark_fig_data_fdr.mat'));
ds     = char(D.dataset);
outdir = fullfile(paths.figs,'S5'); if ~exist(outdir,'dir'), mkdir(outdir); end
mag = magma(256);
COL = {[0.776 0.157 0.157], [0.082 0.396 0.753], [0.180 0.490 0.196]};   % SMOOTH / DK / Schaefer
NAM = {'SMOOTH','ROI: DK (FDR)','ROI: Schaefer-300 (FDR)'};

% ================= BRAINMAPS (individual pial panels) =================
% ROI (vertex) level only -- '_atlas' fields are a bit-identical duplicate of '_roi' in the
% source .mat, so plotting both would just produce two identical panels per method.
methods = {'smooth','SMOOTH'; 'dk_roi','DK, FDR'; 'sch_roi','Schaefer-300, FDR'};
views = {'lh','lateral',[-90 0]; 'lh','medial',[90 0]; 'rh','lateral',[90 0]; 'rh','medial',[-90 0]};
for m = 1:size(methods,1)
    for v = 1:size(views,1)
        h = views{v,1}; vn = views{v,2}; vw = views{v,3};
        c.pos = double(D.([h '_pial_pos'])); c.tri = double(D.([h '_pial_tri'])); c.curv = double(D.([h '_curv'])(:));
        field = double(D.([h '_' methods{m,1}])(:));
        overlay(c, field, isfinite(field), [0 1], mag, vw);
        % no title() -- keeps the raster free of baked-in text so it doesn't need cropping in Illustrator;
        % method/hemi/view are already encoded in the filename.
        savepanel(outdir, sprintf('map_%s_%s_%s', methods{m,1}, h, vn));
    end
end
cbar(outdir, 'map_cbar', mag, [0 1], 'detection sensitivity');

% ================= 2D PLOTS =================
cov = double(D.seed_cov(:)); sm = double(D.seed_smooth(:)); dk = double(D.seed_dk(:)); sch = double(D.seed_sch(:));
smS = double(D.seed_smooth_spec(:)); dkS = double(D.seed_dk_spec(:)); schS = double(D.seed_sch_spec(:));
vir = viridiscmap(256);
curve3(outdir, cov, {sm,dk,sch}, NAM, COL, ds);
hist2pair(outdir,'hist2_smooth_dk',  dk,  sm, 'ROI: DK sensitivity (FDR)',          'SMOOTH sensitivity', 'SMOOTH vs DK (FDR)', vir);
hist2pair(outdir,'hist2_smooth_sch', sch, sm, 'ROI: Schaefer-300 sensitivity (FDR)','SMOOTH sensitivity', 'SMOOTH vs Schaefer-300 (FDR)', vir);
hist2pair(outdir,'hist2_dk_sch',     sch, dk, 'ROI: Schaefer-300 sensitivity (FDR)','ROI: DK sensitivity (FDR)','DK vs Schaefer-300 (FDR)', vir);
distkde(outdir,'distributions_kde_log',   'log',    {sm,dk,sch}, NAM, COL, ds);
% Specificity restricted to vertices where each method's
% own sensitivity exceeds 50%, i.e. "in the regime where we're actually likely to detect something,
% how well-controlled is the false-positive rate?"
specbar(outdir, 'specificity_at_sens50', {sm,dk,sch}, {smS,dkS,schS}, NAM, COL, 0.5);
fprintf('wrote FDR benchmark panels to %s\n', outdir);
end

% ==================================== helpers (identical to genFig4.m) ====================================
function overlay(cortex, field, showmask, cl, cmap, vw)
field = field(:);
figure('color','w','position',[100 100 560 500]); axes('Position',[0.02 0.06 0.96 0.9]);
ft_plot_mesh(cortex, 'vertexcolor','curv', 'edgecolor','none'); hold on;
keep = showmask(:) & isfinite(field);
if any(keep), sub = submesh(cortex, keep); ft_plot_mesh(sub, 'vertexcolor', field(keep), 'edgecolor','none'); end
colormap(cmap); clim(cl); view(vw);
end
function sub = submesh(cortex, keep)
vmap = zeros(size(cortex.pos,1),1); vmap(keep) = 1:nnz(keep);
fk = all(keep(cortex.tri),2);
sub.pos = cortex.pos(keep,:); sub.tri = vmap(cortex.tri(fk,:));
end
function savepanel(outdir,name)
% Transparent-background raster PNG. exportgraphics(...,'BackgroundColor','none') silently fails
% under this R2023a/Rosetta build (writes opaque RGB, no alpha) -- use print's two-background
% compositing instead: render the same figure on white then black, recover true alpha from the
% difference (bg: W-K=255(1-a); opaque brain: W=K -> a=1), unpremultiply colour, write RGBA.
% print (not exportgraphics) so both renders come out the same fixed size (exportgraphics auto-crops
% to non-background content, which crops white/black differently and breaks the differencing).
fw = fullfile(outdir,[name '_w.png']); fkp = fullfile(outdir,[name '_k.png']);
set(gcf,'InvertHardcopy','off','PaperPositionMode','auto'); set(findall(gcf,'type','axes'),'Color','none');
set(gcf,'Color','w'); print(gcf, fw, '-dpng','-r300');
set(gcf,'Color','k'); print(gcf, fkp, '-dpng','-r300'); close(gcf);
W = double(imread(fw)); K = double(imread(fkp));
alpha = min(max(1 - mean(W-K,3)/255, 0), 1);
C = min(max(K./max(alpha,1e-3), 0), 255);
[rr,cc] = find(alpha > 0.02);
if ~isempty(rr), r = min(rr):max(rr); c = min(cc):max(cc); C = C(r,c,:); alpha = alpha(r,c); end
imwrite(uint8(C), fullfile(outdir,[name '.png']), 'Alpha', alpha); delete(fw); delete(fkp);
end
function cbar(outdir,name,cmap,cl,label)
f = figure('color','w','position',[100 100 150 430]);
ax = axes('Parent',f,'Position',[0.05 0.06 0.02 0.88]); axis(ax,'off'); colormap(ax,cmap); set(ax,'CLim',cl);
cb = colorbar(ax,'Position',[0.32 0.06 0.22 0.88]); cb.Limits = cl; set(cb,'FontSize',14); ylabel(cb,label,'FontSize',14);
exportgraphics(f,fullfile(outdir,[name '.pdf']),'ContentType','vector'); close(f);
end
function curve3(outdir,cov,ys,names,cols,ds)
edges = 3:22; ctr = 0.5*(edges(1:end-1)+edges(2:end));
f = figure('color','w','position',[100 100 560 450]); hold on; box off; grid on
for i = 1:numel(ys)
    y = ys{i}; m = nan(1,numel(ctr)); e = m;
    for b = 1:numel(ctr)
        k = cov>=edges(b) & cov<edges(b+1);
        if nnz(k) > 20, m(b) = mean(y(k)); e(b) = 1.96*std(y(k))/sqrt(nnz(k)); end
    end
    errorbar(ctr,m,e,'o-','LineWidth',2,'Color',cols{i},'MarkerFaceColor',cols{i},'MarkerEdgeColor','none','CapSize',3,'DisplayName',names{i});
end
xlabel('coverage (# subjects)'); ylabel('detection sensitivity'); ylim([0 1]); set(gca,'FontSize',12);
legend('Location','northwest','Box','off'); title(sprintf('Detection sensitivity vs coverage, FDR (%s)',ds));
exportgraphics(f,fullfile(outdir,'sens_vs_coverage.pdf'),'ContentType','vector'); close(f);
end
function t = firsttok(s), s = strrep(s,'ROI: ',''); c = strsplit(s,' '); t = c{1}; end
function hist2pair(outdir,name,x,y,xl,yl,ttl,cmap,statxl,statyl)
% statxl/statyl: short labels for the "(Y>X: NN%)" title stat. Default to firsttok(xl/yl) (the method
% name, e.g. 'DK' vs 'Schaefer-300') -- right for pairwise same-metric comparisons. Pass explicitly
% for same-method different-metric panels (e.g. sensitivity vs specificity), where the method-name
% token is identical on both axes and would otherwise print a confusing "X>X" stat.
if nargin < 9, statxl = firsttok(xl); end
if nargin < 10, statyl = firsttok(yl); end
edges = linspace(0,1,41);
f = figure('color','w','position',[100 100 560 470]);
histogram2(x,y,'XBinEdges',edges,'YBinEdges',edges,'DisplayStyle','tile','ShowEmptyBins','off','EdgeColor','none'); hold on
colormap(cmap); set(gca,'ColorScale','log'); clim([1 1e4]);   % shared 10^0..10^4 across all pairs
cb = colorbar; ylabel(cb,'# seeds'); cb.Ticks = [1 10 100 1000 10000];
plot([0 1],[0 1],'--','Color',[0.15 0.15 0.15],'LineWidth',1);
xlim([0 1]); ylim([0 1]); axis square; xlabel(xl); ylabel(yl); set(gca,'FontSize',12);
% Plain comparison, matching genFig4.m and the manuscript text: "more sensitive" = mean(y>x).
title(sprintf('%s   (%s>%s: %.0f%%)', ttl, statyl, statxl, mean(y>x)*100));
exportgraphics(f,fullfile(outdir,[name '.pdf']),'ContentType','vector'); close(f);
end
function specbar(outdir,name,sens,spec,names,cols,thresh)
f = figure('color','w','position',[100 100 380 420]); hold on; box off
nb = numel(sens); m = nan(1,nb); e = m; ns = zeros(1,nb);
for i = 1:nb
    k = sens{i} > thresh; ns(i) = nnz(k);
    v = spec{i}(k); m(i) = mean(v); e(i) = 1.96*std(v)/sqrt(numel(v));
end
for i = 1:nb
    bar(i, m(i), 0.6, 'FaceColor', cols{i}, 'EdgeColor','none');
end
errorbar(1:nb, m, e, 'k.', 'LineWidth',1.5,'CapSize',5);
for i = 1:nb, text(i, m(i)+0.006, sprintf('n=%d',ns(i)), 'HorizontalAlignment','center','FontSize',9,'Color',[.3 .3 .3]); end
set(gca,'XTick',1:nb,'XTickLabel',names,'XTickLabelRotation',15,'FontSize',11);
ylim([0.9 1.0]); ylabel('specificity (1 - FPR)');
title(sprintf('Specificity at sensitivity > %d%%', round(thresh*100)));
exportgraphics(f,fullfile(outdir,[name '.pdf']),'ContentType','vector'); close(f);
end
function cm = viridiscmap(n)
a=[0.267004 0.004874 0.329415; 0.282623 0.140926 0.457517; 0.253935 0.265254 0.529983;
   0.206756 0.371758 0.553117; 0.163625 0.471133 0.558148; 0.127568 0.566949 0.550556;
   0.134692 0.658636 0.517649; 0.266941 0.748751 0.440573; 0.477504 0.821444 0.318195;
   0.741388 0.873449 0.149561; 0.993248 0.906157 0.143936];
x=linspace(0,1,size(a,1)); xi=linspace(0,1,n)';
cm=[interp1(x,a(:,1),xi) interp1(x,a(:,2),xi) interp1(x,a(:,3),xi)];
end
function distkde(outdir,name,scale,ys,names,cols,ds)
grid = linspace(0,1,400); islog = strcmp(scale,'log'); floorv = 1e-3;
lo = 0; if islog, lo = floorv; end
f = figure('color','w','position',[100 100 560 470]); hold on; box off
for i = 1:numel(ys)
    v = min(max(ys{i},1e-4),1-1e-4);
    d = ksdensity(v,grid,'Support',[0 1],'BoundaryCorrection','reflection','Bandwidth',0.03);
    up = d; if islog, up = max(d,floorv); end
    fill([grid fliplr(grid)],[up lo*ones(1,numel(grid))],cols{i},'FaceAlpha',0.18,'EdgeColor','none','HandleVisibility','off');
    plot(grid,d,'LineWidth',2,'Color',cols{i},'DisplayName',sprintf('%s (median %.2f)',names{i},median(ys{i})));
end
xlim([0 1]); set(gca,'YScale',scale); if islog, ylim([floorv inf]); end
xlabel('detection sensitivity'); ylabel('density'); set(gca,'FontSize',12);
legend('Location','north','Box','off'); title(sprintf('Sensitivity densities, FDR (KDE, %s y) - %s',scale,ds));
exportgraphics(f,fullfile(outdir,[name '.pdf']),'ContentType','vector'); close(f);
end
