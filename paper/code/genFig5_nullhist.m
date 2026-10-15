function genFig5_nullhist(style, outroot, R)
% ------------------------------------------------------------------------------------------------
% Figure 5 max-cluster-mass null panels (SMOOTH vs RSF), with the observed clusters and their
% p-values marked. genFig5 draws them through this function during its run; called on its own, it
% redraws them in seconds from the numbers genFig5 saved, so axis or style changes need no re-run.
%
%   genFig5_nullhist                    redraw all three datasets from paper/data/genFig5_results.mat
%                                       into paper/figs/5/<dataset>/null_distribution.pdf
%   genFig5_nullhist('density')         the same nulls as kernel density estimates (null_density.pdf)
%   genFig5_nullhist(style, outroot)    write to <outroot>/<dataset>/ instead
%   genFig5_nullhist(style, outdir, R)  draw one dataset from its results struct R (fields null_SMOOTH,
%                                       null_RSF, cluster_mass, cluster_p) into outdir; genFig5's call
%
% Axes. x runs from 0 to a round limit at least 5% above the largest null value or observed cluster
% mass, with seven evenly spaced ticks, the last at the limit. y runs from 0 to a round limit at least
% 5% above the tallest bar (or density peak), so nothing touches the top, with six evenly spaced ticks
% labelled only at 0 and the limit. The histogram has 48 bins, 8 per x-tick interval, so the ticks
% fall on bin edges. The density is a Gaussian kernel estimate of each null (ksdensity's default
% bandwidth), reflected at 0 because a cluster mass cannot be negative.
%
% A redraw never overwrites: a panel it replaces is moved to figs/legacy/5_<dataset>_<timestamp>/.
% ------------------------------------------------------------------------------------------------
if nargin < 1 || isempty(style), style = 'histogram'; end
style = validatestring(style, {'histogram','density'});
if nargin >= 3
    draw(R, outroot, style); return
end
addpath(fullfile(fileparts(mfilename('fullpath')),'..','..'));
paths = smooth_setup('quiet');
if nargin < 2 || isempty(outroot), outroot = fullfile(paths.figs,'5'); end
RES = load(fullfile(paths.data,'genFig5_results.mat'));
stamp = datestr(now,'yyyymmdd_HHMMSS');
for name = {'primary','berezutskaya','cogitate'}
    if ~isfield(RES, name{1}), continue, end
    outdir = fullfile(outroot, name{1}); if ~isfolder(outdir), mkdir(outdir); end
    old = fullfile(outdir, panelname(style));
    if isfile(old)
        dst = fullfile(paths.figs,'legacy',sprintf('5_%s_%s', name{1}, stamp));
        if ~isfolder(dst), mkdir(dst); end
        movefile(old, dst); fprintf('  archived previous panel -> %s\n', dst);
    end
    draw(RES.(name{1}), outdir, style);
    fprintf('[%s] %s -> %s\n', name{1}, style, fullfile(outdir, panelname(style)));
end
end

function draw(R, outdir, style)
sm = double(R.null_SMOOTH(:)); rf = double(R.null_RSF(:));
mass = double(R.cluster_mass(:))'; prob = double(R.cluster_p(:))';
xmax = nicelim(1.05*max([sm; rf; mass(:); 1]), 6);
f = figure('color','w','position',[100 100 620 460]); hold on; box off
gr = [0.18 0.55 0.34]; or = [0.85 0.37 0.05];
switch style
    case 'histogram'
        edges = linspace(0, xmax, 49);
        histogram(sm, edges, 'FaceColor',gr, 'EdgeColor','none', 'FaceAlpha',0.65);
        histogram(rf, edges, 'FaceColor',or, 'EdgeColor','none', 'FaceAlpha',0.50);
        top = max([histcounts(sm, edges), histcounts(rf, edges)]); ylab = 'count';
    case 'density'
        xi = linspace(0, xmax, 600);
        ds = ksdensity(sm, xi, 'Support','positive', 'BoundaryCorrection','reflection');
        dr = ksdensity(rf, xi, 'Support','positive', 'BoundaryCorrection','reflection');
        patch([xi fliplr(xi)], [ds zeros(size(xi))], gr, 'FaceAlpha',0.45, 'EdgeColor',gr, 'LineWidth',1.2);
        patch([xi fliplr(xi)], [dr zeros(size(xi))], or, 'FaceAlpha',0.35, 'EdgeColor',or, 'LineWidth',1.2);
        top = max([ds dr]); ylab = 'probability density';
end
ymax = nicelim(1.05*top, 5);
xline(prctile(rf,95),'--','RSF 95th','Color',[0.7 0.3 0],'LabelVerticalAlignment','top');
xline(prctile(sm,95),'--','SMOOTH 95th','Color',[0.1 0.4 0.2],'LabelVerticalAlignment','bottom');
for k = 1:min(numel(mass),3)
    xline(mass(k), '-', sprintf('p=%.3f',prob(k)), 'Color','k', 'LineWidth',1.2, 'LabelOrientation','horizontal');
end
xlim([0 xmax]); ylim([0 ymax]);
set(gca,'FontSize',12,'XTick',linspace(0,xmax,7),'YTick',linspace(0,ymax,6), ...
    'YTickLabel',[{'0'}, repmat({''},1,4), {sprintf('%g',ymax)}]);
xlabel('max cluster mass (null)'); ylabel(ylab);
legend({'SMOOTH null','RSF null'}, 'Box','off','Location','northeast');
exportgraphics(f, fullfile(outdir, panelname(style)), 'ContentType','vector'); close(f);
end

function n = panelname(style)
if strcmp(style,'histogram'), n = 'null_distribution.pdf'; else, n = 'null_density.pdf'; end
end

function lim = nicelim(x, n)
% smallest limit above x that divides into n round steps (1, 1.2, 1.5, 2, 2.5, 3, 4, 5, 6 or 8 times
% a power of ten); the list runs past 10 in case log10 rounds an exact power of ten down a decade
p = 10^floor(log10(max(x/n, realmin)));
s = [1 1.2 1.5 2 2.5 3 4 5 6 8 10 12 15 20] * p;
lim = n * s(find(n*s > x, 1));
end
