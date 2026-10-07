function genFigS6
% ------------------------------------------------------------------------------------------------
% Supplementary Figure 6. Comparison of analysis approaches across iEEG datasets. What standard
% FDR-corrected approaches recover, for all three datasets,
% with SMOOTH shown on the same axes for direct comparison.
% 4 x 3 grid rendered on fsaverage5 pial (lh + rh lateral per cell):
%   rows    = method: channel-level | ROI Desikan-Killiany | ROI Schaefer-300 | SMOOTH
%   columns = dataset: socialcog (primary) | Berezutskaya (ds003688) | COGITATE
% First three rows are ONE-TAILED POSITIVE, Benjamini-Hochberg FDR (q<0.05); the SMOOTH row is
% ONE-TAILED POSITIVE, cluster-based FWE (p<0.05) -- each method shown under its own native
% correction. Channel-level flags scattered electrodes; parcellation flags a few coarse regions or
% (finer Schaefer / weaker datasets) nothing, whereas SMOOTH recovers spatially-coherent clusters.
%
% Channel panels: significant electrodes as berries coloured by the dataset DV (grey = non-significant).
% ROI panels: fs5 vertices of FDR-significant parcels filled by the parcel group-t.
% SMOOTH panels: fs5 vertices inside FWE-significant clusters filled by the group t, on a COMPRESSED
%   [0 3] t axis with its own colorbar (the ROI rows keep [0 6]). Thresholded, NOT the raw t-map --
%   the other rows are thresholded too, so an unthresholded SMOOTH map would flatter it unfairly.
%   Cluster-forming threshold is per-dataset (0.05 / 0.01 / 0.05; see export_smooth_sigmaps).
%   NOTE: the SMOOTH row is on a different colour scale from the ROI rows.
% Data: paper/data/methodscomparison_pvals.mat (from export_methodscomparison_pvals.py) and
%   paper/data/smooth_sigmaps.mat (from export_smooth_sigmaps.m -- same cfg/seed as genFig5, so the
%   clusters here are the Figure 5 clusters). Run WITH a display.
% Panels -> paper/figs/S6/*.pdf
% ------------------------------------------------------------------------------------------------
set(0,'DefaultFigureVisible','on');
addpath(fullfile(fileparts(mfilename('fullpath')),'..','..'));
paths  = smooth_setup('quiet');
D      = load(fullfile(paths.data,'methodscomparison_pvals.mat'));
outdir = fullfile(paths.figs,'S6'); if ~exist(outdir,'dir'), mkdir(outdir); end

smoothfile = fullfile(paths.data,'smooth_sigmaps.mat');
if ~isfile(smoothfile)
    error('genFigS6:noSmooth', ...
        ['%s not found -- run export_smooth_sigmaps first (it runs SMOOTH on all three ' ...
         'datasets and caches the FWE-significant clusters on the fs5 mesh).'], smoothfile);
end
S = load(smoothfile);
fprintf('SMOOTH row from %s (NPERM=%d, seed %d)\n', smoothfile, S.NPERM, S.SEED);

% Illustrator-safe text: none + Arial, no tex/latex fonts.
gi = get(groot,{'defaultTextInterpreter','defaultLegendInterpreter'});
set(groot,'defaultTextInterpreter','none','defaultAxesFontName','Arial','defaultTextFontName','Arial');
restore = onCleanup(@() set(groot,'defaultTextInterpreter',gi{1})); %#ok<NASGU>

cortex.pos = double(D.fs5_pos); cortex.tri = double(D.fs5_tri); cortex.curv = double(D.fs5_curv(:));
try rdbu = flipud(slanCM('RdBu')); catch, rdbu = jet(256); end
% Red-blue family throughout. DV is signed -> full diverging map. ROI group-t is one-tailed positive
% (all >0), so it takes the UPPER half only (neutral -> red); a full diverging map on a [0 6] clim
% would paint weakly-positive parcels deep blue, which reads as a negative effect.
redsHalf = rdbu(round(size(rdbu,1)/2):end, :); TCL = [0 6];
SCL = [0 3];   % SMOOTH row only: compressed t axis. SMOOTH's significant t runs well above the ROI
               % range, so [0 6] left most of the cluster in the lower half of the ramp; [0 3] uses
               % the full ramp. NOTE this means the SMOOTH row is NOT on the ROI rows' colour scale
               % -- it needs its own colorbar (below) and a caption line saying so.
DS   = {'socialcog','Primary'; 'berezutskaya','Berezutskaya'; 'cogitate','COGITATE'};
VIEWS = {'lh',[-90 0]; 'rh',[90 0]};

for di = 1:size(DS,1)
    ds = DS{di,1}; dvlim = double(D.([ds '_dvlim']));
    cp  = double(D.([ds '_chanpos'])); dv = double(D.([ds '_dv'])); sig = logical(D.([ds '_chan_sig']));
    for v = 1:size(VIEWS,1)
        hemi = VIEWS{v,1}; vw = VIEWS{v,2};
        % no title() on any panel -- keeps the raster free of baked-in text so it needs no cropping in
        % Illustrator; method/dataset/hemi are all encoded in the filename.
        % --- channel-level (berries) ---
        berry(cortex, cp, dv, sig, [-dvlim dvlim], rdbu, vw);
        savepanel(outdir, sprintf('channel_%s_%s', ds, hemi));
        % --- DK ROI ---
        overlay(cortex, double(D.([ds '_dk_sigt'])), TCL, redsHalf, vw);
        savepanel(outdir, sprintf('dk_%s_%s', ds, hemi));
        % --- Schaefer-300 ROI ---
        overlay(cortex, double(D.([ds '_sch_sigt'])), TCL, redsHalf, vw);
        savepanel(outdir, sprintf('sch_%s_%s', ds, hemi));
        % --- SMOOTH (FWE-significant clusters), compressed [0 3] t axis (see SCL) ---
        overlay(cortex, double(S.([ds '_smooth_sigt'])), SCL, redsHalf, vw);
        savepanel(outdir, sprintf('smooth_%s_%s', ds, hemi));
    end
    cbar(outdir, sprintf('cbar_dv_%s', ds), rdbu, [-dvlim dvlim], sprintf('%s DV', ds));
end
cbar(outdir, 'cbar_roi_t',    redsHalf, TCL, 'ROI group t');       % DK + Schaefer rows
cbar(outdir, 'cbar_smooth_t', redsHalf, SCL, 'SMOOTH group t');    % SMOOTH row (compressed axis)
fprintf('wrote panels to %s\n', outdir);
end

% ==================================== helpers ====================================
function overlay(cortex, field, cl, cmap, vw)
field = field(:);
figure('color','w','position',[100 100 520 470]); axes('Position',[0.02 0.04 0.96 0.9]);
ft_plot_mesh(cortex, 'vertexcolor','curv', 'edgecolor','none'); hold on;
keep = isfinite(field);
if any(keep), sub = submesh(cortex, keep); ft_plot_mesh(sub, 'vertexcolor', field(keep), 'edgecolor','none'); end
colormap(cmap); clim(cl); view(vw); axis off;
end

function berry(cortex, elec, dv, sig, cl, cmap, vw)
% Berry convention: snap electrodes to the nearest fs5 vertex and draw them as
% ft_plot_cloud 'surf' spheres (radius 2.5, no radius scaling) -- NOT flat scatter3 dots.
% Two calls are safe despite the differing colormaps: for cloudtype 'surf', ft_plot_cloud bakes an
% explicit per-sphere FaceColor (ft_plot_cloud.m:671), so the grey spheres are not restyled when the
% second call resets the axes colormap.
dv = dv(:); sig = logical(sig(:));
idx = knnsearch(cortex.pos, elec); e = cortex.pos(idx,:);
figure('color','w','position',[100 100 520 470]); axes('Position',[0.02 0.04 0.96 0.9]);
ft_plot_mesh(cortex, 'vertexcolor','curv', 'edgecolor','none'); hold on;
if any(~sig)   % non-significant: smaller, flat grey spheres so they recede
    grey = repmat([0.62 0.62 0.62], 256, 1);
    ft_plot_cloud(e(~sig,:), zeros(nnz(~sig),1), 'cloudtype','surf', 'radius',1.6, ...
                  'scalerad','no', 'colormap',grey, 'clim',[-1 1]);
end
if any(sig)    % significant: full-size spheres coloured by the dataset DV
    ft_plot_cloud(e(sig,:), dv(sig), 'cloudtype','surf', 'radius',2.5, ...
                  'scalerad','no', 'colormap',cmap, 'clim',cl);
end
view(vw); axis off;
end

function sub = submesh(cortex, keep)
vmap = zeros(size(cortex.pos,1),1); vmap(keep) = 1:nnz(keep);
fk = all(keep(cortex.tri),2);
sub.pos = cortex.pos(keep,:); sub.tri = vmap(cortex.tri(fk,:));
end

function savepanel(outdir,name)
% Transparent-background raster PNG (same trick as genFigS5).
% exportgraphics(...,'BackgroundColor','none') silently writes opaque RGB under this R2023a/Rosetta
% build, so render the scene on white and on black with print (NOT exportgraphics, which auto-crops to
% non-background content and would give the two renders different sizes), recover true alpha from the
% difference (bg: W-K=255(1-a); opaque: W=K -> a=1), unpremultiply colour, write RGBA.
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
cb = colorbar(ax,'Position',[0.32 0.06 0.22 0.88]); cb.Limits = cl; set(cb,'FontSize',13); ylabel(cb,label,'FontSize',13);
exportgraphics(f,fullfile(outdir,[name '.pdf']),'ContentType','vector'); close(f);
end
