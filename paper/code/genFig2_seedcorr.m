function genFig2_seedcorr(sub, seedlabel, climoverride)
% ------------------------------------------------------------------------------------------------
% Figure 2, panel a (Spatial Autocorrelation of Intracranial EEG Signals).
%
% Panel b (FWHM/scatter curves) is genFig2.m; panels c-f are genFig2_surrogates.m.
%
% Seed-correlation berry plot: the subject's own native pial surface with every electrode drawn as
% an ft_plot_cloud sphere coloured by its HFB correlation to a single seed contact in TPJ. This is
% one row of the HFB correlation matrix laid out anatomically -- the same decay that panel b
% plots against distance, shown on the cortex.
%
% Use as:
%   genFig2_seedcorr                          % default subject + seed + clim
%   genFig2_seedcorr('P19')                   % the subject, its default seed
%   genFig2_seedcorr('P19','TG21-TG22')       % a specific seed
%   genFig2_seedcorr('P19','TG21-TG22',[-0.25 0.75])
%
% SUBJECT.  P19 (pseudonym = its demo/source.mat index; corrsource 18), left hemisphere, 154 ECoG
% contacts. P19 was chosen from the cohort's ECoG grid participants for having the cleanest 8 x 8
% grid (its TG array: 112 bipolar channels spanning 85 x 82 mm of temporal-parietal cortex) and a
% pial mesh that is undistorted over that territory. Its frontal lobe renders smooth/under-gyrified,
% so crop the panel to the posterior two-thirds. Only P19's native surface is distributed.
% corrsource, demo/source.mat and the FreeSurfer recons use three different subject orderings; see
% export_natcortex.py for how the mapping was recovered and verified (recomputed dist2surf matches
% the stored values to 0.0000 mm).
%
% SEED.  Defaults to the contact nearest the centre of its subject's grid, so decay is visible in
% every direction: for P19 that is TG36-TG37, 7.0 mm from the TG centroid with 12 grid neighbours
% within 15 mm. The seed is drawn at its true value of r = 1, which makes it the unique maximum.
%
% CLIM.  The colour limits do real work here: off-seed r spans about [-0.05 0.6] with the bulk near
% zero, so [-1 1] parks almost every contact in the middle of the viridis ramp. The upper limit also
% decides whether the seed reads as a distinct maximum ([-0.25 1], nothing clipped) or merely as the
% brightest of several bright contacts ([-0.25 0.75], seed clipped but the decay better resolved).
% CLIM is stamped into every filename, so variants can be swept in one session.
%
% Data: paper/data/corrsource.mat (per-subject LFP/LF/HFB correlation matrices + native electrode
%   positions) and paper/data/natcortex_<SUB>.mat (native pial + curvature; the recon trees are not
%   part of the repository, so they are cached -- rebuild with export_natcortex.py).
% Run WITH a display. Panels -> paper/figs/2_seedcorr/*
% ------------------------------------------------------------------------------------------------

% subject | corrsource idx | hemi | grid | default seed | alternative interior seeds
SUBS = {
 'P19', 18, 'lh', 'TG:112',   'TG36-TG37',   {'TG36-TG44','TG28-TG36','TG29-TG30','TG37-TG45','TG44-TG45'}
};

if nargin < 1 || isempty(sub),         sub = 'P19';    end
if nargin < 2,                         seedlabel = ''; end
if nargin < 3,                         climoverride = []; end
si_sub = find(strcmp(SUBS(:,1), sub));
if isempty(si_sub), error('genFig2_seedcorr:noSub','unknown subject %s (have: %s)', sub, strjoin(SUBS(:,1)',', ')); end
CSUB = SUBS{si_sub,2};
if isempty(seedlabel), seedlabel = SUBS{si_sub,5}; end
if isempty(seedlabel)
    error('genFig2_seedcorr:noDefaultSeed', ...
        '%s (%s) has no default seed -- it was surveyed but not chosen; pass one explicitly', sub, SUBS{si_sub,4});
end

SIGNAL  = 'HFB';                    % 'HFB' | 'LF' | 'LFP'  (field <SIGNAL>corr in corrsource)
CLIM    = [0 1];                    % correlation colour limits: 0 to the seed's own r = 1. Off-seed r
                                    % only reaches ~0.5 and the most negative value on this subject is
                                    % -0.03, so nothing meaningful is clipped at the bottom.
RADIUS  = 2.5;                      % berry radius (mm)
if ~isempty(climoverride), CLIM = climoverride(:).'; end

set(0,'DefaultFigureVisible','on');
addpath(fullfile(fileparts(mfilename('fullpath')),'..','..'));
paths  = smooth_setup('quiet');
outdir = fullfile(paths.figs,'2_seedcorr'); if ~exist(outdir,'dir'), mkdir(outdir); end

% Illustrator-safe text (see genFigS6): interpreter 'none' + Arial.
gi = get(groot,{'defaultTextInterpreter'});
set(groot,'defaultTextInterpreter','none','defaultAxesFontName','Arial','defaultTextFontName','Arial');
restore = onCleanup(@() set(groot,'defaultTextInterpreter',gi{1})); %#ok<NASGU>

% ---------------------------------------------------------------- data
S = load(fullfile(paths.data,'corrsource.mat'));
s = S.corrsource{CSUB};
N = load(fullfile(paths.data,sprintf('natcortex_%s.mat',sub)));
cortex.pos  = double(N.natcortex.pos);
cortex.tri  = double(N.natcortex.tri);
cortex.curv = double(N.natcortex.curv(:));

elec = double(s.elec.nativechanpos);
lab  = string(s.label(:));
R    = double(s.([SIGNAL 'corr']));
si   = find(lab == string(seedlabel));
if isempty(si)
    error('genFig2_seedcorr:noSeed','seed %s not found in %s / corrsource{%d} (%d contacts)', ...
          seedlabel, sub, CSUB, numel(lab));
end
r = R(si,:).';

% Show only the hemisphere the seed sits on; the other one only ever occludes. Native surfaces are
% not split at x = 0, so use the vertex block boundary recorded by the exporter.
isL   = elec(si,1) < 0;
keepv = (1:size(cortex.pos,1))' > double(N.n_lh);
if isL, keepv = ~keepv; end
cortex = submesh(cortex, keepv);
keepe  = (elec(:,1) < 0) == isL;                % drop contacts on the hidden hemisphere
elec = elec(keepe,:); r = r(keepe); lab = lab(keepe); si = find(lab == string(seedlabel));

tag = sprintf('%s_%s_%s_clim%g_%g', sub, SIGNAL, seedlabel, CLIM(1), CLIM(2));
names = {sprintf('seedcorr_%s_lateral.png', tag), sprintf('cbar_%s.pdf', tag)};
archive_existing(outdir, paths.figs, names);
summarize(lab, r, elec, si, SIGNAL, seedlabel, sub, N, cortex, isL);

% ---------------------------------------------------------------- panels
vw = [90*(1-2*isL) 0];              % lateral: rh [90 0], lh [-90 0] (genFigS6)
berry(cortex, elec, r, CLIM, viridis, vw, RADIUS);
savepanel(outdir, sprintf('seedcorr_%s_lateral', tag));
cbar(outdir, sprintf('cbar_%s', tag), viridis, CLIM, sprintf('%s correlation (r)', SIGNAL));
fprintf('wrote %s to %s\n\n', tag, outdir);
end

% ==================================== helpers ====================================
function berry(c, elec, dv, cl, cmap, vw, rad)
% Berry convention: snap contacts to the nearest surface vertex and draw them as
% ft_plot_cloud 'surf' spheres. The seed is included like any other contact and so takes its true
% value, r = 1 -- the unique maximum as long as CLIM's upper limit is 1.
dv = dv(:); idx = knnsearch(c.pos, elec); e = c.pos(idx,:);
figure('color','w','position',[100 100 520 470]); axes('Position',[0.02 0.04 0.96 0.9]);
ft_plot_mesh(c,'vertexcolor','curv','edgecolor','none'); hold on;
ft_plot_cloud(e, dv, 'cloudtype','surf', 'radius',rad, 'scalerad','no', 'colormap',cmap, 'clim',cl);
view(vw); axis off vis3d;
camzoom(0.90);   % spheres sit proud of the surface and can extend past the axes' data limits, where
                 % they get clipped at the frame; zoom out so the whole cloud stays inside
end

function sub = submesh(c, keep)
vmap = zeros(size(c.pos,1),1); vmap(keep) = 1:nnz(keep);
fk = all(keep(c.tri),2);
sub.pos = c.pos(keep,:); sub.tri = vmap(c.tri(fk,:)); sub.curv = c.curv(keep);
end

function summarize(lab, r, elec, si, signal, seedlabel, sub, N, cortex, isL)
% Everything needed to justify the seed choice or pick a different one, printed once per run.
d = sqrt(sum((elec - elec(si,:)).^2, 2));
o = r; o(si) = NaN;
fprintf('\n%s seed correlation | %s (corrsource %d = source %d) | %s: %d contacts, %d vertices\n', ...
        signal, sub, N.corrsource_idx, N.source_idx, hemistr(isL), numel(lab), size(cortex.pos,1));
fprintf('  seed %s  r range [%.2f %.2f]  mean %.3f  n(r>0.3) = %d  n(r>0.5) = %d\n', ...
        seedlabel, min(o), max(o), mean(o,'omitnan'), sum(o>0.3), sum(o>0.5));
edges = [0 10 20 30 40 60 200];
fprintf('  mean r by distance from seed:\n');
for k = 1:numel(edges)-1
    m = d >= edges(k) & d < edges(k+1) & (1:numel(r)).' ~= si;
    if any(m), fprintf('    %3d-%3d mm  n=%3d  mean r = %+0.3f\n', edges(k), edges(k+1), nnz(m), mean(o(m),'omitnan')); end
end
[~,ord] = sort(o,'descend','MissingPlacement','last');
fprintf('  top 8 correlated contacts:\n');
for k = 1:min(8,numel(ord))
    j = ord(k); fprintf('    %-13s r = %+0.3f   d = %5.1f mm\n', lab(j), r(j), d(j));
end
end

function s = hemistr(isL), if isL, s = 'lh'; else, s = 'rh'; end, end

function archive_existing(outdir, figsroot, names)
% Never overwrite: move only the files this run is about to replace into figs/legacy/ (genFig5
% convention). Scoped to `names` rather than the whole folder so several subjects or clims can be
% swept in one session and left side by side for comparison.
hit = names(cellfun(@(n) isfile(fullfile(outdir,n)), names));
if isempty(hit), return; end
dst = fullfile(figsroot,'legacy',sprintf('2_seedcorr_%s',datestr(now,'yyyymmdd_HHMMSS'))); %#ok<TNOW1,DATST>
mkdir(dst);
for k = 1:numel(hit), movefile(fullfile(outdir,hit{k}), fullfile(dst,hit{k})); end
fprintf('archived %d previous file(s) to %s\n', numel(hit), dst);
end

function savepanel(outdir,name)
% Transparent-background raster PNG (house method; see genFigS6).
% exportgraphics(...,'BackgroundColor','none') silently writes opaque RGB in this R2023a/Rosetta
% build, so print the scene on white and on black with print (NOT exportgraphics, which auto-crops
% to non-background content and would give the two renders different sizes), recover true alpha from
% the difference (bg: W-K = 255(1-a); opaque: W = K -> a = 1), unpremultiply colour, write RGBA.
fw = fullfile(outdir,[name '_w.png']); fkp = fullfile(outdir,[name '_k.png']);
set(gcf,'InvertHardcopy','off','PaperPositionMode','auto'); set(findall(gcf,'type','axes'),'Color','none');
set(gcf,'Color','w'); print(gcf, fw, '-dpng','-r300');
set(gcf,'Color','k'); print(gcf, fkp,'-dpng','-r300'); close(gcf);
W = double(imread(fw)); K = double(imread(fkp));
alpha = min(max(1 - mean(W-K,3)/255, 0), 1);
C = min(max(K./max(alpha,1e-3), 0), 255);
[rr,cc] = find(alpha > 0.02);
if ~isempty(rr), rw = min(rr):max(rr); cw = min(cc):max(cc); C = C(rw,cw,:); alpha = alpha(rw,cw); end
imwrite(uint8(C), fullfile(outdir,[name '.png']), 'Alpha', alpha); delete(fw); delete(fkp);
end

function cbar(outdir,name,cmap,cl,label)
f = figure('color','w','position',[100 100 150 430]);
ax = axes('Parent',f,'Position',[0.05 0.06 0.02 0.88]); axis(ax,'off'); colormap(ax,cmap); set(ax,'CLim',cl);
cb = colorbar(ax,'Position',[0.32 0.06 0.22 0.88]); cb.Limits = cl; cb.Ticks = [cl(1) mean(cl) cl(2)];
set(cb,'FontSize',13); ylabel(cb,label,'FontSize',13);
exportgraphics(f,fullfile(outdir,[name '.pdf']),'ContentType','vector'); close(f);
end
