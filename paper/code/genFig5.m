function genFig5(NPERM, only, savedir)
% ------------------------------------------------------------------------------------------------
% Figure 5. SMOOTH recovers real, domain-appropriate effects in THREE independent iEEG datasets.
% One column per dataset (primary social-cognition | Berezutskaya ds003688 | COGITATE), each with:
%   (top)    coverage map (# subjects per vertex, 0-20)
%   (mid)    per-electrode dependent variable (berry), dataset-specific scale
%   (bottom) group-level SMOOTH t-map (+-4)
%   + surrogate-null validation: per-subject spatial-autocorrelation (Moran's I) scatter (empirical vs
%     surrogate mean +- SD, from the SMOOTHdiag values in paper/data/fig5_sa_moran.mat, written by
%     export_fig5_sa.m), PC
%     cumulative-variance of the surrogate ensemble (SMOOTH vs RSF, nPC90), and the max-cluster-mass
%     null distribution (SMOOTH vs RSF) with observed clusters + p-values annotated.
%
% This script ALSO produces Supplementary Figure 7 (RSF reconstruction) as a side effect, rather than
% as a separate entry point: both depend on the exact same per-dataset SMOOTH/RSF run (same NPERM,
% same seed), computed once here (pc_reconstruction, called from the per-dataset loop below). Splitting
% it into a standalone script would mean either duplicating that run (real cluster time, and a second
% place for its seed/parameters to drift out of sync with this one) or reading back the group-level
% surrogate ensembles, which the repo does not ship (multi-GB).
%
% All three datasets are reduced to the same `source` format (per-subject .stat + .elec.chanpos in
% fsaverage coords); SMOOTH is re-run live from those sources (the repo ships the small sources, not the
% multi-GB surrogate ensembles). Run WITH a display. Optional arg NPERM (default 1000; use e.g. 100 for a
% quick preview). Optional arg SAVEDIR: also save each dataset's full SMOOTH and RSF output structs,
% surrogate ensembles included (several GB per dataset), so panels can be redrawn without re-running.
%
% OUTPUT
%   paper/figs/5/<dataset>/   Figure 5 panels
%   paper/figs/S7/<dataset>/  Supplementary Figure 7 panels (RSF-vs-SMOOTH t-map reconstruction)
%   paper/figs/S7/tencomponents.pdf  Supplementary Figure 7 ten-component comparison, all three
%                             datasets (genFigS7_tencomponents, from genFig5_results.mat)
% The null-distribution panels are drawn by genFig5_nullhist, which can also redraw them alone from
% genFig5_results.mat (genFig5_nullhist; genFig5_nullhist('density') for kernel density curves).
%   paper/figs/legacy/        any previous contents of the above, timestamped
% Nothing is overwritten: archive_existing() moves prior output to legacy/ before each run.
% Cluster-forming alpha is PER DATASET (see the DS table) and must match export_smooth_sigmaps.m.
% ------------------------------------------------------------------------------------------------
clc; close all;
if nargin < 1 || isempty(NPERM), NPERM = 1000; end
if nargin < 2, only = {}; end                                % e.g. genFig5(1000,'primary') to redo one
if nargin < 3, savedir = ''; end
if ischar(only) || isstring(only), only = cellstr(only); end
SEED = 42;
addpath(fullfile(fileparts(mfilename('fullpath')),'..','..'));
paths = smooth_setup('quiet');

cortex = load_fsaverage_mesh(paths.freesurfer);
SA     = load(fullfile(paths.data,'fig5_sa_moran.mat'));     % per-participant Moran's I (SMOOTHdiag) for panel F
if ~isfield(SA, 'generated_by')
    error('genFig5:staleSA', 'paper/data/fig5_sa_moran.mat predates export_fig5_sa.m; regenerate it with export_fig5_sa.');
end
grn = slanCM('Greens'); try rdbu = flipud(slanCM('RdBu')); catch, rdbu = jet(256); end

% dataset | source file | DV label | DV colour limit (+-) | cluster-forming alpha
%
% CLUSTER-FORMING THRESHOLD IS PER-DATASET, stated explicitly rather than left to the SMOOTHstat
% default -- it is not uniform across datasets and a silent default misreproduces this figure. These
% values MUST stay in sync with export_smooth_sigmaps.m, which feeds the SMOOTH row of the
% methods-comparison supplement; if the two disagree, the same clusters are drawn two different ways.
% Berezutskaya uses 0.01 (z = 2.33), the threshold the manuscript reports for this dataset. Its three
% clusters (lh superior/middle temporal, lh precentral/pars opercularis, rh insula) are already separate
% at 0.05, where they are larger (13,110 / 7,220 / 1,908 vertices).
DS = { 'primary',      fullfile(paths.demo,'source.mat'),              'z-value',  10, 0.05 ;
       'berezutskaya', fullfile(paths.data,'source_berezutskaya.mat'), 'Fisher-z',  2, 0.01 ;
       'cogitate',     fullfile(paths.data,'source_cogitate.mat'),     'R',         1, 0.05 };

% Every derived number this figure reports is harvested into RES and written to
% paper/data/genFig5_results.mat at the end, so captions/text can be sourced from a file instead of
% re-running this job or reading numbers off a PDF. Provenance (NPERM/SEED/clusteralpha) travels with
% it -- nPC90 in particular is ensemble-size dependent and is meaningless without NPERM.
RES = struct('NPERM',NPERM,'SEED',SEED,'generated',datestr(now,'yyyy-mm-dd HH:MM:SS'));

for di = 1:size(DS,1)
    name = DS{di,1}; srcfile = DS{di,2}; dvlab = DS{di,3}; dvlim = DS{di,4}; ca = DS{di,5};
    if ~isempty(only) && ~ismember(name, only)
        fprintf('[%s] skipped (dataset filter)\n', name); continue
    end
    outdir = fullfile(paths.figs,'5',name); archive_existing(outdir);
    L = load(srcfile); source = L.source;
    fprintf('\n[%s] %d subjects, %d electrodes | NPERM=%d | clusteralpha=%.3f (z=%.2f)\n', ...
            name, numel(source), sum(cellfun(@(s) numel(s.stat), source)), NPERM, ca, norminv(1-ca));

    % ---- run SMOOTH + RSF (surrogates kept at the group level for the PC-variance panel) ----
    cfg = []; cfg.fshome = paths.freesurfer; cfg.numrandomization = NPERM;
    cfg.kernelwidth = 12; cfg.graphsigma = 8; cfg.smooth = 'sphere'; cfg.kernel = 'gaussian';
    cfg.rankrescale = 'exact'; cfg.tail = 'positive'; cfg.randomseed = SEED; cfg.clusteralpha = ca;
    cfg.keepmaps = 'no'; cfg.keepsubsurrogates = 'no'; cfg.keepgroupsurrogates = 'yes'; cfg.normalize = 'no';
    stat    = SMOOTHstat(cfg, source{:});
    statRSF = SMOOTHdummy(cfg, source{:});
    fprintf('  cluster p-values: %s\n', mat2str(round([stat.posclusters.prob],4)));
    R = struct();
    R.clusteralpha  = ca;
    R.nsub          = numel(source);
    R.nchan         = sum(cellfun(@(s) numel(s.stat), source));
    R.cluster_p     = [stat.posclusters.prob];
    R.cluster_mass  = [stat.posclusters.clusterstat];
    if isfield(statRSF,'posclusters') && ~isempty(statRSF.posclusters)
        R.cluster_p_RSF = [statRSF.posclusters.prob];
    else
        R.cluster_p_RSF = [];
    end
    R.null_SMOOTH   = double(stat.posdistribution(:));      % max-cluster-mass nulls, so panel D can be replotted
    R.null_RSF      = double(statRSF.posdistribution(:));
    R.tmap          = double(stat.stat(:));
    R.coverage      = double(stat.coverage(:));

    % ---- (A) coverage ----
    cov = double(stat.coverage); cov(cov < 3) = NaN;
    render(cortex, cov, [0 20], grn, outdir, 'coverage');
    cbar(outdir,'cbar_coverage', grn, [0 20], 'subjects');

    % ---- (B) per-electrode DV berry ----
    berry(cortex, source, [-dvlim dvlim], rdbu, outdir, 'berry');
    cbar(outdir,'cbar_dv', rdbu, [-dvlim dvlim], dvlab);

    % ---- (C) group SMOOTH t-map ----
    render(cortex, double(stat.stat), [-4 4], rdbu, outdir, 'tmap');
    cbar(outdir,'cbar_t', rdbu, [-4 4], 't-value');

    % ---- (D) max-cluster-mass null distribution (SMOOTH vs RSF); genFig5_nullhist redraws it from
    %      genFig5_results.mat in seconds, so axis changes need no re-run ----
    genFig5_nullhist('histogram', outdir, R);

    % ---- (E) PC cumulative variance of the surrogate ensemble ----
    [R.nPC90_SMOOTH, R.nPC90_RSF, R.cumvar_SMOOTH, R.cumvar_RSF] = pc_variance(stat, statRSF, outdir);

    % ---- (E2) reconstruction of the empirical t-map from each ensemble's PCs (RSF supplement) ----
    %      Accuracy curve for every dataset; the example brain panels only for the primary dataset,
    %      since each render is slow.
    recondir = fullfile(paths.figs,'S7',name);   % its own supplementary figure,
    archive_existing(recondir);                                     % not mixed into figs/5
    [R.recon_r, R.recon_k, R.recon_var] = pc_reconstruction(stat, statRSF, cortex, rdbu, recondir, strcmp(name,'primary'));

    % ---- (F) per-subject spatial-autocorrelation scatter (empirical vs surrogate Moran's I) ----
    [R.sa_p, R.sa_t, R.sa_n, R.sa_emp, R.sa_surr, R.sa_surr_sd] = sa_panel(SA.(name), strcmp(name,'primary'), outdir);
    RES.(name) = R;
    if ~isempty(savedir)
        if ~isfolder(savedir), mkdir(savedir); end
        save(fullfile(savedir, sprintf('genFig5_%s_stat.mat', name)), 'stat', 'statRSF', '-v7.3');
    end
    clear stat statRSF                                     % free the surrogate ensembles before the next dataset
end
resfile = fullfile(paths.data,'genFig5_results.mat');
if isfile(resfile)                                            % never clobber a previous harvest
    movefile(resfile, fullfile(paths.data, sprintf('genFig5_results_%s.mat', datestr(now,'yyyymmdd_HHMMSS'))));
end
save(resfile, '-struct', 'RES');
fprintf('\nwrote Figure 5 panels to %s\n', fullfile(paths.figs,'5'));
fprintf('wrote derived numbers to %s\n', resfile);
summarize(RES, DS(:,1));
if all(isfield(RES, DS(:,1)'))                      % cross-dataset panel of Supplementary Figure 7
    genFigS7_tencomponents(10);
end
end

function summarize(RES, names)
% Console digest of everything a caption or Results sentence might need.
fprintf('\n================ genFig5 harvest (NPERM=%d, SEED=%d) ================\n', RES.NPERM, RES.SEED);
for i = 1:numel(names)
    n = names{i}; if ~isfield(RES,n), continue; end
    r = RES.(n);
    sig = r.cluster_p(r.cluster_p <= 0.05);
    fprintf('[%s] clusteralpha=%.3f | %d subjects, %d electrodes\n', n, r.clusteralpha, r.nsub, r.nchan);
    fprintf('   significant cluster p : %s\n', mat2str(round(sig,4)));
    fprintf('   nPC90                 : %d (SMOOTH), %d (RSF)\n', r.nPC90_SMOOTH, r.nPC90_RSF);
    fprintf('   SA paired-t           : p = %.3f (n = %d)\n', r.sa_p, r.sa_n);
    if ~isempty(r.recon_r)
        k = r.recon_k; sel = ismember(k,[1 3 10]);
        fprintf('   recon r @k=1,3,10     : SMOOTH %s | RSF %s\n', ...
                mat2str(round(r.recon_r(sel,1)',3)), mat2str(round(r.recon_r(sel,2)',3)));
    end
end
fprintf('=====================================================================\n');
end

% ==================================== helpers ====================================
function cortex = load_fsaverage_mesh(fshome)
lh = smooth_read_surf(fullfile(fshome,'subjects','fsaverage','surf','lh.pial'));
rh = smooth_read_surf(fullfile(fshome,'subjects','fsaverage','surf','rh.pial'));
cortex.pos = [lh.pos; rh.pos]; cortex.tri = [lh.tri; rh.tri + size(lh.pos,1)]; cortex.unit = lh.unit;
if isfield(lh,'curv') && isfield(rh,'curv'), cortex.curv = [lh.curv; rh.curv];
else, cortex.curv = zeros(size(cortex.pos,1),1); end
end

function archive_existing(d)
% Never overwrite: if the target already holds output, move it to figs/legacy/<name>_<timestamp>/
% before recreating it. Keeps every prior run recoverable.
if isfolder(d)
    e = dir(d); e = e(~ismember({e.name},{'.','..'}));
    if ~isempty(e)
        parts = strsplit(strrep(d,'\\','/'),'/');
        tag = strjoin(parts(max(1,end-1):end),'_');
        dst = fullfile(fileparts(fileparts(d)),'legacy',sprintf('%s_%s',tag,datestr(now,'yyyymmdd_HHMMSS')));
        if ~isfolder(fileparts(dst)), mkdir(fileparts(dst)); end
        movefile(d, dst); fprintf('  archived previous output -> %s\n', dst);
    end
end
if ~isfolder(d), mkdir(d); end
end

function savepanel(outdir,name)
% Transparent-background raster PNG -- the house approach (see genFigS5 / genFigS3). exportgraphics(...,'BackgroundColor','none') silently writes opaque RGB
% for raster output in this R2023a/Rosetta build, so print the same scene on white and on black and
% recover true alpha from the difference. print (not exportgraphics) so both renders come out the
% same fixed size -- exportgraphics auto-crops to non-background content, which would crop the two
% differently and break the differencing.
fw = fullfile(outdir,[name '_w.png']); fkp = fullfile(outdir,[name '_k.png']);
set(gcf,'InvertHardcopy','off','PaperPositionMode','auto'); set(findall(gcf,'type','axes'),'Color','none');
set(gcf,'Color','w'); print(gcf, fw, '-dpng','-r300');
set(gcf,'Color','k'); print(gcf, fkp,'-dpng','-r300'); close(gcf);
W = double(imread(fw)); K = double(imread(fkp));
alpha = min(max(1 - mean(W-K,3)/255, 0), 1);
C = min(max(K./max(alpha,1e-3), 0), 255);
[rr,cc] = find(alpha > 0.02);
if ~isempty(rr), r = min(rr):max(rr); c = min(cc):max(cc); C = C(r,c,:); alpha = alpha(r,c); end
imwrite(uint8(C), fullfile(outdir,[name '.png']), 'Alpha', alpha); delete(fw); delete(fkp);
end

function [R, kvec, Vc] = pc_reconstruction(stat, statRSF, cortex, cmap, outdir, do_maps)
% How well does each surrogate ensemble's own PC basis reconstruct the OBSERVED group t-map?
% RSF surrogates stay closer to the empirical topography, so few RSF PCs already explain the real
% map -- i.e. the RSF null has not broken the structure it is meant to destroy. SMOOTH's ensemble is
% richer and needs far more components, which is the point of the supplement.
%
% Reconstruction uses the ensemble's PC basis only (no refitting): y_hat(k) = B_k (B_k' y), with B_k
% the first k PCs. Because PCA coefficients are orthonormal this is an orthogonal projection, so the
% curve is monotone and directly comparable between ensembles.
% KMAX is the number of components COMPUTED and cached; XMAX is only what the axis shows. The curve
% is accumulated incrementally (yh = yh + b_k*(b_k'y)) instead of rebuilding B_k(B_k'y) at every k,
% making the whole curve O(V*K) rather than O(V*K^2) -- so computing well past the plotted range is
% essentially free, and changing the displayed x-limit later needs NO re-run: just replot from
% genFig5_results.mat (fields recon_r / recon_k). Vc is the fraction of the observed map's variance
% (sum of squares of the centred map) that lies in the first k PCs, |B_k'y|^2 / |y|^2 (field recon_var);
% genFigS7_tencomponents draws its value at k = 10 for all three datasets.
KMAX = 500; XMAX = 250; KSHOW = [1 3 10];
t = double(stat.stat(:));
m = isfinite(t) & ~isnan(stat.surrogate_group(:,1)) & ~isnan(statRSF.surrogate_group(:,1));
y = t(m); y = y - mean(y);

ens = struct('lab',{'SMOOTH','RSF'}, 'E',{stat.surrogate_group(m,:), statRSF.surrogate_group(m,:)}, ...
             'col',{[0.20 0.35 0.65],[0.85 0.35 0.15]});
R = nan(KMAX, numel(ens)); Vc = nan(KMAX, numel(ens));
for e = 1:numel(ens)
    B = pca(double(ens(e).E)');                            % vertices x ncomp, orthonormal columns
    K = min(KMAX, size(B,2));
    yh = zeros(size(y)); cap = 0;
    for k = 1:K
        ck = B(:,k)'*y;
        yh = yh + B(:,k)*ck;                               % cumulative projection onto the first k PCs
        cap = cap + ck^2; Vc(k,e) = cap/(y'*y);            % variance of the observed map captured so far
        R(k,e) = corr(yh, y);
        if do_maps && ismember(k, KSHOW)                   % example reconstruction + residual map
            full = nan(numel(t),1); full(m) = yh;
            render(cortex, full, [-4 4], cmap, outdir, sprintf('recon_%s_k%02d', ens(e).lab, k));
            full(m) = y - yh;                              % residual, same scale as the t-map
            render(cortex, full, [-4 4], cmap, outdir, sprintf('reconerr_%s_k%02d', ens(e).lab, k));
        end
    end
end

f = figure('color','w','position',[100 100 560 450]); hold on; box off; grid on
for e = 1:numel(ens)
    plot(1:KMAX, R(:,e), '-', 'LineWidth',2, 'Color',ens(e).col, 'DisplayName',ens(e).lab);
end
for k = KSHOW, xline(k, ':', 'Color',[.6 .6 .6], 'HandleVisibility','off'); end
xlabel('number of principal components'); ylabel('correlation with empirical t-map');
ylim([0 1]); xlim([1 XMAX]); set(gca,'FontSize',12);
legend('Location','southeast','Box','off');
title('Reconstruction of the observed map from each null ensemble');
exportgraphics(f, fullfile(outdir,'pc_reconstruction.pdf'), 'ContentType','vector'); close(f);
% White-background vector PDF on purpose: this is a 2-D plot, not a brain render. Keeping it
% vector preserves crisp axis text and gridlines at any scale, and matches how every other 2-D
% panel in the repo is exported. Only the brain maps use the transparent savepanel path.

kvec = (1:KMAX)';
fprintf('  reconstruction r at k=%s: SMOOTH %s | RSF %s\n', mat2str(KSHOW), ...
        mat2str(round(R(KSHOW,1)',3)), mat2str(round(R(KSHOW,2)',3)));
fprintf('  variance captured at k=%s: SMOOTH %s | RSF %s\n', mat2str(KSHOW), ...
        mat2str(round(100*Vc(KSHOW,1)',1)), mat2str(round(100*Vc(KSHOW,2)',1)));
end

function [m, keep] = hemimesh(cortex, hemi)
% One hemisphere of the fsaverage mesh (left = first half of the vertices) with its curvature.
nL = size(cortex.pos,1)/2;
if strcmp(hemi,'lh'), keep = (1:nL)'; tri = cortex.tri(all(cortex.tri <= nL, 2), :);
else, keep = (nL+1:2*nL)'; tri = cortex.tri(all(cortex.tri > nL, 2), :) - nL; end
m = struct('pos', cortex.pos(keep,:), 'tri', tri, 'curv', cortex.curv(keep));
end

function render(cortex, y, cl, cmap, outdir, name)
% Lateral view of each hemisphere. Covered vertices are drawn as their own surface over the
% curvature-shaded cortex, so uncovered (NaN) vertices leave no artefacts; saved as a transparent
% PNG (savepanel). Used for the coverage and t-map panels and the Supp Fig 7 reconstructions.
y = y(:);
for v = 1:2
    if v==1, vw=[-90 0]; hemi='lh'; else, vw=[90 0]; hemi='rh'; end
    [m, keep] = hemimesh(cortex, hemi); f = y(keep);
    figure('color','w','position',[100 100 560 500]); axes('Position',[0.02 0.04 0.96 0.9]);
    ft_plot_mesh(m, 'vertexcolor','curv', 'edgecolor','none'); hold on;
    k = isfinite(f);
    if any(k)
        vm = zeros(size(m.pos,1),1); vm(k) = 1:nnz(k); fk = all(k(m.tri),2);
        s.pos = m.pos(k,:); s.tri = vm(m.tri(fk,:));
        ft_plot_mesh(s, 'vertexcolor', f(k), 'edgecolor','none');
    end
    colormap(cmap); clim(cl); view(vw); axis off;
    savepanel(outdir, [name '_' hemi '_lateral']);
end
end

function berry(cortex, source, cl, cmap, outdir, name)
% Electrode values as spheres at each electrode's nearest vertex of the hemisphere shown
% (hemisphere by the sign of x, as in SMOOTHstat); transparent PNG.
for v = 1:2
    if v==1, vw=[-90 0]; hemi='lh'; else, vw=[90 0]; hemi='rh'; end
    m = hemimesh(cortex, hemi); E = []; V = [];
    for i = 1:numel(source)
        cp = double(source{i}.elec.chanpos); dv = double(source{i}.stat(:));
        sel = (cp(:,1) < 0) == strcmp(hemi,'lh');
        E = [E; cp(sel,:)]; V = [V; dv(sel)]; %#ok<AGROW>
    end
    e = m.pos(knnsearch(m.pos, E), :);
    figure('color','w','position',[100 100 560 500]); axes('Position',[0.02 0.04 0.96 0.9]);
    ft_plot_mesh(m, 'vertexcolor','curv', 'edgecolor','none'); hold on;
    ft_plot_cloud(e, V, 'cloudtype','surf','radius',2.5,'scalerad','no','colormap',cmap,'clim',cl);
    view(vw); axis off;
    savepanel(outdir, [name '_' hemi '_lateral']);
end
end

function cbar(outdir,name,cmap,cl,label)
f = figure('color','w','position',[100 100 150 430]);
ax = axes('Parent',f,'Position',[0.05 0.06 0.02 0.88]); axis(ax,'off'); colormap(ax,cmap); set(ax,'CLim',cl);
cb = colorbar(ax,'Position',[0.32 0.06 0.22 0.88]); cb.Limits = cl; set(cb,'FontSize',14); ylabel(cb,label,'FontSize',14);
exportgraphics(f,fullfile(outdir,[name '.pdf']),'ContentType','vector'); close(f);
end

function [nZ, nY, cumZ, cumY] = pc_variance(stat, statRSF, outdir)
mZ = ~isnan(stat.surrogate_group(:,1));    Z = stat.surrogate_group(mZ,:);    [~,~,~,~,eZ] = pca(Z');
mY = ~isnan(statRSF.surrogate_group(:,1)); Y = statRSF.surrogate_group(mY,:); [~,~,~,~,eY] = pca(Y');
cumZ = cumsum(eZ); cumY = cumsum(eY);
nZ = find(cumZ >= 90, 1); nY = find(cumY >= 90, 1);
nZ = def(nZ); nY = def(nY);            % NaN rather than [] so the harvest/summary can't break
f = figure('color','w','position',[100 100 520 440]); hold on; box off; grid on
plot(cumZ, 'LineWidth',2, 'Color',[0.18 0.55 0.34]); plot(cumY, 'LineWidth',2, 'Color',[0.85 0.37 0.05]);
yline(90,':','90%'); if isfinite(nZ), xline(nZ,'-','Color',[0.18 0.55 0.34]); end
if isfinite(nY), xline(nY,'--','Color',[0.85 0.37 0.05]); end
xlim([0 numel(cumZ)]); ylim([0 100]); set(gca,'FontSize',12);
xlabel('PC'); ylabel('% variance explained');
legend({'SMOOTH','RSF'}, 'Box','off','Location','southeast');
title(sprintf('nPC90 = %d (SMOOTH), %d (RSF)', nZ, nY));
exportgraphics(f, fullfile(outdir,'pc_variance.pdf'), 'ContentType','vector'); close(f);
end

function [pp, tt, n, emp, nmu, nsd] = sa_panel(M, primary, outdir)
% Per-participant spatial autocorrelation (Moran's I) of the empirical map vs the mean over that
% participant's SMOOTH surrogates (horizontal bars: surrogate SD), from SMOOTHdiag values saved in
% paper/data/fig5_sa_moran.mat. Axes zoom to the data plus an 8% pad and are labelled at the ends only.
% For the validation datasets, participants whose surrogate Moran's I lies below the largest gap (> 0.02)
% in the lower half of the surrogate values are drawn clamped at the left edge with their value.
emp = M.emp(:); nmu = M.null_mu(:); nsd = M.null_sd(:);
ok = isfinite(emp) & isfinite(nmu); emp = emp(ok); nmu = nmu(ok); nsd = nsd(ok); n = nnz(ok);
[~,pp,~,ts] = ttest(emp, nmu); tt = ts.tstat;
inb = true(size(nmu));
if ~primary
    sv = sort(nmu); g = diff(sv); g(sv(1:end-1) >= median(nmu)) = 0;
    if ~isempty(g)
        [maxgap, gi] = max(g);
        if maxgap > 0.02, inb = nmu >= (sv(gi)+sv(gi+1))/2; end
    end
end
sd0 = nsd; sd0(isnan(sd0)) = 0;
bulk = [emp(inb); nmu(inb)-sd0(inb); nmu(inb)+sd0(inb)]; pad = 0.08*range(bulk);
lim = [min(bulk)-pad, max(bulk)+pad];
f = figure('color','w','units','centimeters','Position',[0 0 11 11]); hold on
plot(lim, lim, '-', 'Color',[0.7 0.7 0.7], 'LineWidth', 1 + 0.5*primary);
errorbar(nmu(inb), emp(inb), nsd(inb), 'horizontal', 'LineStyle','none', 'Color',[0.6 0.6 0.6], 'CapSize',0, 'LineWidth',0.5);
scatter(nmu(inb), emp(inb), 100, 'filled', 'MarkerFaceColor',[53 154 71]/255, 'MarkerFaceAlpha',0.80, 'MarkerEdgeColor','none');
C_off = [0.78 0.20 0.15]; off = find(~inb); r = range(lim);
for j = 1:numel(off)
    k = off(j); yc = min(max(emp(k),lim(1)),lim(2));
    plot([lim(1)+0.16*r, lim(1)+0.03*r], [yc yc], '-', 'Color',C_off, 'LineWidth',1.2, 'Clipping','off');
    plot(lim(1)+0.03*r, yc, '<', 'MarkerSize',8, 'MarkerFaceColor',C_off, 'MarkerEdgeColor','none', 'Clipping','off');
    text(lim(1)+0.18*r, yc, sprintf('surr=%.2f', nmu(k)), 'Color',C_off, 'FontSize',8, 'VerticalAlignment','middle');
end
if primary, xlabel('surrogate Moran''s I (subject mean)'); else, xlabel('surrogate Moran''s I'); end
ylabel('empirical Moran''s I');
ttl = sprintf('t(%d) = %.2f, p = %.3g', ts.df, ts.tstat, pp);
if ~isempty(off), ttl = sprintf('%s  (%d off-scale)', ttl, numel(off)); end
title(ttl, 'FontWeight','normal');
axis square; xlim(lim); ylim(lim); box off;
tk = linspace(lim(1), lim(2), 6); xticks(tk); yticks(tk);
lab = repmat({''},1,6); lab{1} = num2str(lim(1),'%.3f'); lab{6} = num2str(lim(2),'%.3f');
xticklabels(lab); yticklabels(lab);
grid on; set(gca,'FontSize', 11 + primary, 'TickDir','in');
exportgraphics(f, fullfile(outdir,'sa_scatter.pdf'), 'ContentType','vector'); close(f);
end

function v = def(x), if isempty(x), v = NaN; else, v = x; end, end
