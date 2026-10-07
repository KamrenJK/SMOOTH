function export_smooth_sigmaps(NPERM)
% ------------------------------------------------------------------------------------------------
% Cache the SMOOTH group-level result for all three datasets so genFigS6 can add a
% SMOOTH row alongside the channel / DK / Schaefer rows.
%
% cfg is copied verbatim from genFig5 (kernelwidth 12, graphsigma 8, sphere/gaussian, rankrescale
% 'exact', tail 'positive', seed 42) so the clusters shown here are the SAME clusters as Figure 5
% rather than an independently-tuned run.
%
% What is saved is the FWE-corrected SIGNIFICANT cluster t-map (stat.stat masked by stat.mask), which
% is the correct analogue of the other rows: those show FDR-thresholded channels/parcels, so SMOOTH
% must likewise be shown thresholded, not as a raw t-map.
%
% SMOOTH runs on full fsaverage (327684 vertices) but the figure renders on fsaverage5 (20484).
% The transfer is INDEX-based: fsaverage5 vertex i is fsaverage vertex i, by construction of the
% icosahedral hierarchy (ic5 vertices are the first 10242 of ic7, per hemisphere). Their pial
% COORDINATES differ slightly (mean 1.0 mm, max 5.8 mm) because each resolution's pial is
% reconstructed separately -- that is a geometry difference, not a correspondence failure. Do NOT
% substitute a spatial nearest-neighbour match here: it looks reasonable (max ~3.9 mm) but can map a
% vertex onto the opposite bank of a sulcus, which is spatially near and anatomically wrong.
%
% Usage:  export_smooth_sigmaps        % NPERM = 1000, matches genFig5
%         export_smooth_sigmaps(100)   % quick preview
% Writes: paper/data/smooth_sigmaps.mat
% ------------------------------------------------------------------------------------------------
if nargin < 1 || isempty(NPERM), NPERM = 1000; end
SEED = 42;
addpath(fullfile(fileparts(mfilename('fullpath')),'..','..'));
paths = smooth_setup('quiet');

% dataset key must match the naming used in methodscomparison_pvals.mat ('socialcog', not 'primary')
%
% CLUSTER-FORMING THRESHOLD IS PER-DATASET and is stated explicitly here rather than left to the
% SMOOTHstat default, because it is not uniform across datasets and a silent default would
% misreproduce Figure 5. Berezutskaya uses 0.01 (z=2.33): at the 0.05 default the whole left
% temporal-frontal territory merges into ONE 12,691-vertex cluster (p=0.001), whereas 0.01 resolves
% the three anatomically distinct clusters the manuscript reports (all p=0.001; lh superior/middle
% temporal, lh precentral/parsopercularis, rh insula), from a sweep of the threshold over
% 0.010 / 0.025 / 0.050 / 0.100.
DS = { 'socialcog',    fullfile(paths.demo,'source.mat'),              0.05;
       'berezutskaya', fullfile(paths.data,'source_berezutskaya.mat'), 0.01;
       'cogitate',     fullfile(paths.data,'source_cogitate.mat'),     0.05 };

M   = load(fullfile(paths.data,'methodscomparison_pvals.mat'), 'fs5_pos', 'nL');
fs5 = double(M.fs5_pos); nL5 = double(M.nL);

lh = ft_read_headshape(fullfile(paths.freesurfer,'subjects','fsaverage','surf','lh.pial'));
rh = ft_read_headshape(fullfile(paths.freesurfer,'subjects','fsaverage','surf','rh.pial'));
nLf = size(lh.pos,1);
map = [(1:nL5)'; nLf + (1:nL5)'];          % index correspondence, per hemisphere

% sanity-check the correspondence rather than trusting it: fsaverage5's own pial at vertex i should
% be ~1 mm from fsaverage's pial at vertex i (same anatomical vertex, different reconstruction).
% A wrong correspondence would show tens of mm.
lh5 = ft_read_headshape(fullfile(paths.freesurfer,'subjects','fsaverage5','surf','lh.pial'));
rh5 = ft_read_headshape(fullfile(paths.freesurfer,'subjects','fsaverage5','surf','rh.pial'));
dIdx = [sqrt(sum((double(lh.pos(1:nL5,:)) - double(lh5.pos)).^2, 2)); ...
        sqrt(sum((double(rh.pos(1:nL5,:)) - double(rh5.pos)).^2, 2))];
fprintf('fs5 <- fsaverage index map: pial offset mean %.3g mm, max %.3g mm (expect ~1 / ~6)\n', ...
        mean(dIdx), max(dIdx));
assert(mean(dIdx) < 3, 'export_smooth_sigmaps:badIndexMap', ...
    'index correspondence looks wrong (mean pial offset %.3g mm)', mean(dIdx));
assert(max(abs(double(lh5.pos) - fs5(1:nL5,:)), [], 'all') < 1e-6, ...
    'export_smooth_sigmaps:meshMismatch', 'fs5_pos is not fsaverage5 lh/rh.pial');

out = struct();
for di = 1:size(DS,1)
    name = DS{di,1};
    L = load(DS{di,2}); source = L.source;
    fprintf('\n[%s] %d subjects, %d electrodes | NPERM=%d\n', name, numel(source), ...
            sum(cellfun(@(s) numel(s.stat), source)), NPERM);

    cfg = []; cfg.fshome = paths.freesurfer; cfg.numrandomization = NPERM;
    cfg.kernelwidth = 12; cfg.graphsigma = 8; cfg.smooth = 'sphere'; cfg.kernel = 'gaussian';
    cfg.rankrescale = 'exact'; cfg.tail = 'positive'; cfg.randomseed = SEED;
    cfg.clusteralpha = DS{di,3};
    cfg.keepmaps = 'no'; cfg.keepsubsurrogates = 'no'; cfg.keepgroupsurrogates = 'no'; cfg.normalize = 'no';
    fprintf('  clusteralpha = %.3f (z = %.2f)\n', cfg.clusteralpha, norminv(1-cfg.clusteralpha));
    stat = SMOOTHstat(cfg, source{:});

    t = double(stat.stat); msk = logical(stat.mask);
    sig = nan(size(t)); sig(msk) = t(msk);

    out.([name '_smooth_sigt']) = reshape(sig(map), 1, []);   % 1 x 20484, NaN where not significant
    out.([name '_smooth_tmap']) = reshape(t(map),   1, []);   % unthresholded, kept for reference
    % full-resolution copies so the fs5 transfer can be revisited without re-running SMOOTH
    out.([name '_smooth_sigt_fsavg']) = reshape(sig, 1, []);
    out.([name '_smooth_tmap_fsavg']) = reshape(t,   1, []);
    out.([name '_clusteralpha'])      = DS{di,3};
    if isempty(stat.posclusters), pv = []; else, pv = round([stat.posclusters.prob], 4); end
    fprintf('  cluster p-values: %s | %d/%d fs5 vertices significant\n', ...
            mat2str(pv), nnz(isfinite(out.([name '_smooth_sigt']))), numel(map));
    clear stat
end
out.NPERM = NPERM; out.SEED = SEED;
save(fullfile(paths.data,'smooth_sigmaps.mat'), '-struct', 'out');
fprintf('\nwrote %s\n', fullfile(paths.data,'smooth_sigmaps.mat'));
end
