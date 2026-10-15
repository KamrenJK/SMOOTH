function export_fig5_sa(diagdir)
% ------------------------------------------------------------------------------------------------
% Cache the per-participant spatial autocorrelation used by the Figure 5 SA panels: Moran's I (both
% hemispheres) of each participant's empirical map, and the mean and SD of Moran's I over that
% participant's SMOOTH surrogate maps, as computed by SMOOTHdiag.
%
% cfg is copied from genFig5 (kernelwidth 12, graphsigma 8, sphere/gaussian, rankrescale 'exact',
% seed 42). SMOOTHdiag uses the first 100 surrogates of each participant (cfg.nsurrsamples = 100), so
% only 100 permutations are run: with the same seed, these are exactly the first 100 surrogates of
% genFig5's 1000-permutation run. Moran's I does not depend on the cluster-forming threshold.
%
% Memory: the kept surrogates take participants x 327,684 vertices x 100 x 8 bytes (about 10-13 GB
% for these datasets) and SMOOTHdiag needs roughly twice that. Datasets are run one at a time.
%
% Usage:  export_fig5_sa            % writes paper/data/fig5_sa_moran.mat (read by genFig5)
%         export_fig5_sa(diagdir)   % also saves each dataset's full SMOOTHdiag output to diagdir
%                                   % (of this 100-permutation run; only its subject-level blocks
%                                   % correspond to genFig5's 1000-permutation run)
% ------------------------------------------------------------------------------------------------
if nargin < 1, diagdir = ''; end
SEED = 42; NSURR = 100;
addpath(fullfile(fileparts(mfilename('fullpath')),'..','..'));
paths = smooth_setup('quiet');

DS = { 'primary',      fullfile(paths.demo,'source.mat');
       'berezutskaya', fullfile(paths.data,'source_berezutskaya.mat');
       'cogitate',     fullfile(paths.data,'source_cogitate.mat') };

OUT = struct();
for di = 1:size(DS,1)
    name = DS{di,1}; L = load(DS{di,2}); source = L.source;
    cfg = []; cfg.fshome = paths.freesurfer; cfg.numrandomization = NSURR;
    cfg.kernelwidth = 12; cfg.graphsigma = 8; cfg.smooth = 'sphere'; cfg.kernel = 'gaussian';
    cfg.rankrescale = 'exact'; cfg.tail = 'positive'; cfg.randomseed = SEED; cfg.normalize = 'no';
    cfg.keepmaps = 'yes'; cfg.keepsubsurrogates = 'yes'; cfg.subsamplesurr = 'yes';
    cfg.nsurrsamples = NSURR; cfg.keepgroupsurrogates = 'no'; cfg.feedback = 'no';
    fprintf('\n[%s] %d participants | %d surrogates each\n', name, numel(source), NSURR);
    stat = SMOOTHstat(cfg, source{:});
    d    = SMOOTHdiag(stat);
    clear stat                                   % free the surrogates before the next dataset
    p = d.subject.spatial.BOTH.moranI.per_subject;
    OUT.(name) = struct('emp', double(p.emp(:)), 'null_mu', double(p.null_mu(:)), 'null_sd', double(p.null_sd(:)));
    [~, pv, ~, ts] = ttest(OUT.(name).emp, OUT.(name).null_mu);
    fprintf('[%s] empirical vs surrogate Moran''s I: t(%d) = %.2f, p = %.3f\n', name, ts.df, ts.tstat, pv);
    if ~isempty(diagdir)
        if ~isfolder(diagdir), mkdir(diagdir); end
        save(fullfile(diagdir, sprintf('fig5_diag_%s.mat', name)), 'd', '-v7.3');
    end
    clear d
end
OUT.description = ['Per-participant spatial autocorrelation (Moran''s I, both hemispheres) of the ' ...
    'empirical map (emp) and the mean (null_mu) and SD (null_sd) over that participant''s first 100 ' ...
    'SMOOTH surrogate maps, from SMOOTHdiag (export_fig5_sa.m). Drawn as the Figure 5d-f SA scatter by genFig5.m.'];
OUT.generated_by = 'export_fig5_sa.m (surfaces read with smooth_read_surf)';
dst = fullfile(paths.data, 'fig5_sa_moran.mat');
save(dst, '-struct', 'OUT');
fprintf('\nwrote %s\n', dst);
end
