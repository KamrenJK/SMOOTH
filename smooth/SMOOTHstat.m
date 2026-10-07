function stat = SMOOTHstat(cfg, varargin)

% SMOOTHstat - Surrogate map-based group-level analysis of human intracranial EEG data
%
% Use as:
%   stat = SMOOTHstat(cfg, source1, source2, ...)
%
% where source1..N are data structures with fields
%   source.stat         = dependent variable per channel
%   source.label        = channel labels
%   source.elec         = electrode struct with .chanpos in fsavg coords
%
% cfg fields
%   cfg.fshome              = FreeSurfer home (required)
%   cfg.parameter           = field with DV (default 'stat')
%   cfg.smooth              = 'sphere' (default) | 'pial'
%                             'sphere' measures distances on fsaverage sphere.reg
%                             (fast, unfolds cortex so sulcal banks are correctly
%                             separated, but distances are NOT proportional to mm --
%                             see cfg.kernelwidth note below).
%                             'pial' uses true geodesic distances in mm (5-10x slower).
%   cfg.kernel              = 'gaussian' (default) | 'boxcar'
%   cfg.kernelwidth         = FWHM (gaussian) or radius (boxcar) in mm (default 12)
%   cfg.graphsigma          = Laplacian kernel σ for spectral surrogates in mm (default 8)
%   cfg.rankrescale         = 'no' (default) | 'exact' | 'quantile' | 'zscore'
%                             marginal distribution matching of surrogates
%   cfg.equalize_n          = 'no' | 'yes'  (downweight high coverage in t)
%   cfg.normalize           = 'no' | 'yes'  (per-subject z before t)
%   cfg.randmethod          = 'signflip' (default) | 'bandrotate'
%   cfg.numrandomization    = permutations (default 1000)
%   cfg.randomseed          = [] (default, nonreproducible) | integer | 'shuffle'
%                             Seeds the RNG so surrogates -- and hence p-values --
%                             are exactly reproducible. The seed actually used is
%                             always returned in stat.cfg.randomseed.
%   cfg.minnbsub            = mininum number of subjects per vertex (default 3)
%   cfg.clusteralpha        = vertex-level cluster-FORMING α per tail (default 0.05).
%                             Sets the threshold a vertex must exceed to join a
%                             cluster. Larger -> lower bar -> bigger, more diffuse
%                             clusters; smaller -> higher bar -> smaller, more focal
%                             ones. It does NOT control which clusters are declared
%                             significant -- that is cfg.alpha.
%   cfg.alpha               = cluster-level SIGNIFICANCE α (default 0.05). Applied to
%                             the permutation p-value of each cluster, halved when
%                             cfg.tail='both'.
%   cfg.tail                = 'both' (default) | 'positive' | 'negative'
%                             which tail stat.mask reflects, and the α it uses:
%                             one-tailed -> cfg.alpha; 'both' -> cfg.alpha/2
%                             (Bonferroni across the two tails). All clusters are
%                             returned in stat.posclusters/negclusters regardless.
%
%                             The default is two-tailed because it is the safe choice
%                             on unfamiliar data: iEEG task contrasts frequently show
%                             widespread suppression as well as activation (in the
%                             bundled demo dataset 69% of channel-level statistics are
%                             negative), and a one-tailed default would silently hide
%                             half the effects. Set cfg.tail='positive' when the
%                             hypothesis is directional -- this is what the SMOOTH
%                             manuscript reports, so paper/code/genFig*.m set it
%                             explicitly.
%   cfg.keepmaps            = 'no' | 'yes' (default 'no')
%   cfg.keepsubsurrogates   = 'no' | 'yes' (memory expensive, default 'no')
%   cfg.keepgroupsurrogates = 'no' | 'yes' (memory expensive, default 'no')
%   cfg.subsamplesurr       = 'no' | 'yes' (randomly sample 100 surrogates to keep to limit RAM demands, default 'yes')
%   cfg.nsurrsamples        = how many surrogates to save (default 100);
%   cfg.feedback            = 'text' (default) | 'no'  per-permutation progress
%
% Output (stat)
%   stat.stat            = group-level t-map on fsaverage vertices
%   stat.prob            = vertex-wise cluster p-value, ONE-TAILED (not corrected
%                          across tails). Use stat.mask for thresholded inference;
%                          do not threshold stat.prob directly without accounting
%                          for cfg.tail.
%   stat.mask            = significance mask, FWE-controlled per cfg.tail
%   stat.posclusters     = positive clusters, each with .clusterstat and .prob
%   stat.negclusters     = negative clusters, each with .clusterstat and .prob
%   stat.coverage        = per-vertex subject count
%   stat.cfg             = configuration actually used (incl. resolved randomseed)
%
% Requirements: FreeSurfer fsaverage, FieldTrip
% Authors: Kamren Khan & Arjen Stolk, v0.1


% ---- Defaults
cfg                     = ft_checkconfig(cfg, 'required', {'fshome'});
cfg.parameter           = ft_getopt(cfg, 'parameter',         'stat');
cfg.smooth              = ft_getopt(cfg, 'smooth',          'sphere');
cfg.kernel              = ft_getopt(cfg, 'kernel',        'gaussian');
cfg.kernelwidth         = ft_getopt(cfg, 'kernelwidth',           12); % mm FWHM; informed by the HFB autocorrelation analysis
cfg.graphsigma          = ft_getopt(cfg, 'graphsigma',             8); % mm
cfg.rankrescale         = ft_getopt(cfg, 'rankrescale',         'no');
cfg.equalize_n          = ft_getopt(cfg, 'equalize_n',          'no');
cfg.normalize           = ft_getopt(cfg, 'normalize',           'no');
cfg.randmethod          = ft_getopt(cfg, 'randmethod',    'signflip');
cfg.numrandomization    = ft_getopt(cfg, 'numrandomization',    1000);
cfg.randomseed          = ft_getopt(cfg, 'randomseed',           []);
cfg.minnbsub            = ft_getopt(cfg, 'minnbsub',               3);
cfg.clusteralpha        = ft_getopt(cfg, 'clusteralpha',        0.05);
cfg.alpha               = ft_getopt(cfg, 'alpha',               0.05);
cfg.tail                = ft_getopt(cfg, 'tail',             'both');
cfg.keepmaps            = ft_getopt(cfg, 'keepmaps',            'no');
cfg.keepsubsurrogates   = ft_getopt(cfg, 'keepsubsurrogates',   'no');
cfg.keepgroupsurrogates = ft_getopt(cfg, 'keepgroupsurrogates', 'no');
cfg.subsamplesurr       = ft_getopt(cfg, 'subsamplesurr',       'yes');
cfg.nsurrsamples        = ft_getopt(cfg, 'nsurrsamples',         100);
cfg.feedback            = ft_getopt(cfg, 'feedback',          'text');
% Geometry cache. The fsaverage meshes, the electrode->vertex projectors and the
% graph eigenbases depend only on electrode POSITIONS, never on the dependent
% variable. When the same montage is analysed repeatedly (simulations, sensitivity
% analyses, bootstraps) that work can be done once and reused: pass
% cfg.returncache='yes' to get stat.cache back, then hand it in as cfg.cache on
% subsequent calls. Results are identical either way -- the cache only skips
% recomputation. Leave both unset for normal single-analysis use.
cfg.cache               = ft_getopt(cfg, 'cache',                 []);
cfg.returncache         = ft_getopt(cfg, 'returncache',         'no');

% ---- Backward compatibility: cfg.keepsurrogates was split into separate
% subject-level and group-level flags. Honour the old name rather than letting
% it be silently ignored (which would return no surrogates at all).
if isfield(cfg, 'keepsurrogates')
  warning('SMOOTHstat:deprecatedOption', ...
    ['cfg.keepsurrogates is deprecated; it has been split into ' ...
     'cfg.keepsubsurrogates and cfg.keepgroupsurrogates. Applying ''%s'' to both.'], ...
    cfg.keepsurrogates);
  cfg.keepsubsurrogates   = cfg.keepsurrogates;
  cfg.keepgroupsurrogates = cfg.keepsurrogates;
  cfg = rmfield(cfg, 'keepsurrogates');
end

% ---- Reproducibility: seed the RNG and record what was actually used
if ~isempty(cfg.randomseed)
  rng(cfg.randomseed);
else
  % capture the arbitrary state so the run can be replayed even if unseeded
  s = rng();
  cfg.randomseed = s;
  warning('SMOOTHstat:unseeded', ...
    ['cfg.randomseed not set: p-values will differ between runs. ' ...
     'The RNG state used has been stored in stat.cfg.randomseed. ' ...
     'Set cfg.randomseed to an integer for reproducible results.']);
end

% ---- Validate inputs before doing any expensive work
validate_sources(cfg, varargin{:});

% ---- A cache is only valid for the montage and geometry it was built from.
% Reusing one across a different electrode set or smoothing kernel would silently
% project the data through the wrong operators, so refuse rather than guess.
if ~isempty(cfg.cache)
  if ~isfield(cfg.cache, 'signature')
    error('SMOOTHstat:badCache', ...
      'cfg.cache has no signature field. Rebuild it with cfg.returncache=''yes''.');
  end
  sig = cfg.cache.signature;
  now_nchan = cellfun(@(d) numel(d.(cfg.parameter)), varargin(:)');
  chk = {'kernel',cfg.kernel; 'smooth',cfg.smooth};
  for c = 1:size(chk,1)
    if ~strcmpi(sig.(chk{c,1}), chk{c,2})
      error('SMOOTHstat:cacheMismatch', ...
        'cfg.cache was built with cfg.%s=''%s'' but this call uses ''%s''. Rebuild the cache.', ...
        chk{c,1}, sig.(chk{c,1}), chk{c,2});
    end
  end
  if sig.kernelwidth ~= cfg.kernelwidth || sig.graphsigma ~= cfg.graphsigma
    error('SMOOTHstat:cacheMismatch', ...
      ['cfg.cache was built with kernelwidth=%g, graphsigma=%g but this call uses ' ...
       '%g and %g. Rebuild the cache.'], ...
      sig.kernelwidth, sig.graphsigma, cfg.kernelwidth, cfg.graphsigma);
  end
  if sig.nSubjects ~= numel(varargin) || ~isequal(sig.nchan(:), now_nchan(:))
    error('SMOOTHstat:cacheMismatch', ...
      ['cfg.cache was built for a different montage (%d subjects, %d channels) than ' ...
       'this call (%d subjects, %d channels). Rebuild the cache.'], ...
      sig.nSubjects, sum(sig.nchan), numel(varargin), sum(now_nchan));
  end
end

% ---- Meshes (fsaverage) and adjacency (for clustering)
usecache = ~isempty(cfg.cache);
if usecache && isfield(cfg.cache, 'geom')
  g     = cfg.cache.geom;
  lh    = g.lh;    rh    = g.rh;
  sph_l = g.sph_l; sph_r = g.sph_r;
  sph   = g.sph;   adj   = g.adj;
  nvtx  = length(lh.pos);
else
  lh      = ft_read_headshape(fullfile(cfg.fshome,'subjects','fsaverage','surf','lh.pial'));
  rh      = ft_read_headshape(fullfile(cfg.fshome,'subjects','fsaverage','surf','rh.pial'));
  nvtx    = length(lh.pos);  % resolution
  sph_l   = ft_read_headshape(fullfile(cfg.fshome,'subjects','fsaverage','surf','lh.sphere.reg'));
  sph_r   = ft_read_headshape(fullfile(cfg.fshome,'subjects','fsaverage','surf','rh.sphere.reg'));
  sph.pos = [sph_l.pos; sph_r.pos];
  sph.tri = [sph_l.tri; sph_r.tri + size(sph_l.pos,1)];
  i = [sph.tri(:,1); sph.tri(:,1); sph.tri(:,2); sph.tri(:,2); sph.tri(:,3); sph.tri(:,3)];
  j = [sph.tri(:,2); sph.tri(:,3); sph.tri(:,1); sph.tri(:,3); sph.tri(:,1); sph.tri(:,2)];
  adj = sparse(i,j,1,size(sph.pos,1),size(sph.pos,1)); adj = adj | adj'; % undirected
end

% -------------------------------------------------------------------------
% Step 1: Generation of Cortical Activation Maps
% -------------------------------------------------------------------------
fprintf('generating empirical maps\n')

nSubjects = numel(varargin);
nVertL  = size(sph_l.pos,1);
nVertR  = size(sph_r.pos,1);
allmaps = nan(nVertL+nVertR, nSubjects);

% convert FWHM (mm) to σ (mm) for Gaussian kernel (if used)
if strcmpi(cfg.kernel, 'gaussian')
  sigma_mm = cfg.kernelwidth/(2*sqrt(2*log(2)));   % FWHM -> σ
elseif strcmpi(cfg.kernel, 'boxcar')
  sigma_mm = cfg.kernelwidth;                      % boxcar radius
else
  error('cfg.kernel must be ''gaussian'' or ''boxcar''.');
end

% one cache reused across subjects and hemispheres
e2m_cache = struct('pial', struct('L', [], 'R', []), ...
  'sphere', struct('kdtL', [], 'kdtR', []));

for s = 1:nSubjects
  data = varargin{s};
  dv   = data.(cfg.parameter);
  pos  = data.elec.chanpos;
  hemi = pos(:,1) < 0;  % fsavg: L<0

  % electrode to vertex smoothing. The projector depends only on electrode
  % positions, so when a cache is supplied it is reused and only the (cheap)
  % apply is redone against the new dependent variable.
  if usecache && numel(cfg.cache.subject) >= s && ~isempty(cfg.cache.subject(s).smooth_cache)
    cacheL = cfg.cache.subject(s).smooth_cache.L;
    cacheR = cfg.cache.subject(s).smooth_cache.R;
    idxL   = cfg.cache.subject(s).idxL;
    idxR   = cfg.cache.subject(s).idxR;
    mapL   = elec2map_apply(cacheL, dv( hemi));
    mapR   = elec2map_apply(cacheR, dv(~hemi));
  else
    [mapL, cacheL, idxL, e2m_cache] = elec2map_build('L', cfg, lh, sph_l, pos, dv, sigma_mm,  hemi, e2m_cache);
    [mapR, cacheR, idxR, e2m_cache] = elec2map_build('R', cfg, rh, sph_r, pos, dv, sigma_mm, ~hemi, e2m_cache);
  end

  allmaps(1:nVertL,               s) = mapL;
  allmaps(nVertL+1:nVertL+nVertR, s) = mapR;

  % rank-rescale prep
  rescale_cache = struct();
  m = isfinite(mapL);
  rescale_cache.L = struct('mask',m,'tgt_sorted',sort(mapL(m),'ascend'), ...
    'mu',mean(mapL(m)), 'sd',max(std(mapL(m)),eps));
  m = isfinite(mapR);
  rescale_cache.R = struct('mask',m,'tgt_sorted',sort(mapR(m),'ascend'), ...
    'mu',mean(mapR(m)), 'sd',max(std(mapR(m)),eps));

  % keep for permutations
  data.elec.hemi     = hemi;
  data.elec.idxL     = idxL;
  data.elec.idxR     = idxR;
  data.smooth_cache  = struct('L',cacheL,'R',cacheR);
  data.rescale_cache = rescale_cache;
  varargin{s} = data;
end

% -------------------------------------------------------------------------
% Step 2: Surrogate Map Generation via Spectral Resampling
% -------------------------------------------------------------------------
fprintf('generating surrogate maps\n')

for s = 1:nSubjects
  data = varargin{s};

  % electrode sphere coords (arbitrary units) for eigenbasis
  Lh_sph = sph_l.pos(data.elec.idxL, :);
  Rh_sph = sph_r.pos(data.elec.idxR, :);

  % (a)+(b): Graph construction + Eigenbasis (per hemisphere). Depends only on the
  % electrode sphere coordinates and cfg.graphsigma, so it is cacheable.
  if usecache && numel(cfg.cache.subject) >= s && ~isempty(cfg.cache.subject(s).Q)
    Q   = cfg.cache.subject(s).Q;
    lam = cfg.cache.subject(s).eigvals;
  else
    [QL, lamL] = graph_eigenbasis(Lh_sph, mm2sphere(cfg.graphsigma));
    [QR, lamR] = graph_eigenbasis(Rh_sph, mm2sphere(cfg.graphsigma));
    Q   = blkdiag(QL, QR); % (nL+nR) x (nL+nR)
    lam = [lamL; lamR];
  end

  % (c) Projection (modal coefficients)
  y                = data.(cfg.parameter);
  y                = [y(data.elec.hemi); y(~data.elec.hemi)];  % reorder;
  data.modes       = Q;                         % elecs x modes
  % data.modalcoeffs = Q' * data.(cfg.parameter); % modal coefficients
  data.modalcoeffs = Q' * y;
  data.eigvals     = lam;                       % eigenvalues
  data.B           = data.modes * spdiags(data.modalcoeffs,0,length(data.modalcoeffs),length(data.modalcoeffs));
  data.elec.nL     = sum(data.elec.hemi);
  varargin{s} = data;
end

% -------------------------------------------------------------------------
% Steps 3: Group-Level Cluster Statistics
% -------------------------------------------------------------------------
fprintf('performing cluster-based permutation testing\n')

coverage = sum(isfinite(allmaps), 2);  % V × 1
mask     = coverage < cfg.minnbsub; % vertices below min subjects excluded

negtailcritval = norminv(cfg.clusteralpha);
postailcritval = norminv(1 - cfg.clusteralpha);
posdistribution = zeros(cfg.numrandomization,1);
negdistribution = zeros(cfg.numrandomization,1);

if strcmp(cfg.keepsubsurrogates, 'yes')
    if strcmp(cfg.subsamplesurr, 'yes')
      nkeep = cfg.nsurrsamples;
    else
      nkeep = cfg.numrandomization;
    end
    gb = nSubjects * size(sph.pos,1) * nkeep * 8 / 1024^3;
    if gb > 4
      warning('SMOOTHstat:memory', ...
        ['cfg.keepsubsurrogates=''yes'' will allocate %.1f GB (%d subjects x %d vertices x %d surrogates). ' ...
         'Reduce cfg.nsurrsamples, or set cfg.subsamplesurr=''yes'', or keep only group surrogates.'], ...
        gb, nSubjects, size(sph.pos,1), nkeep);
    else
      fprintf('  allocating %.1f GB for subject-level surrogates\n', gb);
    end
    surr_subj = zeros(nSubjects, size(sph.pos,1), nkeep);
end

if strcmp(cfg.keepgroupsurrogates, 'yes')
  fprintf('  allocating %.1f GB for group-level surrogates\n', ...
    size(sph.pos,1) * cfg.numrandomization * 8 / 1024^3);
  surr_group = zeros(size(sph.pos,1), cfg.numrandomization);
end

% main loop: permutations + one empirical pass
dofeedback = ~strcmp(cfg.feedback, 'no');
for iter = 1:cfg.numrandomization+1
  if iter <= cfg.numrandomization
    if dofeedback && (iter==1 || mod(iter,50)==0 || iter==cfg.numrandomization)
      fprintf(' permutation %d/%d\n', iter, cfg.numrandomization);
    end

    permMaps = nan(size(allmaps)); % V × S

    for s = 1:nSubjects
      data = varargin{s};

      % (d) Resampling: spectral surrogate in electrode space
      switch lower(cfg.randmethod)
        case 'signflip' % random sign flips
          sflip = sign(randn(size(data.modalcoeffs)));
          dvn   = data.B * sflip;
        case 'bandrotate' % band-oriented orthogonal rotation
          a2    = spectral_bandrotate(data.modalcoeffs, data.eigvals);
          dvn   = data.modes * a2;  % reconstruct
        otherwise
          error('cfg.randmethod must be ''signflip'' or ''bandrotate''.');
      end

      % (e) Postprocessing: apply same smoothing as empirical
      % dvn is in [L; R] block order (matching Q = blkdiag(QL,QR)), so it must be
      % sliced by position -- NOT by the original-order logical hemisphere mask.
      mapL = elec2map_apply(data.smooth_cache.L, dvn(1:data.elec.nL));
      mapR = elec2map_apply(data.smooth_cache.R, dvn(data.elec.nL+1:end));

      % optional marginal matching
      switch lower(cfg.rankrescale)
        case 'no'
        case 'exact'
          [~, idxL] = sort(mapL); mapL(idxL) = sort(allmaps(1:nVertL,s));
          [~, idxR] = sort(mapR); mapR(idxR) = sort(allmaps(nVertL+1:end,s));
        case 'quantile'
          mapL = rescale_one_quantile(mapL, data.rescale_cache.L);
          mapR = rescale_one_quantile(mapR, data.rescale_cache.R);
        case 'zscore'
          mapL = rescale_one_z(mapL, data.rescale_cache.L);
          mapR = rescale_one_z(mapR, data.rescale_cache.R);
        otherwise
          error('Unknown rankrescale method: %s', cfg.rankrescale);
      end

      permMaps(:,s) = [mapL; mapR];
      if strcmp(cfg.keepsubsurrogates, 'yes')
          if strcmp(cfg.subsamplesurr, 'yes')
              if iter <= cfg.nsurrsamples
                  surr_subj(s,:,iter) = permMaps(:,s);
              end
          else
            surr_subj(s,:,iter) = permMaps(:,s);
          end
      end
    end
  else
    permMaps = allmaps; % empirical pass
  end

  % optional per-subject z-normalization before group t
  if strcmpi(cfg.normalize, 'yes')
    mu  = nanmean(permMaps, 1);
    sd  = nanstd(permMaps, 0, 1);
    permMaps = (permMaps - mu) ./ sd; % standardize
  end

  % vertex-wise one-sample t across subjects
  mu  = nanmean(permMaps, 2);
  sd  = sqrt(nanvar(permMaps, 0, 2));
  n   = double(sum(~isnan(permMaps), 2));
  tmap = mu ./ (sd ./ sqrt(n));
  tmap(isinf(tmap)) = nan;

  if strcmpi(cfg.equalize_n, 'yes')
    neff  = median(n(n>0));
    tmap  = tmap .* sqrt(neff ./ n);
  end
  tmap(mask) = nan;

  % update null or finalize observed stats
  [stat, posdistribution, negdistribution] = clusterstat(...
    tmap, sph, adj, postailcritval, negtailcritval, iter, cfg, posdistribution, negdistribution);

  if strcmp(cfg.keepgroupsurrogates, 'yes') && iter <= cfg.numrandomization
    surr_group(:,iter) = tmap;
  end
end

% ---- Outputs
stat.coverage         = coverage;
% The cache holds meshes and per-subject projectors/eigenbases; it is large and
% carries no information about the result, so it never travels inside stat.cfg.
cfgout = cfg;
cfgout.cache = [];
stat.cfg              = cfgout;
if strcmp(cfg.returncache, 'yes')
  cache = struct();
  cache.geom = struct('lh',lh,'rh',rh,'sph_l',sph_l,'sph_r',sph_r,'sph',sph,'adj',adj);
  subjc = struct('smooth_cache',{},'idxL',{},'idxR',{},'Q',{},'eigvals',{});
  for s = 1:nSubjects
    subjc(s).smooth_cache = varargin{s}.smooth_cache;
    subjc(s).idxL         = varargin{s}.elec.idxL;
    subjc(s).idxR         = varargin{s}.elec.idxR;
    subjc(s).Q            = varargin{s}.modes;
    subjc(s).eigvals      = varargin{s}.eigvals;
  end
  cache.subject = subjc;
  % Guard against silently reusing a cache built for a different montage or a
  % different smoothing geometry: both would invalidate the projectors.
  cache.signature = struct('kernel',cfg.kernel,'kernelwidth',cfg.kernelwidth, ...
    'graphsigma',cfg.graphsigma,'smooth',cfg.smooth,'nSubjects',nSubjects, ...
    'nchan',cellfun(@(d) numel(d.(cfg.parameter)), varargin(:)'));
  stat.cache = cache;
end
if strcmp(cfg.keepmaps, 'yes')
  stat.maps = allmaps;
end
if strcmp(cfg.keepsubsurrogates, 'yes')
  stat.surrogate_subject = surr_subj;
end
if strcmp(cfg.keepgroupsurrogates, 'yes')
  stat.surrogate_group   = surr_group;
end

% =========================================================================
% Subfunctions
% =========================================================================

function validate_sources(cfg, varargin)
% Fail fast, with actionable messages, before any expensive computation.
nSubjects = numel(varargin);
if nSubjects == 0
  error('SMOOTHstat:noData', ...
    'No source structures supplied. Use as: stat = SMOOTHstat(cfg, source{:})');
end

allx = [];
for s = 1:nSubjects
  d  = varargin{s};
  id = sprintf('subject %d of %d', s, nSubjects);

  if ~isstruct(d)
    error('SMOOTHstat:badInput', '%s: expected a struct, got %s.', id, class(d));
  end
  if ~isfield(d, cfg.parameter)
    error('SMOOTHstat:missingParameter', ...
      '%s: field ''%s'' not found (available: %s). Set cfg.parameter to the field holding your per-channel statistic.', ...
      id, cfg.parameter, strjoin(fieldnames(d)', ', '));
  end
  if ~isfield(d, 'elec') || ~isfield(d.elec, 'chanpos')
    error('SMOOTHstat:missingChanpos', ...
      '%s: source.elec.chanpos is required (electrode positions in fsaverage space).', id);
  end

  dv  = d.(cfg.parameter);
  pos = d.elec.chanpos;

  if ~isvector(dv)
    error('SMOOTHstat:badParameter', ...
      '%s: source.%s must be a vector, got size [%s].', id, cfg.parameter, num2str(size(dv)));
  end
  if numel(dv) ~= size(pos,1)
    error('SMOOTHstat:sizeMismatch', ...
      '%s: source.%s has %d values but elec.chanpos has %d rows. These must match.', ...
      id, cfg.parameter, numel(dv), size(pos,1));
  end
  if size(pos,2) ~= 3
    error('SMOOTHstat:badChanpos', ...
      '%s: elec.chanpos must be Nx3, got Nx%d.', id, size(pos,2));
  end
  if any(~isfinite(dv))
    error('SMOOTHstat:nonFiniteData', ...
      '%s: source.%s contains %d NaN/Inf value(s). Remove or interpolate these channels before calling SMOOTHstat -- they would silently corrupt the coverage map and group t-statistic.', ...
      id, cfg.parameter, sum(~isfinite(dv)));
  end
  if any(~isfinite(pos(:)))
    error('SMOOTHstat:nonFiniteChanpos', ...
      '%s: elec.chanpos contains non-finite coordinates.', id);
  end

  % hemisphere assignment is by sign of x; warn about ambiguous contacts
  nMid = sum(abs(pos(:,1)) < 1);
  if nMid > 0
    warning('SMOOTHstat:midlineChannels', ...
      '%s: %d channel(s) lie within 1 mm of the midline (x=0). Hemisphere is assigned by sign(x); x==0 is treated as right. Verify these contacts.', ...
      id, nMid);
  end
  allx = [allx; pos(:,1)]; %#ok<AGROW>
end

% coordinate-space sanity check: fsaverage surface RAS spans roughly +/-70 mm
rng_all = [min(allx) max(allx)];
if max(abs(rng_all)) > 200 || max(abs(rng_all)) < 20
  warning('SMOOTHstat:coordinateSpace', ...
    ['Electrode x-coordinates span [%.1f %.1f], which does not look like fsaverage surface RAS ' ...
     '(expected roughly [-70 70] mm). Hemisphere assignment and all mm-based kernels assume ' ...
     'fsaverage space -- results will be wrong if coordinates are in native, MNI, or voxel space.'], ...
    rng_all(1), rng_all(2));
end

function [map, hemi_cache, idx, e2m_cache] = elec2map_build(side, cfg, surf_pial, surf_sph, pos_all, dv_all, sigma_mm, hemi_mask, e2m_cache)
% Per-hemisphere projector build + single apply, caching heavy aux in e2m_cache.
% side      : 'L' or 'R'
% surf_*    : hemi meshes (pial & sphere.reg)
% pos_all   : Nch x 3 electrode coords (fsavg)
% dv_all    : Nch x 1 dependent variable
% sigma_mm  : Gaussian σ in mm if cfg.kernel='gaussian'; radius(mm) if boxcar
% hemi_mask : logical Nch x 1 (channels for THIS hemi)
% e2m_cache : struct with KD-trees/aux for reuse

nVert = size(surf_pial.pos,1);
map   = nan(nVert,1);
hemi_cache = struct('P',[],'mask',[],'nVert',nVert);
idx   = int32([]);

if ~any(hemi_mask)
  return
end

% subset electrodes & DV for this hemi
posH = pos_all(hemi_mask,:);
dvH  = dv_all(hemi_mask);

% nearest pial vertices (per hemi)
idx  = int32(knnsearch(surf_pial.pos, posH, 'K',1, 'NSMethod','kdtree'));

switch lower(cfg.smooth)
  % ====================== SPHERE / EUCLIDEAN ======================
  % Fast. Builds weights using Euclidean distances in fsaverage sphere.reg
  % coordinates. Distances are in sphere units.
  case 'sphere'
    % KD-tree on sphere verts
    if side=='L'
      if isempty(e2m_cache.sphere.kdtL), e2m_cache.sphere.kdtL = createns(surf_sph.pos,'NSMethod','kdtree'); end
      kdt = e2m_cache.sphere.kdtL;
    else
      if isempty(e2m_cache.sphere.kdtR), e2m_cache.sphere.kdtR = createns(surf_sph.pos,'NSMethod','kdtree'); end
      kdt = e2m_cache.sphere.kdtR;
    end

    I=[]; J=[]; D=[];
    if strcmpi(cfg.kernel,'gaussian')
      % σ in sphere units, support ~ 3σ
      sig_unit = mm2sphere(sigma_mm);
      R        = 3*max(sig_unit, eps);  
      cand_all = rangesearch(kdt, surf_sph.pos(double(idx),:), R);
      for e=1:numel(idx)
        s    = double(idx(e));
        cand = cand_all{e}(:);  if isempty(cand), continue, end
        dvec = sqrt(sum((surf_sph.pos(cand,:) - surf_sph.pos(s,:)).^2,2));
        keep = dvec <= R; cand=cand(keep); dvec=dvec(keep);
        I=[I; cand]; J=[J; repmat(e,numel(cand),1)]; D=[D; dvec];
      end
      Wvals = exp(-(D.^2)/(2*sig_unit^2));

    elseif strcmpi(cfg.kernel,'boxcar')
      % radius in sphere units, hard cutoff
      R = max(mm2sphere(cfg.kernelwidth), eps);
      cand_all = rangesearch(kdt, surf_sph.pos(double(idx),:), R);
      for e=1:numel(idx)
        s    = double(idx(e));
        cand = cand_all{e}(:);  if isempty(cand), continue, end
        dvec = sqrt(sum((surf_sph.pos(cand,:) - surf_sph.pos(s,:)).^2,2));
        keep = dvec <= R; cand=cand(keep);
        I=[I; cand]; J=[J; repmat(e,numel(cand),1)];
      end
      Wvals = ones(size(I));
    else
      error('cfg.kernel must be ''gaussian'' or ''boxcar''.');
    end

    % % COL
    W   = sparse(I,J,Wvals,nVert,numel(idx));
    col = full(sum(W,1));
    C   = spdiags(1./max(col,eps)', 0, size(W,2), size(W,2));
    P   = W * C; 

    % % ROW
    den = full(sum(W,2));
    m   = den>0;
    % P   = spdiags(1./max(den,eps),0,nVert,nVert) * W;

    % ====================== PIAL / GEODESIC (mm) ======================
    % Slower but more anatomically faithful. Uses geodesic distances on
    % the pial mesh (mm) via a Dijkstra traversal truncated at ~3·σ.
    % Recommended when precise folding geometry matters or kernels are
    % large. Expect ~5–10× slower projector builds on typical meshes.
  case 'pial'
    % ensure geodesic aux for this hemi is in cache
    if side=='L'
      if isempty(e2m_cache.pial.L), e2m_cache.pial.L = local_build_hemi_aux(surf_pial.pos, surf_pial.tri); end
      aux = e2m_cache.pial.L;
    else
      if isempty(e2m_cache.pial.R), e2m_cache.pial.R = local_build_hemi_aux(surf_pial.pos, surf_pial.tri); end
      aux = e2m_cache.pial.R;
    end

    I=[]; J=[]; D=[];
    if strcmpi(cfg.kernel,'gaussian')
      sig = max(sigma_mm, eps);  % σ in mm
      R   = 3*sig;
      for e=1:numel(idx)
        s = double(idx(e));
        % truncated Dijkstra to R (mm)
        d  = inf(nVert,1); d(s)=0; seen=false(nVert,1);
        heap_idx=int32(0); heap_key=0; hsz=1; heap_idx(1)=int32(s); heap_key(1)=0;
        while hsz>0
          v  = heap_idx(1); dvh = heap_key(1);
          hsz = hsz-1;
          if hsz>0
            heap_idx(1)=heap_idx(hsz+1); heap_key(1)=heap_key(hsz+1);
            p=1; while true
              l=2*p; r=l+1; m=p;
              if l<=hsz && heap_key(l)<heap_key(m), m=l; end
              if r<=hsz && heap_key(r)<heap_key(m), m=r; end
              if m==p, break; end
              [heap_idx(p),heap_idx(m)] = deal(heap_idx(m),heap_idx(p));
              [heap_key(p),heap_key(m)] = deal(heap_key(m),heap_key(p));
              p=m;
            end
          end
          if seen(v) || dvh>R, continue; end
          seen(v)=true;
          nb = aux.nbr{v}; wg = aux.wgt{v};
          for k=1:numel(nb)
            u = nb(k);
            alt = dvh + wg(k);
            if alt < d(u) && alt <= R
              d(u)=alt; hsz=hsz+1;
              if hsz>numel(heap_idx)
                newsz=max(hsz, ceil(numel(heap_idx)*1.5));
                heap_idx(end+1:newsz)=0; heap_idx=int32(heap_idx);
                heap_key(end+1:newsz)=0;
              end
              heap_idx(hsz)=u; heap_key(hsz)=alt;
              q=hsz; while q>1
                p=floor(q/2); if heap_key(p)<=heap_key(q), break; end
                [heap_idx(p),heap_idx(q)] = deal(heap_idx(q),heap_idx(p));
                [heap_key(p),heap_key(q)] = deal(heap_key(q),heap_key(p)); q=p;
              end
            end
          end
        end
        hit = find(isfinite(d) & d<=R);
        if ~isempty(hit)
          I=[I; hit]; J=[J; repmat(e,numel(hit),1)]; D=[D; d(hit)];
        end
      end
      Wvals = exp(-(D.^2)/(2*(sig^2)));

    elseif strcmpi(cfg.kernel,'boxcar')
      R = max(cfg.kernelwidth, eps);  % radius in mm; hard cutoff
      for e=1:numel(idx)
        s = double(idx(e));
        d  = inf(nVert,1); d(s)=0; seen=false(nVert,1);
        heap_idx=int32(0); heap_key=0; hsz=1; heap_idx(1)=int32(s); heap_key(1)=0;
        while hsz>0
          v  = heap_idx(1); dvh = heap_key(1);
          hsz = hsz-1;
          if hsz>0
            heap_idx(1)=heap_idx(hsz+1); heap_key(1)=heap_key(hsz+1);
            p=1; while true
              l=2*p; r=l+1; m=p;
              if l<=hsz && heap_key(l)<heap_key(m), m=l; end
              if r<=hsz && heap_key(r)<heap_key(m), m=r; end
              if m==p, break; end
              [heap_idx(p),heap_idx(m)] = deal(heap_idx(m),heap_idx(p));
              [heap_key(p),heap_key(m)] = deal(heap_key(m),heap_key(p));
              p=m;
            end
          end
          if seen(v) || dvh>R, continue; end
          seen(v)=true;
          nb = aux.nbr{v}; wg = aux.wgt{v};
          for k=1:numel(nb)
            u = nb(k); alt = dvh + wg(k);
            if alt < d(u) && alt <= R
              d(u)=alt; hsz=hsz+1;
              if hsz>numel(heap_idx)
                newsz=max(hsz, ceil(numel(heap_idx)*1.5));
                heap_idx(end+1:newsz)=0; heap_idx=int32(heap_idx);
                heap_key(end+1:newsz)=0;
              end
              heap_idx(hsz)=u; heap_key(hsz)=alt;
              q=hsz; while q>1
                p=floor(q/2); if heap_key(p)<=heap_key(q), break; end
                [heap_idx(p),heap_idx(q)]=deal(heap_idx(q),heap_idx(p));
                [heap_key(p),heap_key(q)]=deal(heap_key(q),heap_key(p)); q=p;
              end
            end
          end
        end
        hit = find(isfinite(d) & d<=R);
        if ~isempty(hit)
          I=[I; hit]; J=[J; repmat(e,numel(hit),1)];
        end
      end
      Wvals = ones(size(I));
    else
      error('cfg.kernel must be ''gaussian'' or ''boxcar''.');
    end

    W   = sparse(I,J,Wvals,nVert,numel(idx));
    den = full(sum(W,2));
    m   = den>0;
    P   = spdiags(1./max(den,eps),0,nVert,nVert) * W;

  otherwise
    error('Unknown smoothing mode "%s"', cfg.smooth);
end

% apply once, return cache for reuse
tmp = P * dvH; map(m) = tmp(m);
hemi_cache.P    = P;
hemi_cache.mask = m;

% -------------------------------------------------------------------------
function map = elec2map_apply(hemi_cache, dv_hemi)
% Reuse projector/mask for permutations.
map = nan(hemi_cache.nVert,1);
if ~isempty(hemi_cache.P)
  tmp = hemi_cache.P * dv_hemi;
  if ~isempty(hemi_cache.mask)
    map(hemi_cache.mask) = tmp(hemi_cache.mask);
  else
    map = tmp;
  end
end

% -------------------------------------------------------------------------
function aux = local_build_hemi_aux(V, F)
nV = size(V,1);
E  = [F(:,[1 2]); F(:,[2 3]); F(:,[3 1]); F(:,[2 1]); F(:,[3 2]); F(:,[1 3])];
dV = V(E(:,1),:) - V(E(:,2),:);
W  = sqrt(sum(dV.^2,2));
A  = sparse(E(:,1),E(:,2),W,nV,nV); A = min(A, A.');
[i,j,w] = find(A);
nbr = accumarray(i, j, [nV 1], @(x){int32(x)});
wgt = accumarray(i, w, [nV 1], @(x){x});
aux = struct('nbr',{nbr}, 'wgt',{wgt});

% -------------------------------------------------------------------------
function su = mm2sphere(mm)
% Convert a millimetre scale to fsaverage sphere.reg units.
%
% The factor is the empirical pial -> sphere.reg edge-length ratio on fsaverage
% (measured median 1.24 units/mm). Spherical registration preserves topology but
% not metric, so the ratio varies across cortex (5th-95th pct 0.77-2.34); 1.25 is
% the average correspondence, not a uniform one.
su = 1.25 * mm;

% -------------------------------------------------------------------------
function [Q, lam] = graph_eigenbasis(pos, sigma_graph)
% (a) Build Gaussian-weighted graph, normalized Laplacian
% (b) Return eigenvectors (columns) ordered by eigenvalue + eigenvalues
%
% sigma_graph (mm) controls the width of the electrode-level Gaussian
% weighting on the subject’s *electrode graph* (used to form a normalized
% Laplacian and its eigenbasis for spectral surrogates). We first convert
% this mm scale to sphere units per hemi (same local scale factor as
% above), then build W_ij = exp(-d_ij^2 / (2*σ_graph^2)) among electrodes.
% Larger sigma_graph → more low-frequency modes (smoother surrogates).
n = size(pos,1);
if n==0
  Q = zeros(0,0); lam = []; return
elseif n==1
  Q = 1; lam = 0; return
end
r = sqrt(sum(pos.^2, 2));
U = pos ./ r;  % unit vectors
C = U * U.';                 
C = max(-1, min(1, C));
theta = acos(C);             % radians
D = 100 * theta;
W  = exp(-(D.^2) ./ (2*(sigma_graph^2)));
W(1:n+1:end) = 0; % prevent self-loops
d  = sum(W,2); d(d==0) = 1;
Dm = spdiags(1./sqrt(d), 0, n, n);
L  = speye(n) - (Dm * W * Dm);
L  = (L + L')/2;
[Qtmp, Lam] = eig(L);
lam = diag(Lam);
[lam, ord] = sort(lam, 'ascend'); % smooth → rough
Q = Qtmp(:,ord);

% -------------------------------------------------------------------------
function a2 = spectral_bandrotate(a, lam)
% Random orthogonal rotation within eigenvalue bands (preserves band power)
% Inputs:
%   a   : modal coefficients (length m)
%   lam : eigenvalues aligned to a (length m)
% Output:
%   a2  : rotated coefficients (same length)
if isempty(a) || isempty(lam), a2=a; return, end
B=20; edges=quantile(lam, linspace(0,1,B+1)); a2=a;
for b=1:B
  if b<B, idx = find(lam>=edges(b) & lam<edges(b+1));
  else,    idx = find(lam>=edges(b) & lam<=edges(b+1));
  end
  k=numel(idx);
  if k>1
    G=randn(k,k); [Q,~]=qr(G,0); if det(Q)<0, Q(:,1)=-Q(:,1); end
    a2(idx)=Q*a(idx);
  end
end

% -------------------------------------------------------------------------
function out = rescale_one_quantile(in, c)
out=in; m=c.mask;
if any(m)
  vals=in(m);
  [~,ord]=sort(vals,'ascend');
  tmp=vals;
  tmp(ord)=c.tgt_sorted;
  out(m)=tmp;
end

% -------------------------------------------------------------------------
function out = rescale_one_z(in, c)
out=in; m=c.mask;
if any(m)
  x=in(m);
  x=(x-mean(x))./max(std(x),eps);
  out(m)=c.mu + c.sd * x;
end

% -------------------------------------------------------------------------
function [stat, posdistribution, negdistribution] = clusterstat(tmap, sph, adj, postailcritval, negtailcritval, iter, cfg, posdistribution, negdistribution)
if iter <= cfg.numrandomization
  % accumulate null distributions
  posclusrnd = clusterfunc(tmap >= postailcritval, sph, adj);
  negclusrnd = clusterfunc(tmap <= negtailcritval, sph, adj);
  if any(posclusrnd(:))
    poscstats = accumarray(posclusrnd(posclusrnd>0), tmap(posclusrnd>0));
    posdistribution(iter) = max(poscstats);
  end
  if any(negclusrnd(:))
    negcstats = accumarray(negclusrnd(negclusrnd>0), tmap(negclusrnd>0));
    negdistribution(iter) = min(negcstats);
  end
  stat = [];
else
  % finalize observed
  prb_pos = ones(size(tmap)); prb_neg = ones(size(tmap));

  % positive clusters
  posclusobs = clusterfunc(tmap >= postailcritval, sph, adj);
  if any(posclusobs(:))
    poscstats = accumarray(posclusobs(posclusobs>0), tmap(posclusobs>0));
    [sorted_cstats, sort_idx] = sort(poscstats, 'descend');
    map = zeros(max(posclusobs(:)),1); map(sort_idx) = 1:numel(sorted_cstats);
    posremap = nan(size(tmap)); valid = posclusobs > 0;
    posremap(valid) = map(posclusobs(valid)); posclusobs = posremap;
    for j = 1:numel(sorted_cstats)
      posclusters(j).clusterstat = sorted_cstats(j);
      posclusters(j).prob = (sum(posdistribution >= sorted_cstats(j)) + 1) / (numel(posdistribution) + 1);
    end
    prb_pos(valid) = [posclusters(posclusobs(valid)).prob];
  else
    posclusters = [];
  end

  % negative clusters
  negclusobs = clusterfunc(tmap <= negtailcritval, sph, adj);
  if any(negclusobs(:))
    negcstats = accumarray(negclusobs(negclusobs>0), tmap(negclusobs>0));
    [sorted_cstats, sort_idx] = sort(negcstats, 'ascend');
    map = zeros(max(negclusobs(:)),1); map(sort_idx) = 1:numel(sorted_cstats);
    negremap = nan(size(tmap)); valid = negclusobs > 0;
    negremap(valid) = map(negclusobs(valid)); negclusobs = negremap;
    for j = 1:numel(sorted_cstats)
      negclusters(j).clusterstat = sorted_cstats(j);
      negclusters(j).prob = (sum(negdistribution <= sorted_cstats(j)) + 1) / (numel(negdistribution) + 1);
    end
    prb_neg(valid) = [negclusters(negclusobs(valid)).prob];
  else
    negclusters = [];
  end

  % ---- Assemble output. Nulls are accumulated separately per tail.
  %
  % stat.prob is the ONE-TAILED cluster p-value at each vertex (the smaller of the
  % two tails where both apply). It is deliberately NOT corrected across tails --
  % use stat.mask, which applies the alpha appropriate to cfg.tail:
  %
  %   'positive'/'negative' : one-tailed test at cfg.alpha (default 0.05)
  %   'both'                : two-tailed, cfg.alpha/2 per tail (Bonferroni)
  %
  % Thresholding stat.prob directly at 0.05 while considering BOTH tails would
  % double the family-wise error rate; that is what cfg.tail guards against.
  stat = struct();
  stat.stat                = tmap;
  stat.prob                = min(prb_pos, prb_neg);
  stat.prob(stat.prob > 1) = 1;

  switch lower(cfg.tail)
    case 'positive'
      stat.mask = prb_pos <= cfg.alpha;
    case 'negative'
      stat.mask = prb_neg <= cfg.alpha;
    case 'both'
      stat.mask = stat.prob <= cfg.alpha/2;
    otherwise
      error('SMOOTHstat:badTail', ...
        'cfg.tail must be ''positive'', ''negative'', or ''both'' (got ''%s'').', cfg.tail);
  end

  stat.posclusters         = posclusters;
  stat.posclusterslabelmat = posclusobs;
  stat.posdistribution     = posdistribution;
  stat.negclusters         = negclusters;
  stat.negclusterslabelmat = negclusobs;
  stat.negdistribution     = negdistribution;
end

% -------------------------------------------------------------------------
function clustermap = clusterfunc(mask, ctx, adjacency)
sig_idx = find(mask);
if isempty(sig_idx)
  clustermap = nan(size(ctx.pos,1),1);
  return;
end
subgraph = adjacency(sig_idx, sig_idx);
G = graph(subgraph);
bins = conncomp(G);
clustermap = nan(size(ctx.pos,1),1);
clustermap(sig_idx) = bins;
