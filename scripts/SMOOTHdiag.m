function diag = SMOOTHdiag(stat)

% SMOOTHdiag - diagnostic operating on SMOOTHstat output struct. Primarily
% for internal use and data visualization in the paper
%
% Usage:
%   diag = SMOOTHdiag(stat)
%
% Requirements:
%   stat.cfg.fshome must point to FreeSurfer home (for fsaverage sphere.reg).
%
% Notes:
%   - Works with or without surrogates; skips null-based diagnostics if absent.
%   - Uses fast (subsampled) diagnostics; does NOT affect inferential stats.

% default settings
POP_MODE              = 'maxAbs';     % 'maxAbs' | 'minDouble'
RANDOMIZE_TIES        = true;
alpha                 = 0.05;

SPATIAL_VERBOSE       = false;

SUBJ_SPATIAL_VSAMP    = Inf;          % max vertices per block
SUBJ_SPATIAL_NPERM    = 100;          % surrogates/subject for I & C
SUBJ_SPATIAL_MINV     = 400;          % min vertices per block

VALDIST_FAST          = true;         % binned CDF acceleration
VALDIST_NBINS         = 256;
VALDIST_VSAMP         = 5000;         % vertices/subject
VALDIST_NPERM         = 100;          % surrogates/subject

NULLPAIR_SAMPLES      = 100;          % pairs for null self-corr

fprintf('\n===================== SMOOTHdiag =====================\n')

% basic checks
cfg   = stat.cfg;
Vtot  = numel(stat.stat);

if ~isfield(cfg,'fshome') || ~ischar(cfg.fshome)
  error('SMOOTHdiag:MissingFShome', 'stat.cfg.fshome is required to rebuild fsaverage adjacency.');
end

% Read fsaverage spheres & adjacency
sph_l   = ft_read_headshape(fullfile(cfg.fshome,'subjects','fsaverage','surf','lh.sphere.reg'));
sph_r   = ft_read_headshape(fullfile(cfg.fshome,'subjects','fsaverage','surf','rh.sphere.reg'));
nVertL  = size(sph_l.pos,1);
nVertR  = size(sph_r.pos,1);
if (nVertL + nVertR) ~= Vtot
  error('SMOOTHdiag:VertexMismatch', 'stat.stat length (%d) != nVertL+nVertR (%d).', Vtot, nVertL+nVertR);
end
tri = [sph_l.tri; sph_r.tri + nVertL];
i = [tri(:,1); tri(:,1); tri(:,2); tri(:,2); tri(:,3); tri(:,3)];
j = [tri(:,2); tri(:,3); tri(:,1); tri(:,3); tri(:,1); tri(:,2)];
adj = sparse(i,j,1,Vtot,Vtot); adj = adj | adj';  % undirected

% Hemis masks
mLH = false(Vtot,1); mLH(1:nVertL)     = true;
mRH = false(Vtot,1); mRH(nVertL+1:end) = true;
mBT = true(Vtot,1);

% Availability flags
have_maps       = isfield(stat,'maps') && ~isempty(stat.maps);
have_surr_group = isfield(stat,'surrogate_group')   && ~isempty(stat.surrogate_group);
have_surr_subj  = isfield(stat,'surrogate_subject') && ~isempty(stat.surrogate_subject);

if ~have_maps
  warning('SMOOTHdiag:NoSubjectMaps','stat.maps not found; subject-level diagnostics will be limited.');
end

% Pull arrays (if available)
if have_maps
  allmaps = stat.maps;               % V x S
  [~,S]   = size(allmaps);
else
  allmaps = [];
  S       = 0; %#ok<NASGU>
end

if have_surr_group
  surr_group = stat.surrogate_group; % V x P
  P2         = size(surr_group,2);
else
  surr_group = [];
  P2         = 0;
end

% Canonicalize subject surrogates to V×S×P
if have_surr_subj
  sz = size(stat.surrogate_subject);
  if numel(sz)==3 && sz(1) == Vtot
    surr_subj = stat.surrogate_subject;          % V×S×P
  else
    surr_subj = permute(stat.surrogate_subject,[2 1 3]); % S×V×P → V×S×P
  end
else
  surr_subj = [];
end

t_emp    = stat.stat(:);
coverage = stat.coverage(:);

% cal
if isfield(stat,'posdistribution') && isfield(stat,'negdistribution') ...
    && ~isempty(stat.posdistribution) && ~isempty(stat.negdistribution)

  P   = numel(stat.posdistribution);
  posd = stat.posdistribution; posd(~isfinite(posd)) = 0;
  negd = stat.negdistribution; negd(~isfinite(negd)) = 0;

  switch lower(POP_MODE)
    case 'maxabs'
      T = max(posd, -negd);
      [pp, KS_abs, KS_signed] = fast_pop_pvals_ranked(T, RANDOMIZE_TIES);
      modeStr = 'maxAbs (two-sided T = max(pos, -neg))';
    otherwise
      [pp_pos, ~, ~] = fast_pop_pvals_ranked(posd, RANDOMIZE_TIES, 'high');
      [pp_neg, ~, ~] = fast_pop_pvals_ranked(negd, RANDOMIZE_TIES, 'low');
      pp = 2*min(pp_pos, pp_neg); pp(pp>1)=1;
      [~, KS_abs, KS_signed] = fast_pop_pvals_ranked(1-pp, false);
      modeStr = 'minDouble (min one-sided p × 2)';
  end

  FWERat = mean(pp<=alpha);
  D10 = 1.22/sqrt(P); D05 = 1.36/sqrt(P); D01 = 1.63/sqrt(P);
  diag.POP = struct('P',P,'FWER_at_005',FWERat,'KS_abs',KS_abs,'KS_signed',KS_signed, ...
    'mode',modeStr,'randomized_ties',RANDOMIZE_TIES,'D10',D10,'D05',D05,'D01',D01);

  % fprintf('\n[QC] POP (cluster null)\n')
  line_metric('FWER@0.05', FWERat, 'target≈0.05', badge_around(FWERat,0.05,0.02,0.05))
  dir = 'balanced'; if KS_signed>1e-8, dir='liberal (too many small p)'; elseif KS_signed<-1e-8, dir='conservative (too many large p)'; end
  line_metric('KS to Uniform', KS_abs, sprintf('signed %.3f, %s', KS_signed, dir), badge_low(KS_abs,D10,D05))
  line_note(sprintf('Mode: %s; Tie randomization: %s; P=%d; KS refs: D10≈%.3f D05≈%.3f D01≈%.3f', ...
    modeStr, yn(RANDOMIZE_TIES), P, D10, D05, D01))
else
  % fprintf('\n[QC] POP (cluster null): skipped (no stored null distributions)\n')
end

% subj spatial smoothness
if have_maps && have_surr_subj
  % fprintf('\n[QC] Subject spatial smoothness (graph): Moran''s I ↑ smoother, Geary''s C ↓ smoother\n')
  % fprintf('     FAST mode: shared mask per block, subsample V≤%d, %d perms/subject\n', SUBJ_SPATIAL_VSAMP, SUBJ_SPATIAL_NPERM)
  blkLH = subject_spatial_block_fast(allmaps, surr_subj, adj, mLH, 'LH',   SUBJ_SPATIAL_VSAMP, SUBJ_SPATIAL_NPERM, SUBJ_SPATIAL_MINV);
  blkRH = subject_spatial_block_fast(allmaps, surr_subj, adj, mRH, 'RH',   SUBJ_SPATIAL_VSAMP, SUBJ_SPATIAL_NPERM, SUBJ_SPATIAL_MINV);
  blkBT = subject_spatial_block_fast(allmaps, surr_subj, adj, mBT, 'BOTH', SUBJ_SPATIAL_VSAMP, SUBJ_SPATIAL_NPERM, SUBJ_SPATIAL_MINV);
  print_spatial_row('LH', blkLH, SPATIAL_VERBOSE)
  print_spatial_row('RH', blkRH, SPATIAL_VERBOSE)
  print_spatial_row('BOTH', blkBT, SPATIAL_VERBOSE)
  diag.subject.spatial = struct('BOTH',blkBT,'LH',blkLH,'RH',blkRH);
else
  % fprintf('\n[QC] Subject spatial smoothness: skipped (need stat.maps AND stat.surrogate_subject)\n')
end

% subject empirical vs surrogate correlation 
if have_maps && have_surr_subj
  % fprintf('\n[QC] Subject empirical vs subject-surrogate correlation (mean across perms)\n')
  SUBJ_CORR_MAXPERM = 200;   % cap for speed

  % Use the same shared-vertex rule as spatial smoothness (robust w.r.t. NaNs)
  sLH = fast_corr_emp_surr_subject(allmaps, surr_subj, adj, mLH, SUBJ_SPATIAL_MINV, SUBJ_CORR_MAXPERM);
  sRH = fast_corr_emp_surr_subject(allmaps, surr_subj, adj, mRH, SUBJ_SPATIAL_MINV, SUBJ_CORR_MAXPERM);
  sBT = fast_corr_emp_surr_subject(allmaps, surr_subj, adj, mBT, SUBJ_SPATIAL_MINV, SUBJ_CORR_MAXPERM);

  print_corr_row('LH',   sLH);
  print_corr_row('RH',   sRH);
  print_corr_row('BOTH', sBT);

  diag.subject.corr = struct('LH',sLH,'RH',sRH,'BOTH',sBT);
else
  % fprintf('\n[QC] Subject empirical vs subject-surrogate correlation: skipped (need stat.maps AND stat.surrogate_subject)\n')
end

% subject value distributions
if have_maps && have_surr_subj
  % fprintf('\n[QC] Subject value distributions (emp vs pooled surrogates)\n')
  vdALL = subject_value_dist_fast(allmaps, surr_subj, mBT, 'ALL', VALDIST_VSAMP, VALDIST_NPERM, VALDIST_NBINS, VALDIST_FAST);
  print_value_dist_row('ALL', vdALL)
  diag.subject.valuedist = struct('ALL',vdALL);
else
  % fprintf('\n[QC] Subject value distributions: skipped (need stat.maps AND stat.surrogate_subject)\n')
end

% group correlation (empirical vs surrogate) 
if have_surr_group
  % fprintf('\n[QC] Group correlation (empirical t vs surrogate t)\n')

  % Per-hemisphere + both
  [gLH.mu,gLH.sd,gLH.n,gLH.t,gLH.p] = fast_corr_emp_surr_group(t_emp, surr_group, mLH);
  [gRH.mu,gRH.sd,gRH.n,gRH.t,gRH.p] = fast_corr_emp_surr_group(t_emp, surr_group, mRH);
  [gBT.mu,gBT.sd,gBT.n,gBT.t,gBT.p] = fast_corr_emp_surr_group(t_emp, surr_group, mBT);

  print_corr_row('LH',   gLH);
  print_corr_row('RH',   gRH);
  print_corr_row('BOTH', gBT);

  diag.group.corr = struct('LH',gLH,'RH',gRH,'BOTH',gBT);
else
  % fprintf('\n[QC] Group correlation: skipped (no stat.surrogate_group)\n')
end

% coverage correlation 
cvg = coverage;
v   = isfinite(t_emp) & isfinite(cvg);
if have_surr_group
  [mu_r, sd_r, r_emp, p_two] = fast_cov_corr(t_emp, surr_group, cvg, v);
  % fprintf('\n[QC] Coverage correlation (t vs coverage)  ')
  % fprintf('emp r=%.4f | null μ=%.4f sd=%.4f | p(two)=%.3g  %s\n', r_emp, mu_r, sd_r, p_two, badge_hiP(p_two,0.20,0.10))
  diag.group.coverage = struct('empirical_r',r_emp,'null_r_mean',mu_r,'null_r_sd',sd_r,'p_two',p_two);
else
  ok = v;
  r_emp = corr(t_emp(ok), double(cvg(ok)));
  % fprintf('\n[QC] Coverage correlation (t vs coverage)  emp r=%.4f | (null unavailable)\n', r_emp)
  diag.group.coverage = struct('empirical_r',r_emp,'null_r_mean',NaN,'null_r_sd',NaN,'p_two',NaN);
end

% group spatial smoothness (t-map) on graph
% fprintf('\n[QC] Group spatial smoothness (t-map) on graph\n')
[GI, GC] = group_spatial_smoothness_block(t_emp, have_surr_group, surr_group, adj);
line_metric('Moran''s I', GI.empirical, ...
  sprintf('null μ=%.4f sd=%.4f | p=%.3g', GI.null_mu, GI.null_sd, GI.p_two), ...
  badge_hiP(GI.p_two,0.20,0.10))
line_metric('Geary''s C', GC.empirical, ...
  sprintf('null μ=%.4f sd=%.4f | p=%.3g', GC.null_mu, GC.null_sd, GC.p_two), ...
  badge_hiP(GC.p_two,0.20,0.10))
diag.group.spatial = struct('moranI',GI,'gearyC',GC);

% null self-correlation
if have_surr_group && P2 >= 2
  M  = min(NULLPAIR_SAMPLES, max(10, floor(P2/2)));
  a = randi(P2, M, 1);
  b = randi(P2, M, 1);
  same = (b==a); b(same) = mod(b(same), P2) + 1;
  vals = nan(M,1);
  for m = 1:M
    xm = surr_group(:,a(m));
    ym = surr_group(:,b(m));
    ok = isfinite(xm) & isfinite(ym);
    if nnz(ok) >= 3
      xx = xm(ok) - mean(xm(ok));
      yy = ym(ok) - mean(ym(ok));
      den = sqrt(sum(xx.^2) * sum(yy.^2));
      if den > 0
        vals(m) = (xx' * yy) / den;
      end
    end
  end
  mu = nanmean(vals); sd = nanstd(vals); nn = sum(isfinite(vals));
  % fprintf('\n[QC] Null self-corr (two surrogate t-maps)  μ=%.4f ± %.4f (n=%d)  %s\n', ...
  %   mu, sd, nn, badge_around(mu,0,0.02,0.05))
  diag.group.null_selfcorr = struct('mu',mu,'sd',sd,'n',nn);
else
  % fprintf('\n[QC] Null self-corr: skipped (need ≥2 surrogate t-maps)\n')
  diag.group.null_selfcorr = struct('mu',NaN,'sd',NaN,'n',0);
end

% fprintf('======================================================\n')
end % ===== end main function =====


% subfunctions

function [GI, GC] = group_spatial_smoothness_block(t_emp, have_null, surr_group, adj)
mask = isfinite(t_emp);
if ~any(mask)
  GI = struct('empirical',NaN,'null_mu',NaN,'null_sd',NaN,'p_two',NaN);
  GC = GI; return
end
W = adj(mask,mask);
if nnz(W)==0
  GI = struct('empirical',NaN,'null_mu',NaN,'null_sd',NaN,'p_two',NaN);
  GC = GI; return
end
x = t_emp(mask);
I_emp = moran_I(x, W);
C_emp = geary_C(x, W);

I_null = []; C_null = [];
if have_null
  Ptot = size(surr_group,2);
  useP = min(Ptot, max(50, min(400, Ptot)));
  X = surr_group(mask, 1:useP);
  I_null = fast_moran_I_batch(X, W);
  C_null = fast_geary_C_batch(X, W);
  I_null = I_null(isfinite(I_null));
  C_null = C_null(isfinite(C_null));
end

GI = summarize_spatial_null(I_emp, I_null);
GC = summarize_spatial_null_C(C_emp, C_null);
end

function GI = summarize_spatial_null(I_emp, I_null)
if isempty(I_null)
  GI = struct('empirical',I_emp,'null_mu',NaN,'null_sd',NaN,'p_two',NaN);
else
  pI = 2*min(mean(I_null>=I_emp), mean(I_null<=I_emp)); pI=min(max(pI,0),1);
  GI = struct('empirical',I_emp,'null_mu',mean(I_null), 'null_sd',std(I_null), 'p_two',pI);
end
end

function GC = summarize_spatial_null_C(C_emp, C_null)
if isempty(C_null)
  GC = struct('empirical',C_emp,'null_mu',NaN,'null_sd',NaN,'p_two',NaN);
else
  pC = 2*min(mean(C_null<=C_emp), mean(C_null>=C_emp)); pC=min(max(pC,0),1);
  GC = struct('empirical',C_emp,'null_mu',mean(C_null), 'null_sd',std(C_null), 'p_two',pC);
end
end

function I = fast_moran_I_batch(X, W)
[~,P] = size(X); I = nan(1,P);
for j = 1:P
  m = isfinite(X(:,j));
  if nnz(m) < 3, continue, end
  Wc = W(m,m); if nnz(Wc)==0, continue, end
  x = X(m,j); x = x - mean(x);
  den = sum(x.^2); if den <= 0, continue, end
  S0 = full(sum(Wc(:))); if S0 == 0, continue, end
  I(j) = (numel(x)/S0) * (x' * (Wc * x)) / den;
end
end

function C = fast_geary_C_batch(X, W)
[~,P] = size(X); C = nan(1,P);
for j = 1:P
  m = isfinite(X(:,j));
  if nnz(m) < 3, continue, end
  Wc = W(m,m); if nnz(Wc)==0, continue, end
  x = X(m,j); x = x - mean(x);
  den = sum(x.^2); if den <= 0, continue, end
  d  = full(sum(Wc,2));
  L  = spdiags(d,0,nnz(m),nnz(m)) - Wc;
  S0 = full(sum(Wc(:))); if S0 == 0, continue, end
  num = 2 * (x' * (L * x));
  C(j) = ((numel(x)-1)/(2*S0)) * (num / den);
end
end

function block = subject_spatial_block_fast(allmaps_, surr_VSP, adj_, hemiMask_, label_, ~, KPERM, MINV)
% SUBJECT_SPATIAL_BLOCK_FAST (all-vertices version)
% - Uses ALL covered vertices per subject (no subsampling).
% - Deterministic surrogate selection: first K columns.
% - Chunked null evaluation to control memory.
%
% Inputs (kept for compatibility):
%   allmaps_  : V×S empirical maps
%   surr_VSP  : V×S×P surrogate maps
%   adj_      : V×V sparse adjacency (same vertex set)
%   hemiMask_ : V×1 logical (LH/RH/BOTH)
%   label_    : string label
%   ~         : (VMAX) ignored in this version
%   KPERM     : #surrogates per subject to use (cap)
%   MINV      : minimum vertices in mask

V  = size(allmaps_,1);
S  = size(allmaps_,2);
P_ = size(surr_VSP,3);

I_emp = nan(S,1); C_emp = nan(S,1);
I_mu  = nan(S,1); I_sd  = nan(S,1); I_p = nan(S,1);
C_mu  = nan(S,1); C_sd  = nan(S,1); C_p = nan(S,1);

K = min(KPERM, P_);            % deterministic first-K
CHUNK = 200;                    % chunk size for null eval (tune if needed)

for s_ = 1:S
  e = allmaps_(:,s_);
  % vertices where empirical is finite and at least 1 surrogate is finite
  fin_emp = isfinite(e) & hemiMask_;

  % require MINV overlap with surrogates (on average) to be robust
  % count finite surrogate entries across the first-K perms
  fin_surr_count = sum(isfinite(surr_VSP(:,s_,1:K)), 3);
  shared = fin_emp & (fin_surr_count >= max(1, ceil(0.8*1))); % >=1 is OK; adjust if you want ≥10, etc.

  vidx = find(shared);
  if numel(vidx) < MINV, continue, end

  x_emp = e(vidx);
  W = adj_(vidx, vidx);
  if nnz(W)==0 || nanstd(x_emp)<=eps, continue, end

  % Empirical Moran's I / Geary's C
  I_emp(s_) = moran_I(x_emp, W);
  C_emp(s_) = geary_C(x_emp, W);

  % Null (chunked over first-K permutations)
  Inull = nan(1,K); Cnull = nan(1,K);
  k0 = 1;
  while k0 <= K
    k1 = min(k0+CHUNK-1, K);
    X = surr_VSP(vidx, s_, k0:k1);            % (#verts) × chunk
    % drop columns with ~constant or too many NaNs
    bad = false(1, size(X,3));
    for j = 1:size(X,3)
      y = X(:,:,j);
      if ~any(isfinite(y)) || nanstd(y)<=eps, bad(j) = true; end
    end
    good = ~bad;
    if any(good)
      Xg = reshape(X(:,:,good), numel(vidx), sum(good));
      Inull(k0-1 + find(good)) = fast_moran_I_batch(Xg, W);
      Cnull(k0-1 + find(good)) = fast_geary_C_batch(Xg, W);
    end
    k0 = k1 + 1;
  end

  Inull = Inull(isfinite(Inull)); Cnull = Cnull(isfinite(Cnull));
  if ~isempty(Inull)
    I_mu(s_) = mean(Inull); I_sd(s_) = std(Inull);
    I_p(s_)  = 2*min(mean(Inull >= I_emp(s_)), mean(Inull <= I_emp(s_)));
  end
  if ~isempty(Cnull)
    C_mu(s_) = mean(Cnull); C_sd(s_) = std(Cnull);
    C_p(s_)  = 2*min(mean(Cnull <= C_emp(s_)), mean(Cnull >= C_emp(s_)));
  end
end

block.label  = label_;
block.moranI = struct( ...
  'empirical_mu', nanmean(I_emp), 'empirical_sd', nanstd(I_emp), ...
  'null_mu_mu',   nanmean(I_mu),  'null_mu_sd',   nanstd(I_mu), ...
  'p_sig_frac',   mean(I_p<=0.05,'omitnan'), 'n', sum(isfinite(I_emp)), ...
  'per_subject',  struct('emp',I_emp,'null_mu',I_mu,'null_sd',I_sd,'p',I_p) );
block.gearyC = struct( ...
  'empirical_mu', nanmean(C_emp), 'empirical_sd', nanstd(C_emp), ...
  'null_mu_mu',   nanmean(C_mu),  'null_mu_sd',   nanstd(C_mu), ...
  'p_sig_frac',   mean(C_p<=0.05,'omitnan'), 'n', sum(isfinite(C_emp)), ...
  'per_subject',  struct('emp',C_emp,'null_mu',C_mu,'null_sd',C_sd,'p',C_p) );
end

function blk = subject_value_dist_fast(allmaps_, surr_VSP, mask_, label_, VSAMP, NPERM, NBINS, FASTMODE)
S  = size(allmaps_,2);
P_ = size(surr_VSP,3);
ksv = nan(S,1); sdr = nan(S,1); mae = nan(S,1);
vidx = find(mask_ & any(isfinite(allmaps_),2));
if numel(vidx) > VSAMP, vidx = vidx(randperm(numel(vidx), VSAMP)); end
useP = min(P_, NPERM); pidx = 1:useP;
for s = 1:S
  e = allmaps_(vidx, s); e = e(isfinite(e));
  if numel(e) < 10, continue, end
  Y = reshape(surr_VSP(vidx, s, pidx), [], 1);
  y = Y(isfinite(Y));
  if numel(y) < 20, continue, end
  e_q = prctile_safe(e, [2 98]); y_q = prctile_safe(y, [2 98]);
  lo = min(e_q(1), y_q(1)); hi = max(e_q(2), y_q(2));
  if ~isfinite(lo) || ~isfinite(hi) || lo==hi
    lo = prctile_safe(e,2); hi = prctile_safe(e,98); if lo==hi, lo=lo-1; hi=hi+1; end
  end
  edges = linspace(lo, hi, NBINS+1);
  cdf_e = histcounts(e, edges, 'Normalization','cdf');
  cdf_y = histcounts(y, edges, 'Normalization','cdf');
  if FASTMODE, ksv(s) = max(abs(cdf_e - cdf_y));
  else,       ksv(s) = ks_distance_exact(e, y);
  end
  sdr(s) = std(e) / max(std(y), eps);
  qlev = 0.1:0.1:0.9;
  qe = quantile_safe(e, qlev);
  qy = invcdf_from_binned(cdf_y, edges, qlev);
  q25 = quantile_safe(e, 0.25);  q75 = quantile_safe(e, 0.75);
  IQR = max(q75 - q25, eps);
  mae(s) = mean(abs(qe - qy)) / IQR;
end
blk.label          = label_;
blk.KS_mu          = nanmean(ksv);             blk.KS_sd          = nanstd(ksv);
blk.SDratio_mu     = nanmean(sdr);             blk.SDratio_sd     = nanstd(sdr);
blk.DecMAE_IQR_mu  = nanmean(mae);             blk.DecMAE_IQR_sd  = nanstd(mae);
blk.n              = sum(isfinite(ksv) | isfinite(sdr) | isfinite(mae));
end

function [mu,sd,n,t,p] = fast_corr_emp_surr_group(t_emp_, surr_group_, mask_)
e  = t_emp_(mask_);
X  = surr_group_(mask_,:);
vals = fast_vec_corr(e, X);
[mu,sd,n,t,p] = one_sample_t(vals);
end

function out = fast_corr_emp_surr_subject(allmaps_, surr_VSP, ~, mask_, MINV, maxPerm)
% FAST, ROBUST subject emp vs surrogate correlation.
% Signature matches current calls:
%   (allmaps, surr_subj, adj, hemiMask, MINV, maxPerm)
% - Ignores adj (not needed for correlations).
% - Uses per-permutation valid overlap; collects up to maxPerm usable r's.
% - MINV: soft minimum overlap; will fall back to a tiny floor if needed.

% Accept both (V×S×P) and (S×V×P); convert to V×S×P
sz = size(surr_VSP);
if numel(sz) == 3 && sz(1) ~= size(allmaps_,1)
  % assume S×V×P
  surr_VSP = permute(surr_VSP, [2 1 3]);  % -> V×S×P
end

[V, S] = size(allmaps_);
P      = size(surr_VSP, 3);
K      = min(P, max(1, maxPerm));
vals   = nan(S,1);

if nargin < 5 || isempty(MINV), MINV = 50; end
MINV_HARD = 3;   % absolute floor (guardrail)

for s = 1:S
  e_full = allmaps_(:, s);
  if ~any(isfinite(e_full) & mask_), continue, end

  r_buf = nan(K,1);
  got   = 0;

  % Pass 1: stricter overlap (MINV)
  for p = 1:P
    if got >= K, break, end
    y = surr_VSP(:, s, p);
    m = mask_ & isfinite(e_full) & isfinite(y);
    nn = nnz(m);
    if nn >= MINV
      ee = e_full(m); yy = y(m);
      ee = ee - mean(ee); yy = yy - mean(yy);
      den = sqrt(sum(ee.^2) * sum(yy.^2));
      if den > 0
        got = got + 1;
        r_buf(got) = (ee' * yy) / den;
      end
    end
  end

  % Pass 2: relaxed overlap if nothing collected yet
  if got == 0
    for p = 1:P
      if got >= K, break, end
      y = surr_VSP(:, s, p);
      m = mask_ & isfinite(e_full) & isfinite(y);
      nn = nnz(m);
      if nn >= MINV_HARD
        ee = e_full(m); yy = y(m);
        ee = ee - mean(ee); yy = yy - mean(yy);
        den = sqrt(sum(ee.^2) * sum(yy.^2));
        if den > 0
          got = got + 1;
          r_buf(got) = (ee' * yy) / den;
        end
      end
    end
  end

  if got > 0
    vals(s) = mean(r_buf(1:got), 'omitnan');
  else
    vals(s) = NaN;
  end
end

[mu, sd, n, t, p] = one_sample_t(vals);
out = struct('mu',mu,'sd',sd,'n',n,'t',t,'p',p);
end

function vals = fast_vec_corr(e, X)
e = e(:);
ok = isfinite(e);
e = e(ok); X = X(ok,:);
e = e - mean(e);
X = X - mean(X,1);
den = sqrt(sum(e.^2)) * sqrt(sum(X.^2,1));
vals = (e' * X) ./ max(den, eps);
end

function [mu_r, sd_r, r_emp, p_two] = fast_cov_corr(t_emp, surr_group, cvg, vmask)
t = t_emp(vmask); c = double(cvg(vmask));
ok = isfinite(t) & isfinite(c); t=t(ok); c=c(ok);
t = t - mean(t); c = c - mean(c);
r_emp = (t'*c) / max(sqrt(sum(t.^2))*sqrt(sum(c.^2)), eps);
X = surr_group(vmask,:); X = X(ok,:); X = X - mean(X,1);
den = sqrt(sum(X.^2,1)) * sqrt(sum(c.^2));
rr = (c' * X) ./ max(den, eps);
mu_r = mean(rr,'omitnan'); sd_r = std(rr,'omitnan');
p_two = mean(abs(rr) >= abs(r_emp));
end

function [mu,sd,n,t,p] = fast_partial_corr_cov(t_emp, surr_group, cvg, vmask)
t = t_emp(vmask); c = double(cvg(vmask));
ok = isfinite(t) & isfinite(c); t=t(ok); c=c(ok);
t = t - mean(t); c = c - mean(c);
beta_t = (c' * t) / max(c' * c, eps);
rt = t - beta_t*c;
X = surr_group(vmask,:); X = X(ok,:); X = X - mean(X,1);
beta_x = (c' * X) / max(c' * c, eps);
RX = X - c * beta_x;
vals = fast_vec_corr(rt, RX);
[mu,sd,n,t,p] = one_sample_t(vals);
end

function [pp, KS_abs, KS_signed] = fast_pop_pvals_ranked(T, rand_ties, tail)
T = T(:); P = numel(T);
if nargin < 3, tail = 'high'; end
[Ts, ord] = sort(T, 'ascend');
[~,~,grp] = unique(Ts);
cnt = accumarray(grp, 1);
cum = cumsum(cnt);
gt_all = P - cum(grp);
eq_all = cnt(grp) - 1;
gt = zeros(P,1); eq = zeros(P,1);
gt(ord) = gt_all; eq(ord) = eq_all;
u = rand(P,1).*(rand_ties~=0);
if strcmpi(tail,'low')
  lt = cum(grp) - cnt(grp);
  pp = (lt + u.*eq + 1) / P;
else
  pp = (gt + u.*eq + 1) / P;
end
pps = sort(pp); unif = (1:P)'/P; diffs = unif - pps;
[~,imax] = max(abs(diffs));
KS_abs    = abs(diffs(imax));
KS_signed = diffs(imax);
end

function r = fast_pairwise_corr(x, Y, min_n)
if nargin < 3 || isempty(min_n), min_n = 3; end
x = x(:);
[n, k] = size(Y);
if numel(x) ~= n, error('fast_pairwise_corr:dimMismatch'); end
r = nan(k,1);
fin_x = isfinite(x);
for j = 1:k
  y = Y(:,j);
  m = fin_x & isfinite(y);
  nn = sum(m);
  if nn < min_n, r(j) = NaN; continue, end
  xx = x(m) - mean(x(m));
  yy = y(m) - mean(y(m));
  den = sqrt(sum(xx.^2) * sum(yy.^2));
  r(j) = (den>0) * ((xx' * yy) / max(den,eps)) + (den<=0)*NaN;
end
end

function qvals = invcdf_from_binned(cdf_vals, edges, probs)
% Invert a binned CDF (step function) at the requested probabilities.
% cdf_vals : length NBINS, monotonically nondecreasing in [0,1]
% edges    : length NBINS+1 bin edges
% probs    : vector of probabilities in [0,1]

cdf_vals = cdf_vals(:)'; 
edges    = edges(:)';

% work with bin midpoints for a smooth-ish inverse
bin_mid = (edges(1:end-1) + edges(2:end))/2;

% clamp and ensure final CDF reaches 1
cdf_vals = max(0, min(1, cdf_vals));
if cdf_vals(end) < 1, cdf_vals(end) = 1; end

% remove flat duplicates to make interp1 happy
[cdf_vals_u, uniq_idx] = unique(cdf_vals, 'stable');
bin_mid_u = bin_mid(uniq_idx);

% linear interpolation; extrap okay at the ends
qvals = interp1(cdf_vals_u, bin_mid_u, probs, 'linear', 'extrap');
end

function d = ks_distance_exact(x, y)
x = x(isfinite(x)); y = y(isfinite(y));
if isempty(x) || isempty(y), d = NaN; return, end
if exist('kstest2','file')==2
  [~,~,d] = kstest2(x, y); return
end
xx = sort(x); yy = sort(y);
vals = unique([xx; yy]);
Fx = arrayfun(@(v) sum(xx<=v)/numel(xx), vals);
Fy = arrayfun(@(v) sum(yy<=v)/numel(yy), vals);
d = max(abs(Fx - Fy));
end

function [mu,sd,n,t,p] = one_sample_t(x)
mu = nanmean(x); sd = nanstd(x); n = sum(isfinite(x));
if n>1 && sd>0
  t = mu / (sd / sqrt(n));
  if exist('tcdf','file')==2
    p = 2*tcdf(-abs(t), n-1);
  else
    p = 2*normcdf(-abs(t));
  end
else
  t = NaN; p = NaN;
end
end

function I = moran_I(x, W)
x = x(:); x = x - nanmean(x);
den = nansum(x.^2); if den<=0, I=NaN; return, end
S0  = full(sum(W(:))); if S0==0, I=NaN; return, end
I   = (numel(x)/S0) * (x'*(W*x)) / den;
end

function C = geary_C(x, W)
x = x(:); x = x - nanmean(x);
den = nansum(x.^2); if den<=0, C=NaN; return, end
S0  = full(sum(W(:))); if S0==0, C=NaN; return, end
d  = full(sum(W,2));
L  = spdiags(d,0,size(W,1),size(W,2)) - W;
num = 2*(x'*(L*x));
C   = ((numel(x)-1)/(2*S0)) * (num / den);
end

function q = prctile_safe(x, p)
x = x(isfinite(x));
if isempty(x), q = [NaN NaN]; return, end
if exist('prctile','file')==2, q = prctile(x, p);
else,                            q = quantile(x, p/100); end
end

function q = quantile_safe(x, p)
x = x(isfinite(x)); if isempty(x), q = NaN(size(p)); return, end
if exist('quantile','file')==2, q = quantile(x, p);
else,                           q = prctile(x, p*100); end
end

function print_spatial_row(lbl, blk, verbose)
dI = (blk.moranI.empirical_mu - blk.moranI.null_mu_mu) / max(blk.moranI.null_mu_sd, eps);
dC = (blk.gearyC.empirical_mu - blk.gearyC.null_mu_mu) / max(blk.gearyC.null_mu_sd, eps);
frac = blk.moranI.p_sig_frac;
n    = blk.moranI.n;
[ciL,ciU] = wilsonCI(0.05, n, 0.95);
lo = min(ciL,ciU); hi = max(ciL,ciU);
z  = (frac - 0.05) / sqrt(max(0.05*(1-0.05)/max(n,1), eps));
% fprintf(' %-4s  ΔI/SD=%.2f  ΔC/SD=%.2f  %%p≤.05=%.1f%% (target 5%%; 95%% %.1f–%.1f%%)  %s%s%s\n', ...
%   lbl, dI, dC, 100*frac, 100*lo, 100*hi, ...
%   badge_centered(dI,1.0,2.0), badge_centered(dC,1.0,2.0), badge_centered(z,1.5,3.0))

if verbose
  % fprintf('       I: emp μ=%.4f sd=%.4f | null μ=%.4f sd=%.4f\n', ...
    % blk.moranI.empirical_mu, blk.moranI.empirical_sd, blk.moranI.null_mu_mu, blk.moranI.null_mu_sd)
  % fprintf('       C: emp μ=%.4f sd=%.4f | null μ=%.4f sd=%.4f\n', ...
    % blk.gearyC.empirical_mu, blk.gearyC.empirical_sd, blk.gearyC.null_mu_mu, blk.gearyC.null_mu_sd)
end

% optional per-subject QA flagging
if isfield(blk.moranI,'per_subject')
  I_z = (blk.moranI.per_subject.emp - blk.moranI.per_subject.null_mu) ./ ...
        max(blk.moranI.per_subject.null_sd, eps);
  bad = find(abs(I_z) > 2);
  if ~isempty(bad)
    % fprintf('       Note: %d subject(s) |z|>2 for Moran''s I (e.g., #%s)\n', ...
      % numel(bad), strjoin(string(bad(1:min(5,end))), ', '));
  end
end
end


function print_value_dist_row(lbl, vd)
KS_OK=0.10; KS_WARN=0.20; SD_OK=0.10; SD_WARN=0.25; MAE_OK=0.15; MAE_WARN=0.30;
bKS = badge_low(vd.KS_mu, KS_OK, KS_WARN);
bSD = badge_around(vd.SDratio_mu, 1.0, SD_OK, SD_WARN);
bMA = badge_low(vd.DecMAE_IQR_mu, MAE_OK, MAE_WARN);
% fprintf(' %-4s value dist: KS=%.3f±%.3f %s  SD-ratio(emp/null)=%.3f±%.3f %s  decile-MAE/IQR=%.3f±%.3f %s (n=%d)\n', ...
%   lbl, vd.KS_mu, vd.KS_sd, bKS, vd.SDratio_mu, vd.SDratio_sd, bSD, vd.DecMAE_IQR_mu, vd.DecMAE_IQR_sd, bMA, vd.n)
end

function print_corr_row(lbl, s)
% fprintf(' %-4s  μ=%.4f ± %.4f (n=%d), t=%.2f, p=%.3g  %s\n', ...
%   lbl, s.mu, s.sd, s.n, s.t, s.p, badge_around(s.mu, 0, 0.02, 0.05));
end

function b = badge_around(x, target, tol_ok, tol_warn)
dx = abs(x - target);
if dx <= tol_ok, b='[PASS]'; elseif dx <= tol_warn, b='[WARN]'; else, b='[FLAG]'; end
end

function b = badge_low(x, thr_ok, thr_warn)
if x <= thr_ok, b='[PASS]'; elseif x <= thr_warn, b='[WARN]'; else, b='[FLAG]'; end
end

function b = badge_hiP(p, ok, warn)
if p >= ok, b='[PASS]'; elseif p >= warn, b='[WARN]'; else, b='[FLAG]'; end
end

function b = badge_centered(z, thr_ok, thr_warn)
az = abs(z);
if az <= thr_ok, b='[PASS]'; elseif az <= thr_warn, b='[WARN]'; else, b='[FLAG]'; end
end

function line_metric(name, val, extra, badge)
% fprintf(' %-18s %10.4f   %s   %s\n', name, val, extra, badge)
end

function line_note(txt)
% fprintf('   %s\n', txt)
end

function s = yn(tf)
if tf, s='yes'; else, s='no'; end
end

function [ciL, ciU] = wilsonCI(p, n, conf)
if nargin<3 || isempty(conf), conf = 0.95; end
if ~isfinite(n) || n <= 0 || ~isfinite(p)
  ciL = 0; ciU = 1; return
end
p = max(0, min(1, p));
if exist('norminv','file') == 2
  z = norminv(1 - (1 - conf)/2);
else
  z = sqrt(2) * erfcinv(2 * (1 - (1 - conf)/2));
end
den     = 1 + (z^2)/n;
centre  = (p + (z^2)/(2*n)) / den;
halfwid = (z/den) * sqrt( (p*(1-p)/n) + (z^2)/(4*n^2) );
ciL = max(0, centre - halfwid);
ciU = min(1, centre + halfwid);
end
