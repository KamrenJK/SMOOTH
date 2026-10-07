function genFig2_seedcorr_traces(clipsec, sigma_s, chans, zmax)
% ------------------------------------------------------------------------------------------------
% Trace inset for Figure 2A: a short clip of the HGB power time course for four contacts of
% P19's TG grid, each drawn in the SAME colour it carries in the seed-correlation berry map
% (genFig2_seedcorr, viridis on clim [0 1]). The contacts step outward from the seed, so read top to
% bottom the traces go from identical, to visibly co-varying, to unrelated -- the spatial
% autocorrelation the berry map shows as a colour gradient, shown here as signal.
%
%   TG36-TG37   seed          r = 1.000
%   TG28-TG36   d =  7.6 mm   r = 0.517    adjacent
%   TG35-TG36   d = 10.5 mm   r = 0.300    intermediate
%   TG1-TG2     d = 48.7 mm   r = 0.031    far corner of the same 8 x 8 grid
%
% Distance and correlation both fall monotonically. The steps are uneven in distance because r
% collapses fast (HFB FWHM is ~13 mm in this cohort): past ~20 mm every contact sits at r < 0.07, so
% a ladder evenly spaced in distance would have three indistinguishable bottom rows. An alternative
% middle step with wider spacing is TG42-TG43, at 23 mm, r = 0.07.
%
% Use as:
%   genFig2_seedcorr_traces                     % 10 s clip, no smoothing
%   genFig2_seedcorr_traces(10, 0.3)            % 0.3 s Gaussian smoothing
%   genFig2_seedcorr_traces(10, 0, {'TG36-TG37','TG28-TG36','TG42-TG43','TG1-TG2'})
%
% SMOOTHING is off by default and is a display choice only. It also RAISES apparent correlation --
% the noise it removes is independent across contacts -- so the clip's correlations are printed on
% every run next to the whole-recording values, and the colours (and any caption) stay tied to the
% latter. Sigma is in seconds and stamped into the filename.
%
% WINDOW.  A 10 s clip is only ~100 samples, so a single window's correlations are very noisy: across
% windows the distant contact ranges from about -0.85 to +0.8 around a median near zero. Scoring
% windows by "near minus far" therefore does not find a clear example of decay, it finds the window
% where the distant contact is most strongly ANTI-correlated -- structure that is not there. So the
% search keeps only windows where (a) no contact exceeds ZMAX and (b) the intermediate and far
% contacts stay within TOL of their whole-recording correlations, and among those takes the one
% where the near pair is strongest. Treat the clip as an illustration; the quantitative claim rests
% on the full 1652 s recording.
%
% ZMAX rejects windows containing sharp shared transients (they are local, reaching at most 15 of
% 180 channels at once, not a recording-wide artifact). Those excursions
% run to z = 50 and read as artifact even where they are not, so they are kept out of the panel.
% They are NOT what drives the correlation -- excluding every |z| > 3 sample from the WHOLE recording,
% 3.6% of it, moves the seed row only from 0.517 / 0.300 / 0.031 to 0.435 / 0.235 / 0.015 -- so
% dropping them costs the figure nothing it was relying on. Available windows by threshold:
%   ZMAX  3: 59 windows, best near-pair r = 0.72     ZMAX  6: 587, r = 0.76
%   ZMAX  4: 261,        r = 0.73                    ZMAX  8: 746, r = 0.85
%   ZMAX  5: 427,        r = 0.76                    Inf:    798, r = 0.86 (max |z| 11.4)
%
% Data: paper/data/hgb_traces_P19.mat (from export_hgb_traces.m, which reruns the proc_smooth.m
%   pipeline and verifies it reproduces the stored HFBcorr exactly -- it does, to 0.00e+00).
% Output: paper/figs/2_seedcorr/hgbtraces_P19_<seed>_<clip>s_sd<sigma>.pdf   (vector, white
%   background -- 2-D panels are vector PDFs here; only brain maps are transparent PNGs)
% ------------------------------------------------------------------------------------------------
if nargin < 1 || isempty(clipsec), clipsec = 10;  end
if nargin < 2,                     sigma_s = 0;   end
if nargin < 3 || isempty(chans)
    chans = {'TG36-TG37','TG28-TG36','TG35-TG36','TG1-TG2'};
end
if nargin < 4 || isempty(zmax), zmax = 5; end

CLIM = [0 1];                     % must match the berry map's clim or the colours will not agree
TOL  = 0.15;                      % how far the mid/far contacts may drift from their recording r
LW   = 1.0;
GAP  = 1.7;                       % row spacing, in units of the pooled robust trace scale. Generous
                                  % because the shared HGB bursts are tall: tighter spacing lets them
                                  % overflow into the neighbouring row, which is exactly where the eye
                                  % is supposed to be comparing.

addpath(fullfile(fileparts(mfilename('fullpath')),'..','..'));
paths = smooth_setup('quiet');
D = load(fullfile(paths.data,'hgb_traces_P19.mat'));
outdir = fullfile(paths.figs,'2_seedcorr'); if ~exist(outdir,'dir'), mkdir(outdir); end
gi = get(groot,{'defaultTextInterpreter'});
set(groot,'defaultTextInterpreter','none','defaultAxesFontName','Arial','defaultTextFontName','Arial');
restore = onCleanup(@() set(groot,'defaultTextInterpreter',gi{1})); %#ok<NASGU>

all_lab = string(D.chans(:));
[tf, ci] = ismember(string(chans(:)), all_lab);
if ~all(tf), error('genFig2_seedcorr_traces:chan','missing: %s', strjoin(chans(~tf),', ')); end
t0 = double(D.t(:)).'; Y0 = double(D.Y(ci,:)); fs = double(D.fs); n = numel(chans);
r  = double(D.hfbcorr(ci(1), ci));                 % whole-recording r to the seed, as stored

% Contact separation, for the printout -- taken from corrsource so it is the same geometry the
% berry map snaps to the surface.
S = load(fullfile(paths.data,'corrsource.mat'));
ep = double(S.corrsource{D.corrsource_idx}.elec.nativechanpos);
el = string(S.corrsource{D.corrsource_idx}.label(:));
[~, ei] = ismember(string(chans(:)), el);
dmm = vecnorm(ep(ei,:) - ep(ei(1),:), 2, 2).';

% ---------------------------------------------------------------- smooth (off by default)
if sigma_s > 0
    sd = sigma_s*fs; L = 2*ceil(3*sd)+1;
    g  = exp(-((1:L)-(L+1)/2).^2/(2*sd^2)); g = g/sum(g);
    Ys = nan(size(Y0));
    for k = 1:n, Ys(k,:) = conv(Y0(k,:), g, 'same'); end
    Ys(:,[1:floor(L/2), end-floor(L/2)+1:end]) = NaN;   % drop partial-kernel edges
else
    Ys = Y0;
end

% ---------------------------------------------------------------- pick the window
w = round(clipsec*fs); step = max(1, round(1*fs));      % slide in 1 s steps
starts = 1:step:(size(Ys,2)-w);
rw = nan(numel(starts), n); zpk = nan(numel(starts),1);
for k = 1:numel(starts)
    seg = Ys(:, starts(k):starts(k)+w-1); ok = all(isfinite(seg),1);
    if nnz(ok) < 0.95*w, continue; end
    cc = corr(seg(:,ok)'); rw(k,:) = cc(1,:); zpk(k) = max(abs(seg(:,ok)),[],'all');
end
valid    = isfinite(rw(:,2));
quiet    = zpk <= zmax;                                 % no sharp shared transients
faithful = all(abs(rw(:,3:end) - r(3:end)) <= TOL, 2) & valid & quiet;
if any(faithful)
    score = -inf(numel(starts),1); score(faithful) = rw(faithful,2);
    [~, best] = max(score);
else
    warning('genFig2_seedcorr_traces:tol', ...
        ['no %g s window is both quiet (|z| <= %g) and holds mid+far within %.2f of their recording ' ...
         'r; using the most representative'], clipsec, zmax, TOL);
    [~, best] = min(max(abs(rw(:,2:end) - r(2:end)), [], 2));
end
sel = starts(best):starts(best)+w-1;
t = t0(sel) - t0(sel(1)); Y = Ys(:,sel);

% ---------------------------------------------------------------- plot
cmap = viridis; nc = size(cmap,1);                      % same sampling ft_plot_cloud uses
cidx = min(max(round((r - CLIM(1))/(CLIM(2)-CLIM(1)) * (nc-1)) + 1, 1), nc);
col  = cmap(cidx,:);

% One shared vertical scale for all rows -- rescaling each row to fit would hide the amplitude
% differences between contacts and make the distant one look like a copy of the seed.
sc  = median(median(abs(Y - median(Y,2,'omitnan')), 2, 'omitnan')) * 6;
off = (n-1:-1:0)' * GAP * sc;

f  = figure('color','w','units','centimeters','position',[2 2 17.2 8]);
ax = axes('Parent',f,'units','centimeters','position',[5.0 1.5 11.7 6.1]); hold(ax,'on'); box(ax,'off');
for k = 1:n, plot(ax, t, Y(k,:) + off(k), 'Color', col(k,:), 'LineWidth', LW); end
lbl = compose('%s   %s   r = %.2f', string(chans(:)), ...
              [string(sprintf('%s','seed')); compose('%.0f mm', dmm(2:end).')], r(:));
set(ax,'YTick',flipud(off),'YTickLabel',flipud(lbl),'FontSize',9,'TickDir','out','TickLength',[0.006 0.006]);
xlim(ax,[t(1) t(end)]); ylim(ax,[min(off)-1.6*sc, max(off)+1.9*sc]);
xlabel(ax,'time (s)','FontSize',10);

xb = t(end) - 0.03*(t(end)-t(1));
plot(ax, [xb xb], off(1)+[0.3 0.3+sc], 'k-', 'LineWidth',1.2);
text(ax, xb - 0.015*(t(end)-t(1)), off(1)+0.3+sc/2, sprintf('%.0f z', sc), 'FontSize',8, ...
     'HorizontalAlignment','right','VerticalAlignment','middle');

fn = fullfile(outdir, sprintf('hgbtraces_P19_%s_%gs_sd%g_z%g.pdf', chans{1}, clipsec, sigma_s, zmax));
exportgraphics(f, fn, 'ContentType','vector'); close(f);

fprintf('\n%s | %g s clip, sigma %g s, ZMAX %g | window %.0f-%.0f s of %.0f s (peak |z| %.1f) | %d of %d windows usable\n', ...
        D.subject, clipsec, sigma_s, zmax, t0(sel(1)), t0(sel(end)), D.recording_s, zpk(best), nnz(faithful), nnz(valid));
fprintf('  %-12s %7s %9s %9s %9s   colour\n','channel','d (mm)','r(full)','r(clip)','r(median)');
for k = 1:n
    fprintf('  %-12s %7.1f %9.3f %9.3f %9.3f   [%.2f %.2f %.2f]\n', ...
            chans{k}, dmm(k), r(k), rw(best,k), median(rw(:,k),'omitnan'), col(k,1), col(k,2), col(k,3));
end
fprintf('wrote %s\n', fn);

% Lagged correlation against the seed (see seedxcorr_panel.m) -- one line per contact in its
% berry-map colour; the seed's own curve is its autocorrelation.
seedxcorr_panel(Y, fs, col, chans, 1.0, ...
    fullfile(outdir, sprintf('hgbxcorr_P19_%s_%gs_sd%g_z%g.pdf', chans{1}, clipsec, sigma_s, zmax)), ...
    sprintf('10 Hz multitaper, %g s clip', clipsec));
end
