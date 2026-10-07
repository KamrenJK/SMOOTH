function genFig2_seedcorr_traces_hires(method, window, clipsec)
% ------------------------------------------------------------------------------------------------
% Trace inset for Figure 2A, drawn from the high-time-resolution HGB reprocessing
% (export_hgb_traces_hires.m). Same four contacts, same colours, same layout as
% genFig2_seedcorr_traces.m -- only the underlying time series differs.
%
%   method 'mtm'  proc_smooth's own decomposition (DPSS, foi 70:10:150, 300 ms taper window) but
%                 stepped every 10 ms instead of 100 ms, then smoothed with a 25 ms Gaussian.
%                 Isolates the SAMPLING limit: neighbouring windows now share 97% of their data, so
%                 the step-like jaggedness of the 10 Hz series disappears, but the effective
%                 resolution is still the 300 ms window and the 25 ms Gaussian is far narrower than
%                 that, so it barely acts. 100 Hz output.
%   method 'hil'  eight 10 Hz sub-bands, Hilbert envelope, per-band mean normalised, 25 ms Gaussian,
%                 z over the whole recording. No time window at all, so the resolution is set by the
%                 band-pass and the smoother. 500 Hz output.
%
%   window 'z5'   306-316 s   (the quiet window: peak |z| 4.6 in the 10 Hz series)
%   window 'z8'   842-852 s   (peak |z| 6.3; the clearer near-pair covariation)
%
% Colours are sampled from viridis at each contact's STORED HFBcorr correlation on clim [0 1] -- the
% same numbers the berry map is painted with -- not at the correlation of the reprocessed trace.
% Both are printed on every run. The three decompositions agree on the ordering and roughly on the
% magnitudes (seed row: stored 1.000/0.517/0.300/0.031, mtm 1.000/0.569/0.367/0.052, Hilbert
% 1.000/0.542/0.471/0.022); the Hilbert envelope reads the 10.5 mm contact noticeably higher, so any
% quoted correlation should be the stored value, not a trace-derived one.
%
% Use as:
%   genFig2_seedcorr_traces_hires                      % hil, z5, 10 s
%   genFig2_seedcorr_traces_hires('mtm','z8',10)
%
% Data:   paper/data/hgb_traces_P19_hires.mat
% Output: paper/figs/2_seedcorr/hgbtraces_P19_<window>_<method>_<clip>s.pdf  (vector, white bg)
% ------------------------------------------------------------------------------------------------
if nargin < 1 || isempty(method),  method  = 'hil'; end
if nargin < 2 || isempty(window),  window  = 'z5';  end
if nargin < 3 || isempty(clipsec), clipsec = 10;    end

CLIM   = [0 1];                   % must match the berry map's clim or the colours will not agree
LW     = 0.9;
GAP    = 1.7;
XCLAG  = 1.0;                     % largest lag (s), each side, on the seed cross-correlation panel

addpath(fullfile(fileparts(mfilename('fullpath')),'..','..'));
paths = smooth_setup('quiet');
D = load(fullfile(paths.data,'hgb_traces_P19_hires.mat'));
outdir = fullfile(paths.figs,'2_seedcorr'); if ~exist(outdir,'dir'), mkdir(outdir); end
gi = get(groot,{'defaultTextInterpreter'});
set(groot,'defaultTextInterpreter','none','defaultAxesFontName','Arial','defaultTextFontName','Arial');
restore = onCleanup(@() set(groot,'defaultTextInterpreter',gi{1})); %#ok<NASGU>

M = D.(method); C = M.clips.(window);
chans = cellstr(D.chans); n = numel(chans);
r = double(D.r_stored(:)).'; dmm = double(D.d_mm(:)).';
if strcmp(method,'mtm'), rnew = double(D.r_mtm); else, rnew = double(D.r_hilbert); end

% Centre the requested clip inside the cached (padded) segment.
t0 = double(C.t(:)).'; Y0 = double(C.Y); disp_s = double(C.display_s);
mid = mean(disp_s); sel = t0 >= mid-clipsec/2 & t0 < mid+clipsec/2;
t = t0(sel) - t0(find(sel,1)); Y = Y0(:,sel);

cmap = viridis; nc = size(cmap,1);
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
              [ "seed"; compose('%.0f mm', dmm(2:end).') ], r(:));
set(ax,'YTick',flipud(off),'YTickLabel',flipud(lbl),'FontSize',9,'TickDir','out','TickLength',[0.006 0.006]);
xlim(ax,[t(1) t(end)]); ylim(ax,[min(off)-1.6*sc, max(off)+1.9*sc]);
xlabel(ax,'time (s)','FontSize',10);

xb = t(end) - 0.03*(t(end)-t(1));
plot(ax, [xb xb], off(1)+[0.3 0.3+sc], 'k-', 'LineWidth',1.2);
text(ax, xb - 0.015*(t(end)-t(1)), off(1)+0.3+sc/2, sprintf('%.0f z', sc), 'FontSize',8, ...
     'HorizontalAlignment','right','VerticalAlignment','middle');

fn = fullfile(outdir, sprintf('hgbtraces_P19_%s_%s_%gs.pdf', window, method, clipsec));
exportgraphics(f, fn, 'ContentType','vector'); close(f);

fprintf('\n%s | %s | %s | %g s from %.0f-%.0f s | %g Hz, %g ms gaussian\n', ...
        D.subject, method, window, clipsec, t0(find(sel,1)), t0(find(sel,1,'last')), M.fs, D.sigma_ms);
fprintf('  %-12s %7s %11s %11s\n','channel','d (mm)','r(stored)','r(this)');
for k = 1:n, fprintf('  %-12s %7.1f %11.3f %11.3f\n', chans{k}, dmm(k), r(k), rnew(k)); end
fprintf('wrote %s\n', fn);

% Lagged correlation against the seed (see seedxcorr_panel.m). Computed on the FULL cached (padded)
% segment as well as on the displayed clip: the lag profile is a property of the signal, not of the
% window chosen for display, and 50 s gives a far more stable estimate -- but the clip version is
% what the plotted traces actually show, so any disagreement stays visible rather than hidden.
seedxcorr_panel(Y0, M.fs, col, chans, XCLAG, ...
    fullfile(outdir, sprintf('hgbxcorr_P19_%s_%s.pdf', window, method)), ...
    sprintf('%s/%s, %.0f s padded segment', method, window, t0(end)-t0(1)));
seedxcorr_panel(Y,  M.fs, col, chans, XCLAG, ...
    fullfile(outdir, sprintf('hgbxcorr_P19_%s_%s_clip%gs.pdf', window, method, clipsec)), ...
    sprintf('%s/%s, %g s displayed clip', method, window, clipsec));
end
