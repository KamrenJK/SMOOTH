function out = seedxcorr_panel(Y, fs, col, chans, maxlag_s, fn, ttl)
% ------------------------------------------------------------------------------------------------
% Lagged correlation of every contact against the SEED contact (row 1 of Y), all superimposed, one
% line per contact in the colour it carries in the seed-correlation berry map. Companion panel to
% the trace insets (genFig2_seedcorr_traces / _hires).
%
% The seed's own curve is therefore its autocorrelation (1 at zero lag by construction); every other
% curve is a cross-correlation with the seed. Normalisation is xcorr's 'coeff', so the value at zero
% lag is exactly the Pearson correlation between that contact and the seed over this segment -- the
% same quantity the berry map colours, evaluated on this clip.
%
% Lags are two-sided. Cross-correlation is not symmetric, so a peak displaced from zero would mean
% one contact systematically leads the other; keeping both signs makes that visible instead of
% hiding it under a one-sided plot. Sign convention: POSITIVE lag = that contact lags the seed.
%
% Y         nchan x ntime, seed in row 1
% fs        sampling rate (Hz)
% col       nchan x 3 line colours (from the berry map's viridis ramp)
% maxlag_s  largest lag to plot, each side
% ttl       short provenance string, printed to stdout only -- nothing is baked into the panel
%
% Returns a struct with the zero-lag correlation, the peak correlation and the lag at which it
% occurs, per contact.
% Output is a vector PDF on a white background, like the other 2-D panels here.
% ------------------------------------------------------------------------------------------------
n = size(Y,1); ml = round(maxlag_s*fs);
ok = all(isfinite(Y),1); Ys = Y(:,ok);
Ys = Ys - mean(Ys,2);
seed = Ys(1,:).';

C = nan(n, 2*ml+1);
for k = 1:n
    [c, lags] = xcorr(Ys(k,:).', seed, ml, 'coeff');       % +lag: contact k lags the seed
    C(k,:) = c(:).';
end
lag = lags(:).'/fs;

f  = figure('color','w','units','centimeters','position',[2 2 9.5 7.5]);
ax = axes('Parent',f,'units','centimeters','position',[1.9 1.5 7.2 5.6]); hold(ax,'on'); box(ax,'off');
yline(ax, 0, '-', 'Color',[.75 .75 .75], 'LineWidth',0.5);
xline(ax, 0, '-', 'Color',[.75 .75 .75], 'LineWidth',0.5);
for k = 1:n, plot(ax, lag, C(k,:), 'Color', col(k,:), 'LineWidth', 1.4); end
set(ax,'FontSize',9,'TickDir','out','TickLength',[0.012 0.012]);
xlim(ax,[-maxlag_s maxlag_s]); ylim(ax,[min(-0.1, min(C(:))*1.1) 1.02]);
xlabel(ax,'lag (s)','FontSize',10); ylabel(ax,'correlation with seed','FontSize',10);
exportgraphics(f, fn, 'ContentType','vector'); close(f);

[pk, pi_] = max(C, [], 2); pklag = lag(pi_);
out = struct('lag', lag, 'C', C, 'r0', C(:,ml+1).', 'peak', pk.', 'peak_lag', pklag);
fprintf('  seed cross-correlation (%s):\n', ttl);
fprintf('    %-12s %9s %9s %10s\n','channel','r(lag 0)','peak r','peak lag');
for k = 1:n
    fprintf('    %-12s %9.3f %9.3f %9.3f s\n', chans{k}, C(k,ml+1), pk(k), pklag(k));
end
fprintf('  wrote %s\n', fn);
end
