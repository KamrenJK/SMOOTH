function export_hgb_traces_hires(chans, sigma_ms)
% ------------------------------------------------------------------------------------------------
% High-time-resolution HGB traces for the Figure 2A inset, for display only. Computes TWO versions
% of the same four channels so they can be compared side by side.
%
% WHY.  The traces cached by export_hgb_traces.m come straight from the proc_smooth.m pipeline that
% produced HFBcorr: mtmconvol / DPSS, foi 70:10:150, t_ftimwin 300 ms, toi in 100 ms steps. That is
% the right series to compute the correlations on -- it is the one the paper's numbers refer to --
% but it is a poor thing to plot. Two things limit it, and they are independent:
%   * the 300 ms taper window IS the temporal resolution; nothing faster survives
%   * the 100 ms toi step samples that at 10 Hz, giving ~100 knots in a 10 s clip
% Because consecutive 300 ms boxcars overlap by only 200 ms, a burst entering or leaving one
% produces a step-like jump. The jaggedness in the 10 Hz panel is mostly that, not fast dynamics.
%
% METHOD 'mtm'  -- same decomposition, finer sampling.
%   Identical to proc_smooth (DPSS, foi 70:10:150, t_ftimwin 300 ms, tapsmofrq 10, z per channel x
%   frequency over the whole recording, then averaged across frequency), except toi steps 10 ms
%   instead of 100 ms, followed by a SIGMA_MS Gaussian. This isolates the SAMPLING limit: at 10 ms
%   steps neighbouring windows share 97% of their data, so the curve becomes smooth, but the
%   effective resolution is still 300 ms and a 25 ms Gaussian is far narrower than that, so it
%   barely acts. Expect a cleaner-looking line carrying the same information.
%
% METHOD 'hil'  -- sub-band Hilbert envelope, which removes the window entirely.
%   1. same preprocessing (bipolar, 1-198 Hz, 58/118/178 Hz line notches, 500 Hz, 56 trials
%      concatenated to one 1652 s block) -- unchanged, so the input is identical
%   2. eight 10 Hz sub-bands, 70-80 ... 140-150 Hz, zero-phase Butterworth (filtfilt)
%   3. abs(hilbert()) per sub-band
%   4. divide each sub-band by its own whole-recording mean, then average. Same purpose as the
%      per-frequency z-scoring above -- without it the 70-80 Hz band dominates through 1/f -- but
%      involving no time window
%   5. SIGMA_MS Gaussian (default 25 ms; FWHM = 2.355*sigma = 59 ms), then z against the WHOLE
%      recording so the z axis stays comparable to the 10 Hz panel
%   Output at 500 Hz, resolution set by the band-pass and the smoother rather than by a boxcar.
%
% HONESTY.  Neither version is what HFBcorr was computed on, so the r values printed on the figure
% stay the stored multitaper ones (that is what the berry map is painted with). This script
% recomputes the whole-recording correlations under both methods and reports all three: if the
% decompositions disagreed materially, the display traces would be telling a different story from
% the map.
%
% Only the requested channels are processed -- the whole recording at 500 Hz for all 180 would be
% ~1 GB -- and only a padded window around each candidate display segment is cached.
%
% Usage: export_hgb_traces_hires                              % default channels, 25 ms
%        export_hgb_traces_hires({'TG36-TG37',...}, 25)
% Requires the raw data (as export_hgb_traces.m), including $SMOOTH_RAWSUBJECT = P19's raw-data
% subject ID (not distributed). Everything written uses the pseudonym P19.
% Output: paper/data/hgb_traces_P19_hires.mat
% ------------------------------------------------------------------------------------------------
if nargin < 1 || isempty(chans)
    chans = {'TG36-TG37','TG28-TG36','TG35-TG36','TG1-TG2'};
end
if nargin < 2 || isempty(sigma_ms), sigma_ms = 25; end

ECOG   = getenv('SMOOTH_RAWDATA'); assert(~isempty(ECOG), 'Set SMOOTH_RAWDATA to the raw-data root (not distributed).');
RAWSUB = getenv('SMOOTH_RAWSUBJECT'); assert(~isempty(RAWSUB), 'Set SMOOTH_RAWSUBJECT to P19''s raw-data subject ID (not distributed).');
SUB    = 'P19'; CSUB = 18;                              % pseudonym (demo/source.mat 19) for all outputs
FOI    = 70:10:150;  TFTIMWIN = 0.3;  TAPSMOFRQ = 10;   % proc_smooth values, unchanged
TOISTEP = 0.010;                                        % 100 ms -> 10 ms
BANDS  = [70:10:140; 80:10:150].';                      % eight 10 Hz sub-bands
PADSEC = 20;                                            % cached either side of each display window
WINDOWS = struct('name',{'z5','z8'}, 'display_s',{[306 316],[842 852]});

addpath(fullfile(fileparts(mfilename('fullpath')),'..','..'));
PREPROC = getenv('SMOOTH_PREPROC'); assert(~isempty(PREPROC), 'Set SMOOTH_PREPROC to the preprocessing-code directory (not distributed).');
addpath(ECOG); addpath(PREPROC);
paths = smooth_setup('quiet');

% ---------------------------------------------------------------- preprocessing (unchanged)
fprintf('preprocessing %s ...\n', SUB);
data = preproc_tcg(['subject_' RAWSUB], 1, 1, 1);
cfg = []; cfg.resamplefs = 500; data = ft_resampledata(cfg, data);
tmp = rmfield(data, {'trial','time'});
tmp.trial = {cat(2, data.trial{:})};
tmp.time  = {0:1/tmp.fsample:size(tmp.trial{1},2)/tmp.fsample - 1/tmp.fsample};
[tf, ci] = ismember(string(chans(:)), string(tmp.label(:)));
if ~all(tf), error('export_hgb_traces_hires:chan','missing: %s', strjoin(chans(~tf),', ')); end
tmp.label = tmp.label(ci); tmp.trial{1} = tmp.trial{1}(ci,:);
if isfield(tmp,'elec'), tmp = rmfield(tmp,'elec'); end
n  = numel(chans); fs = tmp.fsample; t = tmp.time{1};
fprintf('  %d channels, %.1f s at %g Hz\n', n, t(end), fs);

% ---------------------------------------------------------------- method 1: mtm at 10 ms steps
fprintf('mtmconvol, %g ms steps ...\n', TOISTEP*1000);
cfg = [];
cfg.output = 'pow'; cfg.method = 'mtmconvol'; cfg.taper = 'dpss';
cfg.foi = FOI; cfg.tapsmofrq = TAPSMOFRQ;
cfg.t_ftimwin = TFTIMWIN*ones(numel(FOI),1);
cfg.toi = t(1):TOISTEP:t(end);
cfg.keeptrials = 'yes'; cfg.pad = 'nextpow2';
tfr = ft_freqanalysis(cfg, tmp);
for c = 1:numel(tfr.label)                                  % z per channel x frequency, as proc_smooth
    for f = 1:numel(tfr.freq)
        ts = tfr.powspctrm(:,c,f,:);
        tfr.powspctrm(:,c,f,:) = (ts - mean(ts(:),'omitnan')) / std(ts(:),'omitnan');
    end
end
cfg = []; cfg.avgoverfreq = 'yes'; cfg.nanmean = 'yes';
Ymtm  = squeeze(ft_selectdata(cfg, tfr).powspctrm);
tmtm  = tfr.time(:).'; fsmtm = 1/TOISTEP;
Ymtm  = gsmooth(Ymtm, sigma_ms/1000*fsmtm);
fprintf('  %d x %d at %g Hz\n', size(Ymtm,1), size(Ymtm,2), fsmtm);

% ---------------------------------------------------------------- method 2: sub-band Hilbert
fprintf('sub-band hilbert ...\n');
X = tmp.trial{1}; env = zeros(n, size(X,2));
for b = 1:size(BANDS,1)
    [bb, aa] = butter(4, BANDS(b,:)/(fs/2), 'bandpass');
    e = abs(hilbert(filtfilt(bb, aa, X.'))).';              % filtfilt/hilbert work down columns
    env = env + e ./ mean(e, 2);                            % per-band mean normalisation
    fprintf('  %3d-%3d Hz\n', BANDS(b,1), BANDS(b,2));
end
Yhil = gsmooth(env/size(BANDS,1), sigma_ms/1000*fs);
Yhil = (Yhil - mean(Yhil,2,'omitnan')) ./ std(Yhil,0,2,'omitnan');
thil = t; fshil = fs;

% ---------------------------------------------------------------- compare against the stored HFBcorr
S = load(fullfile(paths.data,'corrsource.mat'));
ref = double(S.corrsource{CSUB}.HFBcorr); rl = string(S.corrsource{CSUB}.label(:));
[~, ri] = ismember(string(chans(:)), rl);
rstored = ref(ri(1), ri);
rmtm = rowcorr(Ymtm); rhil = rowcorr(Yhil);
ep   = double(S.corrsource{CSUB}.elec.nativechanpos);
dmm  = vecnorm(ep(ri,:) - ep(ri(1),:), 2, 2).';
fprintf('\n  whole-recording correlation to the seed, by decomposition\n');
fprintf('  %-12s %7s %11s %11s %11s\n','channel','d (mm)','stored','mtm 10ms','hilbert');
for k = 1:n
    fprintf('  %-12s %7.1f %11.3f %11.3f %11.3f\n', chans{k}, dmm(k), rstored(k), rmtm(k), rhil(k));
end
fprintf('  ("stored" is the 100 ms multitaper series HFBcorr and the berry map come from)\n\n');

% ---------------------------------------------------------------- cache padded windows, both methods
out = struct('chans', {chans}, 'd_mm', dmm, 'r_stored', rstored, 'r_mtm', rmtm, 'r_hilbert', rhil, ...
             'sigma_ms', sigma_ms, 'bands', BANDS, 'foi', FOI, 't_ftimwin', TFTIMWIN, ...
             'toistep_s', TOISTEP, 'pad_s', PADSEC, 'subject', SUB, 'corrsource_idx', CSUB, ...
             'recording_s', t(end));
for m = {{'mtm', tmtm, Ymtm, fsmtm}, {'hil', thil, Yhil, fshil}}
    nm = m{1}{1}; tt = m{1}{2}; YY = m{1}{3}; ff = m{1}{4};
    C = struct();
    for k = 1:numel(WINDOWS)
        d = WINDOWS(k).display_s; sel = tt >= d(1)-PADSEC & tt <= d(2)+PADSEC;
        C.(WINDOWS(k).name) = struct('t', tt(sel), 'Y', YY(:,sel), 'display_s', d);
        fprintf('  %s/%s: %.0f-%.0f s, %d samples at %g Hz\n', nm, WINDOWS(k).name, ...
                tt(find(sel,1)), tt(find(sel,1,'last')), nnz(sel), ff);
    end
    out.(nm) = struct('clips', C, 'fs', ff);
end
save(fullfile(paths.data,'hgb_traces_P19_hires.mat'), '-struct', 'out');
fprintf('wrote %s\n', fullfile(paths.data,'hgb_traces_P19_hires.mat'));
end

% ==================================== helpers ====================================
function Y = gsmooth(Y, sd)
% Gaussian smoothing along time; edges where the kernel would run off are set NaN rather than
% silently tapered, so nothing partial is ever plotted or correlated.
if sd <= 0, return; end
L = 2*ceil(3*sd)+1; g = exp(-((1:L)-(L+1)/2).^2/(2*sd^2)); g = g/sum(g);
for k = 1:size(Y,1), Y(k,:) = conv(Y(k,:), g, 'same'); end
Y(:,[1:floor(L/2), end-floor(L/2)+1:end]) = NaN;
end

function r = rowcorr(Y)
ok = all(isfinite(Y),1); C = corr(Y(:,ok)'); r = C(1,:);
end
