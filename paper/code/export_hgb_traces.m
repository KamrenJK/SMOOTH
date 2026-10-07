function export_hgb_traces
% ------------------------------------------------------------------------------------------------
% Regenerate the HGB (high-frequency broadband) power time courses for participant P19 and cache them for the
% trace inset of Figure 2A (see genFig2_seedcorr).
%
% corrsource.mat stores only the channel-by-channel correlation MATRICES, not the time courses they
% were computed from, so the traces have to be recomputed. This script reruns the exact pipeline
% from proc_smooth.m that produced `HFBcorr` -- bipolar referencing, resample to 500 Hz, concatenate
% trials, mtmconvol dpss 70:10:150 Hz (tapsmofrq 10, 300 ms window, 100 ms steps), z-score each
% channel x frequency over the whole recording, then average across frequency -- and VERIFIES the
% result by recomputing corr(pow_hg') and comparing it against the stored HFBcorr. If that check
% does not pass, the traces are not the ones behind the berry map and the script errors out.
%
% ALL 180 channels are cached for the WHOLE recording (1652 s at 10 Hz, ~24 MB), not a pre-cut clip
% for a fixed channel set: which contacts to show, which window, and whether to smooth are all
% plotting decisions, and caching everything means changing any of them never requires rerunning the
% ~3 min preprocessing + TFR. Channel and window selection live in genFig2_seedcorr_traces.m.
%
% Requires the raw data under $SMOOTH_RAWDATA/data/<raw subject ID> and preproc_tcg.m; not runnable
% from the published repository alone, which is why the output is cached. The raw-data subject ID
% for P19 is not distributed: set it in $SMOOTH_RAWSUBJECT. Everything written uses the pseudonym P19.
% Output: paper/data/hgb_traces_P19.mat
% ------------------------------------------------------------------------------------------------
ECOG = getenv('SMOOTH_RAWDATA'); assert(~isempty(ECOG), 'Set SMOOTH_RAWDATA to the raw-data root (not distributed).');
RAWSUB = getenv('SMOOTH_RAWSUBJECT'); assert(~isempty(RAWSUB), 'Set SMOOTH_RAWSUBJECT to P19''s raw-data subject ID (not distributed).');
SUB  = 'P19'; CSUB = 18;                                    % pseudonym (demo/source.mat 19) for all outputs

addpath(fullfile(fileparts(mfilename('fullpath')),'..','..'));
PREPROC = getenv('SMOOTH_PREPROC'); assert(~isempty(PREPROC), 'Set SMOOTH_PREPROC to the preprocessing-code directory (not distributed).');
addpath(ECOG); addpath(PREPROC);
paths = smooth_setup('quiet');

% ---------------------------------------------------------------- preprocessing (proc_smooth.m)
fprintf('preprocessing %s ...\n', SUB);
data = preproc_tcg(['subject_' RAWSUB], 1, 1, 1);           % bipolar referencing, full trials
cfg = []; cfg.resamplefs = 500; data = ft_resampledata(cfg, data);
tmp = rmfield(data, {'trial','time'});
tmp.trial = {cat(2, data.trial{:})};
tmp.time  = {0:1/tmp.fsample:size(tmp.trial{1},2)/tmp.fsample - 1/tmp.fsample};
fprintf('  %d channels, %.1f s concatenated\n', numel(tmp.label), tmp.time{1}(end));

% ---------------------------------------------------------------- HGB (proc_smooth.m settings)
cfg = [];
cfg.output     = 'pow';   cfg.method    = 'mtmconvol';  cfg.taper = 'dpss';
cfg.foi        = 70:10:150;
cfg.tapsmofrq  = 10;
cfg.t_ftimwin  = .3*ones(numel(cfg.foi),1);
cfg.toi        = tmp.time{1}(1):.100:tmp.time{1}(end);
cfg.keeptrials = 'yes';   cfg.pad = 'nextpow2';
tfr = ft_freqanalysis(cfg, tmp);
for c = 1:numel(tfr.label)                                   % z-score per channel x frequency
    for f = 1:numel(tfr.freq)
        ts = tfr.powspctrm(:,c,f,:);
        tfr.powspctrm(:,c,f,:) = (ts - mean(ts(:),'omitnan')) / std(ts(:),'omitnan');
    end
end
cfg = []; cfg.avgoverfreq = 'yes'; cfg.nanmean = 'yes';
hg = ft_selectdata(cfg, tfr);
pow = squeeze(hg.powspctrm);                                 % nchan x ntime
toi = tfr.time(:).';
fprintf('  pow_hg %d x %d at %g Hz\n', size(pow,1), size(pow,2), 1/mean(diff(toi)));

% ---------------------------------------------------------------- verify against the stored HFBcorr
S = load(fullfile(paths.data,'corrsource.mat'));
ref = double(S.corrsource{CSUB}.HFBcorr); reflab = string(S.corrsource{CSUB}.label(:));
ok  = ~isnan(pow(1,:));                                      % proc_smooth drops these before corr
C   = corr(pow(:,ok)');
[~, ia] = ismember(reflab, string(tfr.label));
if any(ia == 0)
    error('export_hgb_traces:labels','%d stored labels absent from the recomputed data', nnz(ia==0));
end
err = max(abs(C(ia,ia) - ref), [], 'all');
fprintf('  HFBcorr agreement: max |recomputed - stored| = %.2e\n', err);
if err > 1e-6
    error('export_hgb_traces:mismatch', ...
        ['recomputed HFBcorr differs from corrsource{%d} by %.3e -- the pipeline has drifted, so ' ...
         'these traces are NOT the ones behind the berry map'], CSUB, err);
end

% ---------------------------------------------------------------- cache everything
% Reorder to the stored label order so Y, `chans`, and HFBcorr all share one indexing.
fs = 1/mean(diff(toi));
Y  = pow(ia,:);
fprintf('  cached %d channels x %d samples (%.0f s at %g Hz)\n\n', size(Y,1), size(Y,2), toi(end), fs);

out = struct('t', toi, 'Y', Y, 'chans', {cellstr(reflab)}, 'hfbcorr', ref, ...
             'fs', fs, 'subject', SUB, 'corrsource_idx', CSUB, 'hfbcorr_err', err, ...
             'recording_s', toi(end));
save(fullfile(paths.data,'hgb_traces_P19.mat'), '-struct', 'out', '-v7.3');
fprintf('wrote %s\n', fullfile(paths.data,'hgb_traces_P19.mat'));
end
