function test_smooth
% TEST_SMOOTH - functional test suite for the SMOOTH toolbox
%
% Run from anywhere:
%   test_smooth
%
% Exercises every public function with small permutation counts, checks
% reproducibility under a fixed seed, and confirms that input validation
% rejects malformed data. Takes a couple of minutes.

addpath(fullfile(fileparts(mfilename('fullpath')),'..'));
paths = smooth_setup('quiet');
S = load(fullfile(paths.demo,'source.mat')); source = S.source;
FS = paths.freesurfer;

pass = 0; fail = 0;
function check(name, fn)
    try
        fn();
        fprintf('  PASS  %s\n', name); pass = pass + 1;
    catch ME
        fprintf('  FAIL  %s\n          %s\n', name, ME.message); fail = fail + 1;
    end
end

fprintf('\nSMOOTH test suite\n=================\n');

%% --- SMOOTHstat
fprintf('\nSMOOTHstat\n');

check('runs and returns documented fields', @() assertFields( ...
    runstat(FS, source, 5, 1), {'stat','prob','mask','posclusters','negclusters','coverage','cfg'}));

check('is reproducible under a fixed seed', @() assertEqualRuns( ...
    runstat(FS, source, 5, 42), runstat(FS, source, 5, 42)));

check('differs under a different seed', @() assertDifferentRuns( ...
    runstat(FS, source, 5, 42), runstat(FS, source, 5, 7)));

check('records the seed in stat.cfg', @() assert( ...
    isequal(getfield(runstat(FS, source, 3, 99),'cfg').randomseed, 99)));

check('cfg.tail does not change the t-map', @() assertSameStat( ...
    runstat(FS, source, 5, 3, 'both'), runstat(FS, source, 5, 3, 'positive')));

%% --- validation
fprintf('\ninput validation\n');

check('rejects a missing parameter field', @() assertErrors( ...
    @() SMOOTHstat(basecfg(FS), rmfield(source{1},'stat')), 'SMOOTHstat:missingParameter'));

check('rejects a size mismatch', @() assertErrors( ...
    @() SMOOTHstat(basecfg(FS), setfield(source{1},'stat',source{1}.stat(1:5))), ...
    'SMOOTHstat:sizeMismatch'));

check('rejects NaN in the statistic', @() assertErrors( ...
    @() SMOOTHstat(basecfg(FS), setfield(source{1},'stat',[NaN; source{1}.stat(2:end)])), ...
    'SMOOTHstat:nonFiniteData'));

check('rejects missing chanpos', @() assertErrors( ...
    @() SMOOTHstat(basecfg(FS), rmfield(source{1},'elec')), 'SMOOTHstat:missingChanpos'));

check('rejects an empty call', @() assertErrors( ...
    @() SMOOTHstat(basecfg(FS)), 'SMOOTHstat:noData'));

%% --- SMOOTHsub
fprintf('\nSMOOTHsub\n');

check('runs with default cfg', @() assertFields( ...
    runsub(FS, source{10}, 3, false), {'data','cfg'}));

check('returns sflipmat of the right size', @() assert( ...
    size(getfield(runsub(FS, source{10}, 4, false),'data').sflipmat, 2) == 4));

check('returns surrogates when requested', @() assertFields( ...
    runsub(FS, source{10}, 3, true), {'surrogate_subject','maps'}));

check('does not return a phantom surrogate_group', @() assert( ...
    ~isfield(runsub(FS, source{10}, 3, true), 'surrogate_group')));

%% --- SMOOTHdummy
fprintf('\nSMOOTHdummy\n');

check('runs (RSF) and returns a t-map', @() assertFields( ...
    rundummy(FS, source, 5), {'stat','mask','posclusters','cfg'}));

%% --- backward compatibility
fprintf('\ncompatibility\n');

check('honours deprecated cfg.keepsurrogates', @() assertFields( ...
    runlegacy(FS, source), {'surrogate_subject','surrogate_group'}));

%% --- cross-file consistency
% SMOOTHstat, SMOOTHsub and SMOOTHdummy are each self-contained and carry their
% OWN private copy of mm2sphere and of the kernel defaults. If those drift apart,
% the same cfg silently produces different smoothing in different functions --
% which is exactly what happened once during development.
fprintf('\ncross-file consistency\n');

check('mm2sphere factor identical in all three functions', @() assertSameFactor());
check('kernelwidth/graphsigma defaults identical', @() assertSameDefaults());

fprintf('\n=================\n%d passed, %d failed\n\n', pass, fail);
if fail > 0, error('test_smooth:failures','%d test(s) failed', fail); end
end

% ---------------- helpers ----------------
function cfg = basecfg(FS)
cfg = []; cfg.fshome = FS; cfg.numrandomization = 3;
cfg.randomseed = 1; cfg.feedback = 'no';
end

function s = runstat(FS, source, n, seed, tail)
cfg = basecfg(FS); cfg.numrandomization = n; cfg.randomseed = seed;
if nargin > 4, cfg.tail = tail; end
w = warning('off','all'); s = SMOOTHstat(cfg, source{:}); warning(w);
end

function s = runsub(FS, subj, n, keep)
cfg = basecfg(FS); cfg.numrandomization = n;
if keep, cfg.keepsurrogates = 'yes'; cfg.keepmaps = 'yes'; end
w = warning('off','all'); s = SMOOTHsub(cfg, subj); warning(w);
end

function s = rundummy(FS, source, n)
cfg = basecfg(FS); cfg.numrandomization = n;
w = warning('off','all'); s = SMOOTHdummy(cfg, source{:}); warning(w);
end

function s = runlegacy(FS, source)
cfg = basecfg(FS); cfg.keepsurrogates = 'yes';
cfg.subsamplesurr = 'yes'; cfg.nsurrsamples = 2;
w = warning('off','all'); s = SMOOTHstat(cfg, source{:}); warning(w);
end

function assertFields(s, flds)
for k = 1:numel(flds)
    assert(isfield(s, flds{k}), 'missing field: %s', flds{k});
end
end

function assertEqualRuns(a, b)
assert(isequaln(a.stat, b.stat), 't-maps differ');
assert(isequaln(a.prob, b.prob), 'p-values differ');
assert(isequaln(a.mask, b.mask), 'masks differ');
end

function assertDifferentRuns(a, b)
assert(~isequaln(a.prob, b.prob), 'different seeds produced identical results');
end

function assertSameStat(a, b)
assert(isequaln(a.stat, b.stat), 'cfg.tail altered the t-map');
end

function assertErrors(fn, id)
try
    fn();
    error('assert:noError', 'expected error %s, none raised', id);
catch ME
    assert(strcmp(ME.identifier, id), 'expected %s, got %s', id, ME.identifier);
end
end

function assertSameFactor()
% Each function carries a private mm2sphere; they must agree.
files = {'SMOOTHstat','SMOOTHsub','SMOOTHdummy'};
vals  = nan(size(files));
for k = 1:numel(files)
    src = fileread(which(files{k}));
    tok = regexp(src, 'function\s+su\s*=\s*mm2sphere.*?su\s*=\s*([0-9.]+)\s*\*\s*mm', 'tokens', 'once');
    assert(~isempty(tok), 'could not locate mm2sphere factor in %s', files{k});
    vals(k) = str2double(tok{1});
end
assert(all(vals == vals(1)), ...
    'mm2sphere factors differ: %s = %s', strjoin(files,'/'), mat2str(vals));
end

function assertSameDefaults()
% kernelwidth and graphsigma defaults must agree across the three functions.
files = {'SMOOTHstat','SMOOTHsub','SMOOTHdummy'};
for opt = {'kernelwidth','graphsigma'}
    o = opt{1}; vals = nan(size(files));
    for k = 1:numel(files)
        src = fileread(which(files{k}));
        pat = sprintf('ft_getopt\\(cfg,\\s*''%s''\\s*,\\s*([0-9.]+)\\s*\\)', o);
        tok = regexp(src, pat, 'tokens', 'once');
        assert(~isempty(tok), 'could not locate %s default in %s', o, files{k});
        vals(k) = str2double(tok{1});
    end
    assert(all(vals == vals(1)), ...
        '%s defaults differ across %s: %s', o, strjoin(files,'/'), mat2str(vals));
end
end
