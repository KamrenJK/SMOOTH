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

%% --- surface geometry
% A read_surf.m that adds 1 mm to every vertex (as the copy in iELVis does) silently moves the
% cortex relative to the electrodes when it precedes FieldTrip's copy on the path.
fprintf('\nsurface geometry\n');

check('surfaces do not depend on read_surf.m on the path', @() assertPathIndependent(FS));

check('SMOOTHstat results do not depend on read_surf.m on the path', @() assertStatPathIndependent(FS, source(1:5)));

check('warns when electrodes sit on a +1 mm shifted surface', @() assertWarning( ...
    @() runquiet(FS, shift(source(1:5), 1)), 'SMOOTHstat:shiftedFrame', true));

check('does not warn for the demo electrodes', @() assertWarning( ...
    @() runquiet(FS, source(1:5)), 'SMOOTHstat:shiftedFrame', false));

check('does not warn for electrodes 1.5 mm off the surface', @() assertWarning( ...
    @() runquiet(FS, offsurface(FS, source(1:10), 1.5)), 'SMOOTHstat:shiftedFrame', false));

%% --- cross-file consistency
% SMOOTHstat, SMOOTHsub and SMOOTHdummy are each self-contained and carry their
% OWN private copy of mm2sphere and of the kernel defaults. If those drift apart,
% the same cfg silently produces different smoothing in different functions.
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

function runquiet(FS, source)
cfg = basecfg(FS); cfg.numrandomization = 2;
SMOOTHstat(cfg, source{:});
end

function src = shift(src, d)
for k = 1:numel(src), src{k}.elec.chanpos = src{k}.elec.chanpos + d; end
end

function src = offsurface(FS, src, c)
% move every electrode c mm along the outward normal of its nearest fsaverage pial vertex
for h = {'lh','rh'}
    s = smooth_read_surf(fullfile(FS,'subjects','fsaverage','surf',[h{1} '.pial']));
    fn = cross(s.pos(s.tri(:,2),:) - s.pos(s.tri(:,1),:), s.pos(s.tri(:,3),:) - s.pos(s.tri(:,1),:), 2);
    n = zeros(size(s.pos));                      % vertex normal = sum of adjacent face normals
    for k = 1:3
        for d = 1:3, n(:,d) = n(:,d) + accumarray(s.tri(:,k), fn(:,d), [size(s.pos,1) 1]); end
    end
    n = n ./ vecnorm(n, 2, 2);
    for i = 1:numel(src)
        p = src{i}.elec.chanpos; sel = (p(:,1) < 0) == strcmp(h{1},'lh');
        if any(sel), src{i}.elec.chanpos(sel,:) = p(sel,:) + c*n(knnsearch(s.pos, p(sel,:)),:); end
    end
end
end

function assertPathIndependent(FS)
% Put a read_surf.m that shifts every vertex by +1 mm first on the path: ft_read_headshape
% follows it, smooth_read_surf must not.
f    = fullfile(FS,'subjects','fsaverage5','surf','lh.pial');
ref  = smooth_read_surf(f);
shim = makeShim(FS);
cleanup = onCleanup(@() removeShim(shim));
b = smooth_read_surf(f);
assert(isequal(b, ref), 'smooth_read_surf changed with a different read_surf.m on the path');
end

function assertStatPathIndependent(FS, source)
cfg = basecfg(FS); cfg.numrandomization = 3;
w = warning('off','all'); restore = onCleanup(@() warning(w));
ref  = SMOOTHstat(cfg, source{:});
shim = makeShim(FS);
cleanup = onCleanup(@() removeShim(shim));
got  = SMOOTHstat(cfg, source{:});
assert(isequaln(ref.stat, got.stat) && isequaln(ref.prob, got.prob) && isequaln(ref.coverage, got.coverage), ...
    'SMOOTHstat results changed with a shifted read_surf.m on the path');
end

function shim = makeShim(FS)
% a read_surf.m that adds 1 mm to every vertex, first on the path (as iELVis's copy can be)
ft_hastoolbox('freesurfer', 1);   % add FieldTrip's copy now, so it cannot land in front of the shim later
shim = tempname; mkdir(shim);
fid  = fopen(fullfile(shim,'read_surf.m'), 'w');
fprintf(fid, 'function [v, f] = read_surf(fname)\ns = smooth_read_surf(fname); v = s.pos + 1; f = s.tri - 1;\n');
fclose(fid);
addpath(shim, '-begin'); rehash;
f = fullfile(FS,'subjects','fsaverage5','surf','lh.pial');
a = ft_read_headshape(f); b = smooth_read_surf(f);
if ~isequal(a.pos, b.pos + 1)
    removeShim(shim);
    error('test:setup', 'test setup: the shifted read_surf.m was not picked up');
end
end

function removeShim(shim)
rmpath(shim); rmdir(shim, 's'); rehash;
end

function assertWarning(fn, id, expected)
% expected = true: fn must issue warning id; false: it must not. The warning is turned into
% an error for the duration, so no other warning can mask it.
w = warning('off','all'); warning('error', id);
restore = onCleanup(@() warning(w));
got = '';
try
    evalc('fn()');
catch ME
    got = ME.identifier;
    if ~strcmp(got, id), rethrow(ME); end
end
if expected
    assert(strcmp(got, id), 'expected warning %s, none issued', id);
else
    assert(isempty(got), 'unexpected warning %s', id);
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
