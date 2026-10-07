function paths = smooth_setup(varargin)
% SMOOTH_SETUP - one-time path configuration for the SMOOTH toolbox
%
% EDIT THE TWO PATHS BELOW, then everything else (demo, paper figures, tests)
% works without further changes.
%
% Use as:
%   paths = smooth_setup();          % configure paths and add them to the MATLAB path
%   paths = smooth_setup('quiet');   % same, without the banner
%
% Returns a struct with:
%   paths.root        SMOOTH repository root
%   paths.core        core functions (smooth/)
%   paths.demo        demo dataset directory
%   paths.data        paper data directory
%   paths.figs        paper figure output directory (created on demand)
%   paths.fieldtrip   FieldTrip root
%   paths.freesurfer  FreeSurfer home (must contain subjects/fsaverage)
%
% Requirements:
%   MATLAB    R2022a or newer (tested on R2024b), Statistics and Machine Learning Toolbox
%   FieldTrip 20220000 or newer  -- https://www.fieldtriptoolbox.org
%   FreeSurfer 7.x with fsaverage and fsaverage5 subjects installed

% =========================================================================
%                    >>>  EDIT THESE TWO LINES  <<<
% =========================================================================
paths.fieldtrip  = '/path/to/fieldtrip';
paths.freesurfer = '/Applications/freesurfer/7.3.2';
% =========================================================================
% (Alternatively, set the FIELDTRIP_HOME / FREESURFER_HOME environment variables,
% which take precedence over the two lines above.)
if ~isempty(getenv('FIELDTRIP_HOME')),  paths.fieldtrip  = getenv('FIELDTRIP_HOME');  end
if ~isempty(getenv('FREESURFER_HOME')), paths.freesurfer = getenv('FREESURFER_HOME'); end

quiet = any(strcmpi(varargin, 'quiet'));

% repository layout is derived from this file's own location
paths.root   = fileparts(mfilename('fullpath'));
paths.core   = fullfile(paths.root, 'smooth');
paths.demo   = fullfile(paths.root, 'demo');
paths.data   = fullfile(paths.root, 'paper', 'data');
paths.figs   = fullfile(paths.root, 'paper', 'figs');

% ---- validate before doing anything else, with actionable messages
if ~isfolder(paths.fieldtrip)
  error('smooth_setup:noFieldTrip', ...
    ['FieldTrip not found at:\n    %s\n' ...
     'Edit paths.fieldtrip in %s.'], paths.fieldtrip, [mfilename '.m']);
end
if ~isfolder(paths.freesurfer)
  error('smooth_setup:noFreeSurfer', ...
    ['FreeSurfer not found at:\n    %s\n' ...
     'Edit paths.freesurfer in %s.'], paths.freesurfer, [mfilename '.m']);
end
fsavg = fullfile(paths.freesurfer, 'subjects', 'fsaverage', 'surf', 'lh.pial');
if ~isfile(fsavg)
  error('smooth_setup:noFsaverage', ...
    ['FreeSurfer found, but fsaverage surfaces are missing:\n    %s\n' ...
     'SMOOTH requires the fsaverage subject (and fsaverage5 for the Figure 3-4 ' ...
     'simulations).'], fsavg);
end

% ---- add to path
addpath(genpath(paths.core));
addpath(genpath(fullfile(paths.root, 'external')));
addpath(paths.fieldtrip);
evalc('ft_defaults');   % suppress FieldTrip's banner

if ~quiet
  fprintf('SMOOTH configured\n');
  fprintf('  root       : %s\n', paths.root);
  fprintf('  FieldTrip  : %s\n', paths.fieldtrip);
  fprintf('  FreeSurfer : %s\n', paths.freesurfer);
end
end
