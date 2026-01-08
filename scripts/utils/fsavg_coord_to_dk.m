% -------------------------------------------------------------------------
function [roi_name, roi_index, hemi, vert_idx] = fsavg_coord_to_dk(pos_mm, fshome)
% Map fsaverage RAS coords (mm) to Desikan–Killiany (aparc) ROI names.
%
% Inputs
%   pos_mm  : N x 3 coordinates in fsaverage space (RAS, millimeters)
%   fshome  : path to FreeSurfer home (string), e.g., '/usr/local/freesurfer'
%
% Outputs
%   roi_name : N x 1 cell array of DK ROI names (now prefixed with 'lh-' or 'rh-')  %% NEW
%   roi_index: N x 1 double; index into the colortable.struct_names for the chosen hemi
%              (NaN if unknown)
%   hemi     : N x 1 char; 'L' or 'R' for the hemisphere used
%   vert_idx : N x 1 int32; vertex index on the chosen pial surface
%
% Requirements
%   - FieldTrip (for ft_read_headshape) or you can replace with freesurfer_read_surf, etc.
%   - FreeSurfer MATLAB tools on the path (for read_annotation.m):
%       addpath(fullfile(getenv('FREESURFER_HOME'),'matlab'))
%
% Notes
%   - Uses nearest-vertex on each hemi pial, picks the closer hemi.
%   - DK atlas files used:
%       $SUBJECTS_DIR/fsaverage/label/lh.aparc.annot
%       $SUBJECTS_DIR/fsaverage/label/rh.aparc.annot
%
% Example
%   [roi, idx, h, v] = fsavg_coord_to_dk(pos, '/usr/local/freesurfer');

  arguments
    pos_mm (:,3) double
    fshome (1,1) string
  end

  % Resolve SUBJECTS_DIR
  subj_dir = fullfile(fshome, 'subjects');
  fsavg    = fullfile(subj_dir, 'fsaverage');

  % Check for read_annotation
  if ~exist('read_annotation','file')
    addpath(fullfile(char(fshome), 'matlab')); % try to add FS matlab
    if ~exist('read_annotation','file')
      error(['read_annotation.m not found. Add FreeSurfer MATLAB tools to the path, e.g., ' ...
             'addpath(fullfile(getenv(''FREESURFER_HOME''),''matlab''))']);
    end
  end

  % Persistent caches for speed
  persistent L_pial R_pial kdtL kdtR L_annot L_ctab R_annot R_ctab
  if isempty(L_pial) || isempty(R_pial)
    % Load pial meshes (mm)
    L_pial = ft_read_headshape(fullfile(fsavg,'surf','lh.pial'));
    R_pial = ft_read_headshape(fullfile(fsavg,'surf','rh.pial'));

    % KD-trees for nearest-vertex search
    kdtL   = createns(L_pial.pos, 'NSMethod','kdtree');
    kdtR   = createns(R_pial.pos, 'NSMethod','kdtree');

    % Read DK annotations (aparc) for each hemi
    [~, L_annot, L_ctab] = read_annotation(fullfile(fsavg,'label','lh.aparc.annot'));
    [~, R_annot, R_ctab] = read_annotation(fullfile(fsavg,'label','rh.aparc.annot'));
  end

  N = size(pos_mm,1);
  roi_name  = repmat({''}, N, 1);
  roi_index = nan(N,1);
  hemi      = repmat('L', N, 1);
  vert_idx  = nan(N,1,'double');

  % Nearest vertex in each hemi
  [idxL, dL] = knnsearch(kdtL, pos_mm);  % left hemi
  [idxR, dR] = knnsearch(kdtR, pos_mm);  % right hemi

  useLeft  = dL <= dR;
  useRight = ~useLeft;

  % Assign hemisphere & vertex
  hemi(useLeft)     = 'L';
  vert_idx(useLeft) = int32(idxL(useLeft));
  hemi(useRight)    = 'R';
  vert_idx(useRight)= int32(idxR(useRight));

  % Convert vertex -> label id -> ROI name via colortable
  % FreeSurfer .annot label vector is color codes; match against ctab.table(:,5)
  % (RGBA packed int).

  % Left hemi
  if any(useLeft)
    labL = double(L_annot(vert_idx(useLeft)));  % color code at vertex
    tL   = double(L_ctab.table(:,5));
    % Map color code to index into struct_names
    [~, mapIdxL] = ismember(labL, tL);
    roi_index(useLeft) = mapIdxL;
    namesL = L_ctab.struct_names;
    for k = find(useLeft).'
      mi = roi_index(k);
      if mi>=1 && mi<=numel(namesL)
        roi_name{k} = ['lh-' namesL{mi}];    
      else
        roi_name{k} = 'lh-unknown';            
        roi_index(k)= NaN;
      end
    end
  end

  % Right hemi
  if any(useRight)
    labR = double(R_annot(vert_idx(useRight))); % color code at vertex
    tR   = double(R_ctab.table(:,5));
    [~, mapIdxR] = ismember(labR, tR);
    roi_index(useRight) = mapIdxR;
    namesR = R_ctab.struct_names;
    for k = find(useRight).'
      mi = roi_index(k);
      if mi>=1 && mi<=numel(namesR)
        roi_name{k} = ['rh-' namesR{mi}];         
      else
        roi_name{k} = 'rh-unknown';         
        roi_index(k)= NaN;
      end
    end
  end
end