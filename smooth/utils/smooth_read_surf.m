function shape = smooth_read_surf(filename)
% SMOOTH_READ_SURF - read a FreeSurfer triangle surface (e.g. lh.pial, lh.sphere.reg)
%
% Use as:
%   shape = smooth_read_surf(filename)
%
% Returns the struct that ft_read_headshape returns for a FreeSurfer surface:
%   shape.pos    vertex coordinates (FreeSurfer surface RAS, mm), Nx3 double
%   shape.tri    triangles (1-based), Mx3 double
%   shape.unit   'mm'
%   shape.sulc, .curv, .area, .thickness   per-vertex values, when the matching
%                files (e.g. lh.sulc next to lh.pial) exist
%
% SMOOTH reads every surface through this self-contained parser instead of
% ft_read_headshape, because ft_read_headshape calls whichever read_surf.m comes
% first on the MATLAB path, and several toolboxes ship their own copy. At least
% one widely installed copy (iELVis, FSURF_MATLAB/read_surf.m) adds 1 to every
% vertex coordinate, which silently moves the cortex by +[1 1 1] mm relative to
% the electrodes. Reading the files directly makes results independent of the path.
%
% Unlike ft_read_headshape, inflated surfaces are returned as stored (FieldTrip
% moves the two inflated hemispheres apart for display).
%
% The file formats are read as in FreeSurfer's read_surf.m and read_curv.m (Bruce
% Fischl, MGH; FreeSurfer Software License), including the header fallback that
% FieldTrip's copy adds.

TRIANGLE_FILE_MAGIC_NUMBER = 16777214;

filename = char(filename);          % accept string scalars as well as char vectors
fid = fopen(filename, 'rb', 'b');   % FreeSurfer surfaces are big-endian
if fid < 0
  error('smooth_read_surf:open', 'Could not open surface file %s.', filename);
end
closer = onCleanup(@() fclose(fid));

magic = fread3(fid);
if magic ~= TRIANGLE_FILE_MAGIC_NUMBER
  error('smooth_read_surf:format', ...
    '%s is not a FreeSurfer triangle surface (magic number %d).', filename, magic);
end

% "created by <user> on <date>\n\n", then the vertex and face counts
fgets(fid); fgets(fid);
vnum = fread(fid, 1, 'int32');
fnum = fread(fid, 1, 'int32');
if isempty(vnum) || isempty(fnum) || round(abs(vnum)) ~= vnum || round(abs(fnum)) ~= fnum
  % comment line not followed by a blank line (same fallback as FieldTrip's read_surf)
  frewind(fid); fread3(fid); fgets(fid);
  vnum = fread(fid, 1, 'int32');
  fnum = fread(fid, 1, 'int32');
end
if isempty(vnum) || isempty(fnum) || vnum < 0 || fnum < 0
  error('smooth_read_surf:truncated', '%s has no valid vertex/face counts.', filename);
end

pos = fread(fid, vnum*3, 'float32');
tri = fread(fid, fnum*3, 'int32');
if numel(pos) ~= vnum*3 || numel(tri) ~= fnum*3
  error('smooth_read_surf:truncated', '%s is truncated.', filename);
end

shape      = struct();
shape.pos  = reshape(pos, 3, vnum)';
shape.tri  = reshape(tri, 3, fnum)';
if min(shape.tri(:)) == 0
  shape.tri = shape.tri + 1;        % FreeSurfer counts vertices from 0
end
shape.unit = 'mm';

% per-vertex morphometry stored next to the surface, as ft_read_headshape does
[p, name] = fileparts(filename);
morph = {'sulc', 'curv', 'area', 'thickness'};
for k = 1:numel(morph)
  f = fullfile(p, [name '.' morph{k}]);
  if exist(f, 'file')
    shape.(morph{k}) = read_morph(f);
  end
end

% -------------------------------------------------------------------------
function vals = read_morph(filename)
% FreeSurfer curvature-format file (new format; old int16 format as fallback)
fid = fopen(filename, 'rb', 'b');
if fid < 0
  error('smooth_read_surf:open', 'Could not open file %s.', filename);
end
closer = onCleanup(@() fclose(fid));
vnum = fread3(fid);
if vnum == 16777215
  vnum = fread(fid, 1, 'int32');
  fread(fid, 1, 'int32');           % number of faces
  fread(fid, 1, 'int32');           % values per vertex
  vals = fread(fid, vnum, 'float');
else
  fread3(fid);                      % number of faces
  vals = fread(fid, vnum, 'int16') ./ 100;
end

% -------------------------------------------------------------------------
function n = fread3(fid)
% 3-byte big-endian integer
b = fread(fid, 3, 'uchar');
if numel(b) < 3
  error('smooth_read_surf:truncated', 'Unexpected end of file in %s.', fopen(fid));
end
n = bitshift(b(1), 16) + bitshift(b(2), 8) + b(3);
