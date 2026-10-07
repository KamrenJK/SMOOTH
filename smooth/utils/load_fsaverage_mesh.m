function cortex = load_fsaverage_mesh(fshome)
% load cortex mesh given fs path

lh = ft_read_headshape(fullfile(fshome,'subjects','fsaverage','surf','lh.pial'));
rh = ft_read_headshape(fullfile(fshome,'subjects','fsaverage','surf','rh.pial'));

cortex      = [];
cortex.pos  = [lh.pos; rh.pos];
cortex.tri  = [lh.tri; rh.tri + size(lh.pos,1)];
cortex.unit = lh.unit;

if isfield(lh,'curv') && isfield(rh,'curv')
    cortex.curv = [lh.curv; rh.curv];
else
    cortex.curv = zeros(size(cortex.pos,1),1);
end
end