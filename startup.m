root = fileparts(mfilename('fullpath'));
addpath(root);
dirs = {'models', 'src', 'analysis', 'experiments', 'tests', 'tools', 'visualization'};
for k = 1:numel(dirs)
    addpath(genpath(fullfile(root, dirs{k})));
end
clear root dirs k
