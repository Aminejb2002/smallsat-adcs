function list_spacecraft_blocks
% Lists the masked blocks in the Aerospace Blockset spacecraft libraries
% and writes results/spacecraft_blocks.txt.

root = fileparts(fileparts(mfilename('fullpath')));
outDir = fullfile(root, 'results');
if ~exist(outDir, 'dir'), mkdir(outDir); end
reportFile = fullfile(outDir, 'spacecraft_blocks.txt');
fid = fopen(reportFile, 'w');
closer = onCleanup(@() fclose(fid));

libs = {'aerolibcubesat', 'aerolibcubesatveh', 'aerolibsatdyn', ...
    'aerolibattitudesys', 'asbCubeSatBlockLib'};

for k = 1:numel(libs)
    emit(fid, '\n%s', libs{k});
    try
        load_system(libs{k});
        b = find_system(libs{k}, 'SearchDepth', 4, 'LookUnderMasks', 'all', ...
            'FollowLinks', 'on', 'Variants', 'AllVariants', ...
            'BlockType', 'SubSystem', 'Mask', 'on');
        for j = 1:numel(b)
            emit(fid, '  %s', strrep(b{j}, newline, ' '));
        end
        close_system(libs{k}, 0);
    catch err
        emit(fid, '  FAILED %s', err.message);
    end
end

emit(fid, '\nReport written to %s', reportFile);
end

function emit(fid, fmt, varargin)
msg = sprintf(fmt, varargin{:});
fprintf('%s\n', msg);
fprintf(fid, '%s\n', msg);
end
