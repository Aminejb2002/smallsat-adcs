function probe_blocks
% Checks which Aerospace Blockset blocks and Aerospace Toolbox functions exist
% in this MATLAB install and writes results/probe_report.txt.

root = fileparts(fileparts(mfilename('fullpath')));
outDir = fullfile(root, 'results');
if ~exist(outDir, 'dir'), mkdir(outDir); end
reportFile = fullfile(outDir, 'probe_report.txt');
fid = fopen(reportFile, 'w');
closer = onCleanup(@() fclose(fid));

emit(fid, 'MATLAB %s', version('-release'));

blocks = {
    'Quaternion Multiplication'
    'Quaternion Normalize'
    'Quaternion Conjugate'
    'Quaternion Rotation'
    'Quaternions to Direction Cosine Matrix'
    'Direction Cosine Matrix to Quaternions'
    'Orbit Propagator'
    'Zonal Harmonic Gravity Model'
    'Spherical Harmonic Gravity Model'
    'World Magnetic Model'
    'COESA Atmosphere Model'
    'NRLMSISE-00 Atmosphere Model'
    'Planetary Ephemeris'
    'Julian Date Conversion'
    '6DOF (Quaternion)'
    };

emit(fid, '\nAerospace folders under toolbox:');
d = dir(fullfile(matlabroot, 'toolbox', 'aero*'));
for k = 1:numel(d)
    emit(fid, '  %s', fullfile(d(k).folder, d(k).name));
end
if isempty(d), emit(fid, '  none'); end

emit(fid, '\nLibrary files found:');
libs = findLibraries();
if isempty(libs), emit(fid, '  none'); end
loaded = {};
for k = 1:numel(libs)
    try
        load_system(libs(k).path);
        loaded{end+1} = libs(k).name; %#ok<AGROW>
        emit(fid, '  loaded   %s', libs(k).path);
    catch err
        emit(fid, '  FAILED   %s (%s)', libs(k).path, err.message);
    end
end

emit(fid, '\nBlocks (searched by name in the loaded libraries):');
for k = 1:numel(blocks)
    hit = '';
    for j = 1:numel(loaded)
        try
            h = find_system(loaded{j}, 'LookUnderMasks', 'all', 'FollowLinks', 'on', ...
                'Variants', 'AllVariants', 'Name', blocks{k});
        catch
            h = {};
        end
        if ~isempty(h)
            hit = h{1};
            break
        end
    end
    if isempty(hit)
        emit(fid, '  MISSING  %s', blocks{k});
    else
        emit(fid, '  found    %s -> %s', blocks{k}, hit);
    end
end
for j = 1:numel(loaded)
    try close_system(loaded{j}, 0); catch, end %#ok<NOCOMMA>
end

funcs = {'igrfmagm', 'wrldmagm', 'atmosnrlmsise00', 'atmoscoesa', 'planetEphemeris', ...
    'juliandate', 'quatmultiply', 'quatnormalize', 'quat2dcm', 'dcm2quat', 'quat2eul', ...
    'eul2quat', 'eci2ecef', 'satelliteScenario'};

emit(fid, '\nFunctions:');
for k = 1:numel(funcs)
    if exist(funcs{k}, 'file') || exist(funcs{k}, 'builtin')
        emit(fid, '  found    %s', funcs{k});
    else
        emit(fid, '  MISSING  %s', funcs{k});
    end
end

emit(fid, '\nProgrammatic build check:');
mdl = 'probe_tmp_model';
try
    new_system(mdl);
    add_block('simulink/User-Defined Functions/MATLAB Function', [mdl '/fcn']);
    add_block('simulink/Continuous/Integrator', [mdl '/int']);
    emit(fid, '  ok       MATLAB Function and Integrator blocks added');
catch err
    emit(fid, '  FAILED   %s', err.message);
end
if bdIsLoaded(mdl), close_system(mdl, 0); end

emit(fid, '\nReport written to %s', reportFile);
end

function libs = findLibraries()
libs = struct('name', {}, 'path', {});
base = fullfile(matlabroot, 'toolbox', 'aeroblks');
if ~exist(base, 'dir'), return; end
patterns = {'*.slx', '*.mdl'};
for p = 1:numel(patterns)
    files = dir(fullfile(base, '**', patterns{p}));
    for k = 1:numel(files)
        [~, n] = fileparts(files(k).name);
        if contains(lower(n), 'lib')
            libs(end+1) = struct('name', n, 'path', fullfile(files(k).folder, files(k).name)); %#ok<AGROW>
        end
    end
end
end

function emit(fid, fmt, varargin)
msg = sprintf(fmt, varargin{:});
fprintf('%s\n', msg);
fprintf(fid, '%s\n', msg);
end
