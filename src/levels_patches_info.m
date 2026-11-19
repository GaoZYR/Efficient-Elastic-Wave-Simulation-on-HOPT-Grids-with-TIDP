function [levels, patches, patch_catalog, patch_offset] = levels_patches_info(patch_velocity_dir)
    % Levels
    level_files = dir(fullfile(patch_velocity_dir,'level*'));
    level_index = nan(size(level_files));
    for k = 1:numel(level_files)
        tok = regexp(level_files(k).name, 'level(\d+)', 'tokens', 'once');
        if ~isempty(tok)
            level_index(k) = str2double(tok{1});
        end
    end
    [level_index, order] = sort(level_index);
    level_files = level_files(order);
    levels = length(level_index); % The number of levels
    
    % Patches
    for i = 1:levels
        level_fname = fullfile(level_files(i).folder, level_files(i).name);
        patch_files = dir(fullfile(level_fname,'patch_*.txt'));
        patch_index = nan(size(patch_files));
        for j = 1:numel(patch_files)
            tok = regexp(patch_files(j).name, 'patch_(\d+)\.txt', 'tokens', 'once');
            if ~isempty(tok)
                patch_index(j) = str2double(tok{1});
            end
        end
        [patch_index, order] = sort(patch_index);
        patch_files = patch_files(order);
        patches(i) = length(patch_index);
    end
    
    patch_catalog = build_patch_catalog(level_files, patches, levels);
    patch_offset = [0, cumsum(patches(1:end-1))];
end

function patch_catalog = build_patch_catalog(level_files, patches, levels)

    k = 0;
    patch_catalog = struct('level', {}, 'patch_id', {}, ...
                           'x_lo', {}, 'x_hi', {}, 'z_lo', {}, 'z_hi', {}, ...
                           'nx', {}, 'nz', {}, 'dx', {}, 'dz', {},...
                           'i_lo', {}, 'i_hi', {}, 'j_lo', {}, 'j_hi', {});

    for il = 1:levels
        level_fname = fullfile(level_files(il).folder, level_files(il).name);
        for ip = 1:patches(il)
            patch_fname = fullfile(level_fname, sprintf('patch_%d.txt', ip - 1));
            patch_par = read_patch_meta(patch_fname);

            k = k + 1;
            patch_catalog(k).level = il - 1;
            patch_catalog(k).patch_id = ip - 1;

            patch_catalog(k).x_lo = patch_par.x_lo;
            patch_catalog(k).x_hi = patch_par.x_hi;
            patch_catalog(k).z_lo = patch_par.z_lo;
            patch_catalog(k).z_hi = patch_par.z_hi;

            patch_catalog(k).nx = patch_par.nx;
            patch_catalog(k).nz = patch_par.nz;
            patch_catalog(k).dx = patch_par.dx;
            patch_catalog(k).dz = patch_par.dz;

            patch_catalog(k).i_lo = patch_par.i_lo;
            patch_catalog(k).i_hi = patch_par.i_hi;
            patch_catalog(k).j_lo = patch_par.j_lo;
            patch_catalog(k).j_hi = patch_par.j_hi;
        end
    end
end

function meta = read_patch_meta(fname)
% Read key-value text file written by patch velocity exporter.
% Expected keys:
% level, patch_id, i_lo, i_hi, j_lo, j_hi, nx, nz,
% x_lo, x_hi, z_lo, z_hi, dx, dz, ...

    fid = fopen(fname, 'r');
    assert(fid > 0, 'Failed to open %s', fname);

    meta = struct();
    while true
        tline = fgetl(fid);
        if ~ischar(tline); break; end
        tline = strtrim(tline);
        if isempty(tline) || startsWith(tline, '#')
            continue;
        end

        tok = regexp(tline, '^(\S+)\s+(.*)$', 'tokens', 'once');
        if isempty(tok)
            continue;
        end

        key = tok{1};
        val = tok{2};

        numv = str2double(val);
        if ~isnan(numv)
            meta.(key) = numv;
        else
            meta.(key) = val;
        end
    end
    fclose(fid);

    required = {'level','patch_id', ...
            'i_lo','i_hi','j_lo','j_hi', ...
            'nx','nz', ...
            'x_lo','x_hi','z_lo','z_hi', ...
            'dx','dz'};
    
    for k = 1:numel(required)
        assert(isfield(meta, required{k}), ...
            'Missing key "%s" in %s', required{k}, fname);
    end
end