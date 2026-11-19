%% This function is used to find the patches that source in
function info = locate_point_in_highest_patch(patch_catalog, x, z)

    info = struct( ...
        'level', [], 'patch_id', [], ...
        'x_lo', [], 'x_hi', [], 'z_lo', [], 'z_hi', [], ...
        'nx', [], 'nz', [], 'dx', [], 'dz', [], ...
        'i_lo', [], 'i_hi', [], 'j_lo', [], 'j_hi', []);

    best_level = -inf;

    for k = 1:numel(patch_catalog)
        pk = patch_catalog(k);

        if x >= pk.x_lo && x <= pk.x_hi && z >= pk.z_lo && z <= pk.z_hi
            if pk.level > best_level
                best_level = pk.level;
                info = pk;
            end
        end
    end
end