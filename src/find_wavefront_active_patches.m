function [FMstart] = find_wavefront_active_patches(Tfir_dir, patch_catalog, levels, patches, patch_offset, src_x, src_z, dt, t_safe_offset)

    for il = 1:levels
        FMstart{il} = cell(patches(il), 1);

        Tfir_filename = fullfile(Tfir_dir, sprintf('T_level%d.mat', il - 1));
        T = load(Tfir_filename).T;

        for ip = 1:patches(il)
            patch = patch_catalog(patch_offset(il) + ip);
            ix_start = patch.i_lo + 1; ix_end = patch.i_hi + 1;
            iz_start = patch.j_lo + 1; iz_end = patch.j_hi + 1;
            T0 = T(iz_start:iz_end, ix_start:ix_end);
            t_start = min(T0(:));

            % add safe time offset
            t_start = t_start - t_safe_offset;
            if t_start < 0
                t_start = 0;
            end

            % repair source
            if src_x >= patch.x_lo && src_x <= patch.x_hi && src_z >= patch.z_lo && src_z <= patch.z_hi
                t_start = 0;
            end
            
            % save min start time
            it_start = floor(t_start/dt);
            FMstart{il}{ip} = it_start + 1;
        end

    end
end
