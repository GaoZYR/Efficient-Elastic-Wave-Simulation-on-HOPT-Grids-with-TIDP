function [Vx_gather] = gather_wavefield_patch(Vx_patches, Vx_gather, patch_catalog, levels, patches, patch_offset, N)

    for il = 1:levels
        for ip = 1:patches(il)

            patch_global_index = ip + patch_offset(il);
            patch_par = patch_catalog(patch_global_index);
            nx = patch_par.nx; nz = patch_par.nz;
            ix_start_src = patch_par.i_lo + 1; ix_end_src = patch_par.i_hi + 1;
            iz_start_src = patch_par.j_lo + 1; iz_end_src = patch_par.j_hi + 1;
        
            ol = levels - il;
            dd = 2^ol;
            ix_start_dst =  dd*ix_start_src - (dd - 1);
            ix_end_dst   =  dd*ix_end_src   - (dd - 1);
            iz_start_dst =  dd*iz_start_src - (dd - 1);
            iz_end_dst   =  dd*iz_end_src   - (dd - 1);
        
            % Wavefield patches gather
            Vx = Vx_patches{il}{ip};
            Vx = Vx(N + 1:N + nz, N + 1:N + nx);
            for iz_shift = 1:dd
                rows = (iz_start_dst + iz_shift - 1):dd:(iz_end_dst + iz_shift - 1);
                for ix_shift = 1:dd
                    cols = (ix_start_dst + ix_shift - 1):dd:(ix_end_dst + ix_shift - 1);
                    Vx_gather(rows, cols) = Vx;
                end
            end
            
        end
    end