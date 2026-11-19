function [material_patches, vmax] = ...
    read_model_parameters(level_vp_dir, level_vs_dir, level_rho_dir, ...
    levels, patches, patch_catalog, patch_offset, Nx, Nz)

    for il = 1:levels
        level_name = sprintf('level%d', il - 1);
        vs_level_fname = fullfile(level_vs_dir, sprintf('level%d_vs.bin', il - 1));
        vp_level_fname = fullfile(level_vp_dir, sprintf('level%d_vp.bin', il - 1));
        rho_level_fname = fullfile(level_rho_dir, sprintf('level%d_rho.bin', il - 1));

        %%%%%%%%%%%%%%%%%%% This part should be modified %%%%%%%%%%%%%%%%%%%
        fid = fopen(vs_level_fname, 'rb');
            vs_levels.(level_name) = fread(fid, Nx(il) * Nz(il), 'float32=>double');
        fclose(fid);
        vs_levels.(level_name) = reshape(vs_levels.(level_name), [Nx(il), Nz(il)]).';

        fid = fopen(vp_level_fname, 'rb');
            vp_levels.(level_name) = fread(fid, Nx(il) * Nz(il), 'float32=>double');
        fclose(fid);
        vp_levels.(level_name) = reshape(vp_levels.(level_name), [Nx(il), Nz(il)]).';

        fid = fopen(rho_level_fname, 'rb');
            rho_levels.(level_name) = fread(fid, Nx(il) * Nz(il), 'float32=>double');
        fclose(fid);
        rho_levels.(level_name) = reshape(rho_levels.(level_name), [Nx(il), Nz(il)]).';

        % vp_levels.(level_name) = vs_levels.(level_name).*1.5;
        % rho_levels.(level_name) = ones(Nz(il), Nx(il));
        %%%%%%%%%%%%%%%%%%% This part should be modified %%%%%%%%%%%%%%%%%%%
    
        vp = vp_levels.(level_name);
        vs = vs_levels.(level_name);
        rho = rho_levels.(level_name);
    
        % Full grid
        lamda_levels.(level_name) = rho.*(vp.^2 - 2.*vs.^2);
        niu_levels.(level_name) = rho.*(vs.^2);
    
        % Half grid
        lamda = lamda_levels.(level_name);
        niu = niu_levels.(level_name);
        for ix = 1:Nx(il)
            for iz = 1:Nz(il)
                if ix ~= Nx(il) && iz ~= Nz(il)
                    lamda_half_x_levels.(level_name)(iz, ix) = 0.5*(lamda(iz, ix + 1) + lamda(iz, ix));
                    niu_half_x_levels.(level_name)(iz, ix) = 0.5*(niu(iz, ix + 1) + niu(iz, ix));
                    niu_half_z_levels.(level_name)(iz, ix) = 0.5*(niu(iz + 1, ix) + niu(iz, ix));
                    rho_half_x_half_z_levels.(level_name)(iz, ix) = 0.25*(rho(iz, ix) + rho(iz + 1, ix) + rho(iz + 1, ix + 1) + rho(iz, ix + 1));
                else
                    lamda_half_x_levels.(level_name)(iz, ix) = lamda(iz, ix);
                    niu_half_x_levels.(level_name)(iz, ix) = niu(iz, ix);
                    niu_half_z_levels.(level_name)(iz, ix) = niu(iz, ix);
                    rho_half_x_half_z_levels.(level_name)(iz, ix) = rho(iz, ix);
                end
            end
        end
    end
    
    material_patches = cell(levels, 1);
    
    for il = 1:levels
        level_name = sprintf('level%d', il - 1);
    
        % Cache whole-level arrays once
        lamdaL    = lamda_levels.(level_name);
        niuL      = niu_levels.(level_name);
        rhoL      = rho_levels.(level_name);
        lamdaHxL  = lamda_half_x_levels.(level_name);
        niuHxL    = niu_half_x_levels.(level_name);
        niuHzL    = niu_half_z_levels.(level_name);
        rhoHxHzL  = rho_half_x_half_z_levels.(level_name);
    
        material_patches{il} = cell(patches(il), 1);
    
        for ip = 1:patches(il)
            patch_global_index = ip + patch_offset(il);
            patch_par = patch_catalog(patch_global_index);
    
            ix_start = patch_par.i_lo + 1;
            ix_end   = patch_par.i_hi + 1;
            iz_start = patch_par.j_lo + 1;
            iz_end   = patch_par.j_hi + 1;
    
            mp = struct();
            mp.lamda               = lamdaL(iz_start:iz_end, ix_start:ix_end);
            mp.niu                 = niuL(iz_start:iz_end, ix_start:ix_end);
            mp.rho                 = rhoL(iz_start:iz_end, ix_start:ix_end);
            mp.lamda_half_x        = lamdaHxL(iz_start:iz_end, ix_start:ix_end);
            mp.niu_half_x          = niuHxL(iz_start:iz_end, ix_start:ix_end);
            mp.niu_half_z          = niuHzL(iz_start:iz_end, ix_start:ix_end);
            mp.rho_half_x_half_z   = rhoHxHzL(iz_start:iz_end, ix_start:ix_end);
    
            material_patches{il}{ip} = mp;
        end
    end

    vmax = max(vp_levels.level2(:));