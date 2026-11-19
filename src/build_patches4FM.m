function [patches4FM] = build_patches4FM(patch_catalog, levels, patches, patch_offset, surface_id, PML_id, Nx, Nz, N)
    patches4FM = cell(levels, 1);
    
    for il = 1:levels
        patches4FM{il} = cell(patches(il), 1);
    
        for ip = 1:patches(il)
            patch_global_index = ip + patch_offset(il);
            patch_par = patch_catalog(patch_global_index);
    
            ps = struct();
    
            % original patch_par kept if later functions still need it
            ps.patch_par = patch_par;
    
            % geometry / indices
            ps.nx = patch_par.nx;
            ps.nz = patch_par.nz;
            ps.dx = patch_par.dx;
            ps.dz = patch_par.dz;
            ps.ix_start = patch_par.i_lo + 1;
            ps.ix_end   = patch_par.i_hi + 1;
            ps.iz_start = patch_par.j_lo + 1;
            ps.iz_end   = patch_par.j_hi + 1;
    
            % flags
            lp = [il - 1, ip - 1];
            ps.is_in_surface = any(ismember(surface_id, lp, 'rows'));
            ps.is_in_PML     = any(ismember(PML_id, lp, 'rows'));
    
            % precompute indices for exchange_points_on_same_levels
            ex_ix_start = patch_par.i_lo + 1 - N;
            ex_ix_end   = patch_par.i_hi + 1 + N;
            ex_iz_start = patch_par.j_lo + 1 - N;
            ex_iz_end   = patch_par.j_hi + 1 + N;
    
            ex_jx_start = 1; ex_jx_end   = patch_par.nx + 2*N;
            ex_jz_start = 1; ex_jz_end   = patch_par.nz + 2*N;
    
            if ex_ix_start < 1; ex_ix_start = 1; ex_jx_start = N + 1; end
            if ex_ix_end > Nx(il); ex_ix_end = Nx(il); ex_jx_end = patch_par.nx + N; end
            if ex_iz_start < 1; ex_iz_start = 1; ex_jz_start = N + 1; end
            if ex_iz_end > Nz(il); ex_iz_end = Nz(il); ex_jz_end = patch_par.nz + N; end
    
            ps.ex_ix_start = ex_ix_start; ps.ex_ix_end   = ex_ix_end;
            ps.ex_iz_start = ex_iz_start; ps.ex_iz_end   = ex_iz_end;
            ps.ex_jx_start = ex_jx_start; ps.ex_jx_end   = ex_jx_end;
            ps.ex_jz_start = ex_jz_start; ps.ex_jz_end   = ex_jz_end;
    
            patches4FM{il}{ip} = ps;
        end
    end