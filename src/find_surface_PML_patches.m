function [surface_id, PML_id] = find_surface_PML_patches(patch_catalog, levels, Nx, Nz, force_retained_pml_tags)
    np = size(patch_catalog, 2);
    surface_id = [];
    PML_id = [];
    for ip = 1:np
        patch = patch_catalog(ip);
        level_id = patch.level;
        patch_id = patch.patch_id;
        ix_start = patch.i_lo + 1;
        iz_start = patch.j_lo + 1;
        ix_end = patch.i_hi + 1;
        iz_end = patch.j_hi + 1;


        if iz_start == 1
            surface_id = [surface_id; level_id patch_id];
        end

        is_pml_patch = (ix_start == 1 || ix_end == Nx(level_id + 1) || iz_end == Nz(level_id + 1));
        if  is_pml_patch && (force_retained_pml_tags == 0 || (force_retained_pml_tags == 1 && level_id == levels - 1))
            PML_id = [PML_id; level_id patch_id];
        end

     end

