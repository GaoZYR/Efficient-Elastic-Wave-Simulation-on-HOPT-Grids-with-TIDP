function [Memory_DxVx, Memory_DzVz, Memory_DxVz, Memory_DzVx, ...
    Memory_DxTxz, Memory_DzTxz, Memory_DxTxx, Memory_DzTzz] = ...
    make_memory(patch_catalog, PML_id, N)

    % Find PML patches
    np_PML = size(PML_id, 1);
    PML_patch = repmat(struct('level', [], 'patch_id', [], ...
                              'x_lo', [], 'x_hi', [], 'z_lo', [], 'z_hi', [], ...
                              'nx', [], 'nz', [], 'dx', [], 'dz', [], ...
                              'i_lo', [], 'i_hi', [], 'j_lo', [], 'j_hi', []), 1, np_PML);
    index = 1;
    for ip = 1:numel(patch_catalog)
        lp = [patch_catalog(ip).level patch_catalog(ip).patch_id];
        is_in_PML = any(ismember(PML_id, lp, 'rows'));
        if is_in_PML
            PML_patch(index) = patch_catalog(ip);
            index = index + 1;
        end
    end

    % Make ceils
    Memory_DxVx = make_single_memory(PML_patch, np_PML, N); 
    Memory_DzVz = make_single_memory(PML_patch, np_PML, N);
    Memory_DxVz = make_single_memory(PML_patch, np_PML, N);
    Memory_DzVx = make_single_memory(PML_patch, np_PML, N);
    Memory_DxTxz = make_single_memory(PML_patch, np_PML, N);
    Memory_DzTxz = make_single_memory(PML_patch, np_PML, N);
    Memory_DxTxx = make_single_memory(PML_patch, np_PML, N);
    Memory_DzTzz = make_single_memory(PML_patch, np_PML, N);

end


function [Memory_DxVx] = make_single_memory(PML_patch, np_PML, N)

    % Make ceil
    Memory_DxVx = cell(1, np_PML); 

    % Initialize
    for ip = 1:np_PML
        patch = PML_patch(ip);
        nz = patch.nz;
        nx = patch.nx;
        Memory_DxVx{ip} = zeros(nz + 2*N, nx + 2*N);
    end
end