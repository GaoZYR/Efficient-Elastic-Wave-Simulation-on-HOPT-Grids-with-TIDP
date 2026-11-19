%% this function is used to add cpml boundarysurface
function [Memory_DxVx0, Memory_DzVz0, Memory_DxVz0, Memory_DzVx0, DxVx, DzVz, DxVz, DzVx] = ...
        add_cpml(Memory_DxVx0, Memory_DzVz0, Memory_DxVz0, Memory_DzVx0, DxVx, DzVz, DxVz, DzVx, ...
        ix_start, iz_start, nx, nz, N, index, ...
        a_x0, b_x0, a_z0, b_z0, dkx0, dkz0, ...
        a_x_half0, b_x_half0, a_z_half0, b_z_half0, dkx_half0, dkz_half0)

    a_x = a_x0(ix_start:ix_start + nx + 2*N - 1);
    b_x = b_x0(ix_start:ix_start + nx + 2*N - 1);
    dkx = dkx0(ix_start:ix_start + nx + 2*N - 1);
    a_x_half = a_x_half0(ix_start:ix_start + nx + 2*N - 1);
    b_x_half = b_x_half0(ix_start:ix_start + nx + 2*N - 1);
    dkx_half = dkx_half0(ix_start:ix_start + nx + 2*N - 1);

    a_z = a_z0(iz_start:iz_start + nz + 2*N - 1);
    b_z = b_z0(iz_start:iz_start + nz + 2*N - 1);
    dkz = dkz0(iz_start:iz_start + nz + 2*N - 1);
    a_z_half = a_z_half0(iz_start:iz_start + nz + 2*N - 1);
    b_z_half = b_z_half0(iz_start:iz_start + nz + 2*N - 1);
    dkz_half = dkz_half0(iz_start:iz_start + nz + 2*N - 1);

    Memory_DxVz = cell2mat(Memory_DxVz0(index));
    Memory_DzVz = cell2mat(Memory_DzVz0(index));
    Memory_DxVx = cell2mat(Memory_DxVx0(index));
    Memory_DzVx = cell2mat(Memory_DzVx0(index));

    % 等价关系 DxVz~DxTxx DzVz~DzTxz DxVx~DxTxz DzVx~DzTzz
    Memory_DxVz = b_x.*Memory_DxVz + a_x.*DxVz;
    Memory_DzVz = b_z'.*Memory_DzVz + a_z'.*DzVz;
    Memory_DxVx = b_x_half.*Memory_DxVx + a_x_half.*DxVx;
    Memory_DzVx = b_z_half'.*Memory_DzVx + a_z_half'.*DzVx;  

    DxVz = DxVz./dkx + Memory_DxVz;
    DzVz = DzVz./dkz' + Memory_DzVz;
    DxVx = DxVx./dkx_half + Memory_DxVx;
    DzVx = DzVx./dkz_half' + Memory_DzVx;

    % Memory_DxTxx(iz, ix) = b_x(ix)*Memory_DxTxx(iz, ix) + a_x(ix)*DxTxx(iz, ix);
    % Memory_DzTxz(iz, ix) = b_z(iz)*Memory_DzTxz(iz, ix) + a_z(iz)*DzTxz(iz, ix);
    % Memory_DxTxz(iz, ix) = b_x_half(ix)*Memory_DxTxz(iz, ix) + a_x_half(ix)*DxTxz(iz, ix);
    % Memory_DzTzz(iz, ix) = b_z_half(iz)*Memory_DzTzz(iz, ix) + a_z_half(iz)*DzTzz(iz, ix);
    % 
    % DxTxx(iz, ix) = DxTxx(iz, ix)/dkx(ix) + Memory_DxTxx(iz, ix);
    % DzTxz(iz, ix) = DzTxz(iz, ix)/dkz(iz) + Memory_DzTxz(iz, ix);
    % DxTxz(iz, ix) = DxTxz(iz, ix)/dkx_half(ix) + Memory_DxTxz(iz, ix);
    % DzTzz(iz, ix) = DzTzz(iz, ix)/dkz_half(iz) + Memory_DzTzz(iz, ix);

    Memory_DxVz0{index} = Memory_DxVz;
    Memory_DzVz0{index} = Memory_DzVz;
    Memory_DxVx0{index} = Memory_DxVx;
    Memory_DzVx0{index} = Memory_DzVx;