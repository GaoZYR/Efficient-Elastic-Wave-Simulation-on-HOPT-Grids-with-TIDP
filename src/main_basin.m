clc, clear, close all
%% Basic parameters

% Path 
base_dir = fullfile(fileparts(mfilename('fullpath')), '..', 'data');
patch_velocity_dir = fullfile(base_dir, 'model', 'patch_vs');
level_vp_dir = fullfile(base_dir, 'model', 'level_vp');
level_vs_dir = fullfile(base_dir, 'model', 'level_vs');
level_rho_dir = fullfile(base_dir, 'model', 'level_rho');

% F.D.M
dt = 0.0005; Nt = 40001;
Lx = 1000; Lz = 500;
N = 4; cn=solve_cn(N);

% Frequency of Ricker
fm = 3; 

% Whether set to free surface
free_surface = 1;

% Whether use wavefront activation patch
wavefront_flag = 1;
Tfir_dir = fullfile(base_dir, 'T_first_arrive');
t_safe_offset = 0;
%% Get numbers of levels and patches
[levels, patches, patch_catalog, patch_offset] = levels_patches_info(patch_velocity_dir);

patch00_par = patch_catalog(1); %%%%%%%%%%%%%% 
nx00 = 272; nz00 = 96; %%%%%%%%%%%%%%%%%%%%%%%%
dx00 = patch00_par.dx; dz00 = patch00_par.dz;

start_offset_x = patch00_par.x_lo - patch00_par.i_lo*dx00;
start_offset_z = patch00_par.z_lo - patch00_par.j_lo*dz00;

%% Initialize wavefields, build interp, restriction and T plans

[Vx_patches, Vz_patches, Txx_patches, Tzz_patches, Txz_patches, ...
    Vx_gather, Vz_gather, Txx_gather, Tzz_gather, Txz_gather, ...
    copy_plan, interp_plan, restriction_plan, Nx, Nz] = ...
    Initialize_wavefield_and_plan(patch_catalog, levels, patches, patch_offset, N, nx00, nz00);
%% Count total number of grid points

NNx = [patch_catalog.nx]; NNz = [patch_catalog.nz];
NN = 0;
for i = 1:length(NNx)
    NN = NN + NNx(i) * NNz(i);
end
fprintf("Total number of patch points in all levels is %d\n",NN);

%% Read model parameters

[material_patches, vmax] = ...
    read_model_parameters(level_vp_dir, level_vs_dir, level_rho_dir, ...
    levels, patches, patch_catalog, patch_offset, Nx, Nz);
%% Source & receiver
% The crood of src&rec is in raw model

% Set source location
src_amp = 0.5;
src_x = start_offset_x + 7791.666666304000 - 0.5*dx00/2^(levels-1);
src_z = start_offset_z + 379.5572916490000 - 0.5*dz00/2^(levels-1);
% src_x = 13500.0;
% src_z = 100.0;

% Set receivers location
% nrec = 3;
% rec_x = start_offset_x + [7791.666666304000 7791.666666304000 2814.453124869000] - 0.5*dx00/2^(levels-1);
% rec_z = start_offset_z + [2263.020833228000 3774.088541491000 3774.088541491000] - 0.5*dz00/2^(levels-1);
nrec = 280;
rec_x = 1000 + (50:50:14000) - 0.5*dx00/2^(levels-1);
rec_z = 0 + zeros(1, nrec);

% % Set source location
% src_amp = 0.5;
% src_x = 8000.000000000000 - 0.5*7.561436672967864;
% src_z = 377.7472527472527 - 0.5*7.554945054945055;
% 
% % Set receivers location
% nrec = 3;
% rec_x = [8000.000000000000 8000.000000000000 3024.574669187145] - 0.5*7.561436672967864;
% rec_z = [2266.483516483516 3777.472527472527 3777.472527472527] - 0.5*7.554945054945055;

% Find the patches that source and receivers in, and calculate Source
[source_patch, receiver_patch, rec_bilinear_plan, s, ricker, ricker_diff] = ...
    build_source_receivers(patch_catalog, src_x, src_z, rec_x, rec_z, fm, Nt, dt);

% [rec_index_x, rec_index_z] = record_location_index(receiver_patch, nrec, rec_x, rec_z, N);

% Define record
record_Vx = zeros(Nt-1, length(rec_x));
record_Vz = zeros(Nt-1, length(rec_x));

%% Check source & receivers

nx_raw = 161; nz_raw = 111; 
dx_raw = 100; dz_raw = 50;
xmax = (nx_raw - 1)*dx_raw;
zmax = (nz_raw - 1)*dz_raw;
xx_raw = (0:(nx_raw - 1))*dx_raw;
zz_raw = (0:(nz_raw - 1))*dz_raw;

vs_raw_fname = fullfile(base_dir, 'model', 'vs.bin');
fid = fopen(vs_raw_fname, 'rb');
vs_raw = fread(fid, [nx_raw, nz_raw], 'float32');
fclose(fid);

figure(1)
font_size = 12;
imagesc(xx_raw, zz_raw, vs_raw');
hold on

po = plot(rec_x, rec_z, 'Marker','v','LineStyle','none','MarkerFaceColor','g','MarkerEdgeColor','g','MarkerSize',6);
ps = plot(src_x, src_z, 'Marker','p','LineStyle','none','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',8);

npd_top = 0; npd_blr = 10;
xline(npd_blr*dx_raw, 'k', LineWidth=1)
xline((nx_raw-npd_blr)*dx_raw, 'k', LineWidth=1)
yline((nz_raw-npd_blr)*dz_raw, 'k', LineWidth=1)

legend([po, ps],'Receiver','Source','Location','NorthWest','FontSize',12);
axis equal tight
xlabel('Distance(m)');ylabel('Depth(m)');
set(gca,'FontName','Airl','FontSize',font_size);

close


%% Free surface & CPML boundary
force_retained_pml_tags = 0;
mpml_eta = 0.05; % 0: original CPML; 0.05: cross-direction stabilization


% Find the patch id of surface and CPML
[surface_id, PML_id] = find_surface_PML_patches(patch_catalog, levels, Nx, Nz, force_retained_pml_tags);

% Coefficients of CPML
if force_retained_pml_tags == 1
    nLx = ceil((Lx - patch00_par.x_lo)/patch_catalog(end).dx) + 1;
    nLz = ceil((Lz - (zmax - nz00 * dz00))/patch_catalog(end).dz) + 1;
    [ddx, ddx_half, dkx, dkx_half, dalphax, dalphax_half, ...
        ddz, ddz_half, dkz, dkz_half, dalphaz, dalphaz_half, ...
        a_x, a_x_half, b_x, b_x_half, a_z, a_z_half, b_z, b_z_half] = ...
        cpml_coeff(nLx + N, nLz + N, fm, vmax, Nx(levels) + 2*N, Nz(levels) + 2*N, patch_catalog(end).dx, dt, free_surface);
elseif force_retained_pml_tags == 0

    PML_level = unique(PML_id(:,1));

    dkx_levels = cell(numel(PML_level), 1);
    dkx_half_levels = cell(numel(PML_level), 1);
    dkz_levels = cell(numel(PML_level), 1);
    dkz_half_levels = cell(numel(PML_level), 1);
    a_x_levels = cell(numel(PML_level), 1);
    a_x_half_levels = cell(numel(PML_level), 1);
    a_z_levels = cell(numel(PML_level), 1);
    a_z_half_levels = cell(numel(PML_level), 1);
    b_x_levels = cell(numel(PML_level), 1);
    b_x_half_levels = cell(numel(PML_level), 1);
    b_z_levels = cell(numel(PML_level), 1);
    b_z_half_levels = cell(numel(PML_level), 1);

    for k = 1:numel(PML_level)
        il = PML_level(k);
        nLx = ceil((Lx - patch00_par.x_lo)/patch_catalog(end).dx) + 1;
        nLx = nLx/(2^((levels-1)-il));
        nLz = ceil((Lz - (zmax - nz00 * dz00))/patch_catalog(end).dz) + 1;
        nLz = nLz/(2^((levels-1)-il)); nLz = ceil(nLz);
        [ddx, ddx_half, dkx, dkx_half, dalphax, dalphax_half, ...
            ddz, ddz_half, dkz, dkz_half, dalphaz, dalphaz_half, ...
            a_x, a_x_half, b_x, b_x_half, a_z, a_z_half, b_z, b_z_half] = ...
            cpml_coeff(nLx + N, nLz + N, fm, vmax, Nx(levels)/(2^((levels-1)-il)) + 2*N, Nz(levels)/(2^((levels-1)-il)) + 2*N, patch_catalog(end).dx*(2^((levels-1)-il)), dt, free_surface);

        dkx_levels{k} = dkx;
        dkx_half_levels{k} = dkx_half;
        dkz_levels{k} = dkz;
        dkz_half_levels{k} = dkz_half;
        a_x_levels{k} = a_x;
        a_x_half_levels{k} = a_x_half;
        a_z_levels{k} = a_z;
        a_z_half_levels{k} = a_z_half;
        b_x_levels{k} = b_x;
        b_x_half_levels{k} = b_x_half;
        b_z_levels{k} = b_z;
        b_z_half_levels{k} = b_z_half;

    end
end

% Memory variable 
[Memory_DxVx, Memory_DzVz, Memory_DxVz, Memory_DzVx, ...
    Memory_DxTxz, Memory_DzTxz, Memory_DxTxx, Memory_DzTzz] = ...
    make_memory(patch_catalog, PML_id, N);


%%

patches4FM = build_patches4FM(patch_catalog, levels, patches, patch_offset, surface_id, PML_id, Nx, Nz, N);

if wavefront_flag == 1
    FMstart = find_wavefront_active_patches(Tfir_dir, patch_catalog, levels, patches, patch_offset, src_x, src_z, dt, t_safe_offset);
end
%% Finite difference Modeling in patches

time_stat = struct();

% Vx/Vz half
time_stat.V_step01 = 0;   
time_stat.V_step02 = 0;   
time_stat.V_step03 = 0;   
time_stat.V_step04 = 0;   

% Txx/Tzz/Txz half
time_stat.T_step01 = 0; 
time_stat.T_step02 = 0; 
time_stat.T_step03 = 0; 
time_stat.T_step04 = 0;  


% receiver record
time_stat.record = 0;

% optional: patch-level details
time_stat.V_step04_patch = 0;
time_stat.T_step04_patch = 0;

%% Snapshot settings
snapshot_dir = fullfile(fileparts(mfilename('fullpath')), '..', 'snapshot_basin');

if ~exist(snapshot_dir, 'dir')
    mkdir(snapshot_dir);
end

snapshot_t1 = 0.0;
snapshot_t2 = 10.0;
snapshot_interval = 0.02;

snapshot_step = round(snapshot_interval / dt);
snapshot_it1 = round(snapshot_t1 / dt);
snapshot_it2 = round(snapshot_t2 / dt);
%% Forward Modeling
t_total = tic;
for it = 1:Nt

    %% Vx Vz
    PMLpatchindex = 1;
    
    % Step 01: Get ghost points from same, coarser and finer levels
    tt = tic;
    [Vx_patches, Vz_patches] = apply_restriction_plan_V(restriction_plan, Vx_patches, Vz_patches);
    [Vx_patches, Vz_patches] = apply_copy_plan_V(copy_plan, Vx_patches, Vz_patches);
    [Vx_patches, Vz_patches] = apply_interp_plan_V(interp_plan, Vx_patches, Vz_patches);
    time_stat.V_step01 = time_stat.V_step01 + toc(tt);

    % Record receivers
    tt = tic;
    if it > 1
        record_Vx = record_receivers(Vx_patches, record_Vx, receiver_patch, rec_bilinear_plan, nrec, N, it - 1);
        record_Vz = record_receivers(Vz_patches, record_Vz, receiver_patch, rec_bilinear_plan, nrec, N, it - 1);
        % for ig = 1:nrec 
        %     rl = receiver_patch(ig).level;
        %     rp = receiver_patch(ig).patch_id;
        %     Vxrec = Vx_patches{rl+1}{rp+1};
        %     record_Vx(it-1, ig) = Vxrec(rec_index_z(ig), rec_index_x(ig));
        % end
    end
    time_stat.record = time_stat.record + toc(tt);


    for il = 1:levels
        for ip = 1:patches(il)

            % Step 02: Extract patch parameters
            tt = tic;
            ps = patches4FM{il}{ip};
            patch_par = ps.patch_par;
            nx = ps.nx; nz = ps.nz;
            dx = ps.dx; dz = ps.dz;
            ix_start = ps.ix_start; ix_end = ps.ix_end;
            iz_start = ps.iz_start; iz_end = ps.iz_end;
            is_in_surface = ps.is_in_surface;
            is_in_PML = ps.is_in_PML;

            Vx = Vx_patches{il}{ip};
            Vz = Vz_patches{il}{ip};
            Txx = Txx_patches{il}{ip};
            Tzz = Tzz_patches{il}{ip};
            Txz = Txz_patches{il}{ip};

            time_stat.V_step02 = time_stat.V_step02 + toc(tt);

            if wavefront_flag == 1
                if it < FMstart{il}{ip}
                    if is_in_PML
                        PMLpatchindex = PMLpatchindex + 1;
                    end
                    continue
                end
            end

            % Step 03: Read patch parameters & velocity model
            tt = tic;
            mp = material_patches{il}{ip};
            lamda = mp.lamda;
            niu = mp.niu;
            rho = mp.rho;
            lamda_half_x = mp.lamda_half_x;
            niu_half_x = mp.niu_half_x;
            niu_half_z = mp.niu_half_z;
            rho_half_x_half_z = mp.rho_half_x_half_z;
            time_stat.V_step03 = time_stat.V_step03 + toc(tt);

            % Step 04: Finite difference Modeling
            tt_patch = tic;
            DxVx = zeros(nz + 2*N, nx + 2*N); DzVz = zeros(nz + 2*N, nx + 2*N);
            DxVz = zeros(nz + 2*N, nx + 2*N); DzVx = zeros(nz + 2*N, nx + 2*N); 

            CV_Txx = zeros(nz + 2*N, nx + 2*N); CV_Tzz = zeros(nz + 2*N, nx + 2*N); CV_Txz = zeros(nz + 2*N, nx + 2*N);

            for iz = N + 1:nz + N
                for ix = N + 1:nx + N
                    DxVx(iz, ix) = 0; DzVz(iz, ix) = 0; DxVz(iz, ix) = 0; DzVx(iz, ix) = 0;
                    for m = 1:N
                        DxVx(iz, ix) = DxVx(iz, ix) + cn(m)*(Vx(iz, ix + m) - Vx(iz, ix - m + 1));
                        DzVz(iz, ix) = DzVz(iz, ix) + cn(m)*(Vz(iz - 1 + m, ix) - Vz(iz - m, ix));
                        DxVz(iz, ix) = DxVz(iz, ix) + cn(m)*(Vz(iz, ix - 1 + m) - Vz(iz, ix - m));
                        DzVx(iz, ix) = DzVx(iz, ix) + cn(m)*(Vx(iz + m, ix) - Vx(iz + 1 - m, ix));
                    end
                    DxVx(iz, ix) = DxVx(iz, ix)/dx;
                    DzVz(iz, ix) = DzVz(iz, ix)/dz;
                    DxVz(iz, ix) = DxVz(iz, ix)/dx;
                    DzVx(iz, ix) = DzVx(iz, ix)/dz;
                end
            end
            

            % Add CPML boundary
            if is_in_PML
                if force_retained_pml_tags == 0
                    ipl = find(PML_level == il - 1);
                    dkx = dkx_levels{ipl};
                    dkx_half = dkx_half_levels{ipl};
                    dkz = dkz_levels{ipl};
                    dkz_half = dkz_half_levels{ipl};
                    a_x = a_x_levels{ipl};
                    a_x_half = a_x_half_levels{ipl};
                    a_z = a_z_levels{ipl};
                    a_z_half = a_z_half_levels{ipl};
                    b_x = b_x_levels{ipl};
                    b_x_half = b_x_half_levels{ipl};
                    b_z = b_z_levels{ipl};
                    b_z_half = b_z_half_levels{ipl};
                end
                [Memory_DxVx, Memory_DzVz, Memory_DxVz, Memory_DzVx, DxVx, DzVz, DxVz, DzVx] = ...
                    add_mpml(Memory_DxVx, Memory_DzVz, Memory_DxVz, Memory_DzVx, DxVx, DzVz, DxVz, DzVx, ...
                    ix_start, iz_start, nx, nz, N, PMLpatchindex,...
                    a_x, b_x, a_z, b_z, dkx, dkz, ...
                    a_x_half, b_x_half, a_z_half, b_z_half, dkx_half, dkz_half, dt, mpml_eta, 'stress');
                PMLpatchindex = PMLpatchindex + 1;
            end

            for iz = N + 1:nz + N
                for ix = N + 1:nx + N
                    if free_surface == 1 && is_in_surface && iz == N + 1
                        CV_Txx(iz, ix) = 4*niu_half_x(iz - N, ix - N)*(lamda_half_x(iz - N, ix - N) + niu_half_x(iz - N, ix - N))/(lamda_half_x(iz - N, ix - N) + 2*niu_half_x(iz - N, ix - N))*DxVx(iz, ix);
                        CV_Tzz(iz, ix) = 0;
                        CV_Txz(iz, ix) = niu_half_z(iz - N, ix - N)*DzVx(iz, ix) + niu_half_z(iz - N, ix - N)*DxVz(iz, ix);    
                    else
                        CV_Txx(iz, ix) = (lamda_half_x(iz - N, ix - N) + 2*niu_half_x(iz - N, ix - N))*DxVx(iz, ix) + lamda_half_x(iz - N, ix - N)*DzVz(iz, ix);
                        CV_Tzz(iz, ix) = lamda_half_x(iz - N, ix - N)*DxVx(iz, ix) + (lamda_half_x(iz - N, ix - N) + 2*niu_half_x(iz - N, ix - N))*DzVz(iz, ix);
                        CV_Txz(iz, ix) = niu_half_z(iz - N, ix - N)*DzVx(iz, ix) + niu_half_z(iz - N, ix - N)*DxVz(iz, ix);
                    end
                end
            end

            Txx = Txx + dt*CV_Txx;
            Tzz = Tzz + dt*CV_Tzz;
            Txz = Txz + dt*CV_Txz;

            Txx_patches{il}{ip} = Txx;
            Txz_patches{il}{ip} = Txz;
            Tzz_patches{il}{ip} = Tzz;

            dt_patch = toc(tt_patch);
            time_stat.V_step04 = time_stat.V_step04 + dt_patch;
            time_stat.V_step04_patch = time_stat.V_step04_patch + dt_patch;
        end
    end

    %% Txx Tzz Txz
    PMLpatchindex = 1;

    % Step 01: Get ghost points from coarser levels
    tt = tic;
    [Txx_patches, Tzz_patches, Txz_patches] = apply_restriction_plan_T(restriction_plan, Txx_patches, Tzz_patches, Txz_patches);
    [Txx_patches, Tzz_patches, Txz_patches] = apply_copy_plan_T(copy_plan, Txx_patches, Tzz_patches, Txz_patches);
    [Txx_patches, Tzz_patches, Txz_patches] = apply_interp_plan_T(interp_plan, Txx_patches, Tzz_patches, Txz_patches);
    time_stat.T_step01 = time_stat.T_step01 + toc(tt);
 
    for il = 1:levels
        for ip = 1:patches(il)

            % Step 02: Extract patch parameters
            tt = tic;
            ps = patches4FM{il}{ip};
            patch_par = ps.patch_par;
            nx = ps.nx; nz = ps.nz;
            dx = ps.dx; dz = ps.dz;
            ix_start = ps.ix_start; ix_end = ps.ix_end;
            iz_start = ps.iz_start; iz_end = ps.iz_end;
            is_in_surface = ps.is_in_surface;
            is_in_PML = ps.is_in_PML;

            Vx = Vx_patches{il}{ip};
            Vz = Vz_patches{il}{ip};
            Txx = Txx_patches{il}{ip};
            Tzz = Tzz_patches{il}{ip};
            Txz = Txz_patches{il}{ip};
            time_stat.T_step02 = time_stat.T_step02 + toc(tt);

            if wavefront_flag == 1
                if it < FMstart{il}{ip}
                    if is_in_PML
                        PMLpatchindex = PMLpatchindex + 1;
                    end
                    continue
                end
            end

            % Step 03: Read patch parameters & velocity model
            tt = tic;
            mp = material_patches{il}{ip};
            lamda = mp.lamda;
            niu = mp.niu;
            rho = mp.rho;
            lamda_half_x = mp.lamda_half_x;
            niu_half_x = mp.niu_half_x;
            niu_half_z = mp.niu_half_z;
            rho_half_x_half_z = mp.rho_half_x_half_z;
            time_stat.T_step03 = time_stat.T_step03 + toc(tt);

            % Step 04: Finite difference Modeling
            tt_patch = tic;
            DxTxx = zeros(nz + 2*N, nx + 2*N); DzTzz = zeros(nz + 2*N, nx + 2*N); 
            DxTxz = zeros(nz + 2*N, nx + 2*N); DzTxz = zeros(nz + 2*N, nx + 2*N); 
    
            BS_Vx = zeros(nz + 2*N, nx + 2*N); BS_Vz = zeros(nz + 2*N, nx + 2*N);

            if free_surface == 1 && is_in_surface
                for in = 1:N 
                    Txx(in, :) = -Txx(2*N + 1 - in, :);
                    Tzz(in, :) = -Tzz(2*N + 1 - in, :);
                    Txz(in, :) = -Txz(2*N + 1 - in, :);
                end
            end

            % Add source
            if (il - 1) == source_patch.level && (ip - 1) == source_patch.patch_id
                Txx(N + 1:N + nz, N + 1:N + nx) = Txx(N + 1:N + nz, N + 1:N + nx) + src_amp * ricker(it) * s;
                Tzz(N + 1:N + nz, N + 1:N + nx) = Tzz(N + 1:N + nz, N + 1:N + nx) + src_amp * ricker(it) * s;
                
            end

            for iz = N + 1:nz + N
                for ix = N + 1:nx + N
                    DxTxx(iz, ix) = 0; DzTxz(iz, ix) = 0; DxTxz(iz, ix) = 0; DzTzz(iz, ix) = 0;
                    for m = 1:N
                        DxTxx(iz, ix) = DxTxx(iz, ix) + cn(m)*(Txx(iz, ix - 1 + m) - Txx(iz, ix - m));
                        DzTxz(iz, ix) = DzTxz(iz, ix) + cn(m)*(Txz(iz - 1 + m, ix) - Txz(iz - m, ix));
                        DxTxz(iz, ix) = DxTxz(iz, ix) + cn(m)*(Txz(iz, ix + m) - Txz(iz, ix - m + 1));
                        DzTzz(iz, ix) = DzTzz(iz, ix) + cn(m)*(Tzz(iz + m, ix) - Tzz(iz - m + 1, ix));
                    end
                    DxTxx(iz, ix) = DxTxx(iz, ix)/dx;
                    DzTxz(iz, ix) = DzTxz(iz, ix)/dz;
                    DxTxz(iz, ix) = DxTxz(iz, ix)/dx;
                    DzTzz(iz, ix) = DzTzz(iz, ix)/dz;
                end
            end

            % Add CPML boundary
            if is_in_PML
                if force_retained_pml_tags == 0
                    ipl = find(PML_level == il - 1);
                    dkx = dkx_levels{ipl};
                    dkx_half = dkx_half_levels{ipl};
                    dkz = dkz_levels{ipl};
                    dkz_half = dkz_half_levels{ipl};
                    a_x = a_x_levels{ipl};
                    a_x_half = a_x_half_levels{ipl};
                    a_z = a_z_levels{ipl};
                    a_z_half = a_z_half_levels{ipl};
                    b_x = b_x_levels{ipl};
                    b_x_half = b_x_half_levels{ipl};
                    b_z = b_z_levels{ipl};
                    b_z_half = b_z_half_levels{ipl};
                end
                [Memory_DxTxz, Memory_DzTxz, Memory_DxTxx, Memory_DzTzz, DxTxz, DzTxz, DxTxx, DzTzz] = ...
                    add_mpml(Memory_DxTxz, Memory_DzTxz, Memory_DxTxx, Memory_DzTzz, DxTxz, DzTxz, DxTxx, DzTzz, ...
                    ix_start, iz_start, nx, nz, N, PMLpatchindex,...
                    a_x, b_x, a_z, b_z, dkx, dkz, ...
                    a_x_half, b_x_half, a_z_half, b_z_half, dkx_half, dkz_half, dt, mpml_eta, 'velocity');
                PMLpatchindex = PMLpatchindex + 1;
            end

            for iz = N + 1:nz + N
                for ix = N + 1:nx + N    
                    BS_Vx(iz, ix) = 1/rho(iz - N, ix - N)*(DxTxx(iz, ix) + DzTxz(iz, ix));
                    BS_Vz(iz, ix) = 1/rho_half_x_half_z(iz - N, ix - N)*(DzTzz(iz, ix) + DxTxz(iz, ix));
                end
            end
        
            Vx = Vx + dt.*BS_Vx;
            Vz = Vz + dt.*BS_Vz;
            
            Vx_patches{il}{ip} = Vx;
            Vz_patches{il}{ip} = Vz;
            Txx_patches{il}{ip} = Txx;
            Txz_patches{il}{ip} = Txz;
            Tzz_patches{il}{ip} = Tzz;

            dt_patch = toc(tt_patch);
            time_stat.T_step04 = time_stat.T_step04 + dt_patch;
            time_stat.T_step04_patch = time_stat.T_step04_patch + dt_patch;

        end
    end

    if mod(it,5) == 0
        % Vx_gather = gather_wavefield_patch(Vx_patches, Vx_gather, patch_catalog, levels, patches, patch_offset, N);
        Vz_gather = gather_wavefield_patch(Vz_patches, Vz_gather, patch_catalog, levels, patches, patch_offset, N);
        imagesc((0:(nx00-1)*dx00)./1000, (0:(nz00-1)*dz00)./1000, Vz_gather);
        colormap('gray');
        colorbar;
        axis equal tight
        xlabel('Distance (km)','FontName','Airl','FontSize',13);
        ylabel('Depth (km)','FontName','Airl','FontSize',13);
        title("Time = " + it);
        drawnow
    end


    %% Save wavefield snapshots
    % if it >= snapshot_it1 && it <= snapshot_it2 && ...
    %         mod(it, snapshot_step) == 0
    % 
    %     % Calculate Vx and Vz gathered wavefields separately
    %     Vx_gather = gather_wavefield_patch( ...
    %         Vx_patches, Vx_gather, patch_catalog, ...
    %         levels, patches, patch_offset, N);
    % 
    %     Vz_gather = gather_wavefield_patch( ...
    %         Vz_patches, Vz_gather, patch_catalog, ...
    %         levels, patches, patch_offset, N);
    % 
    %     % Snapshot time
    %     snapshot_time = it * dt;
    % 
    %     % Save both components and grid information
    %     snapshot_file = fullfile(snapshot_dir, ...
    %         sprintf('snapshot_%06d.mat', it));
    % 
    %     save(snapshot_file, ...
    %         'Vx_gather', 'Vz_gather', ...
    %         'snapshot_time', 'it', ...
    %         'dx00', 'dz00', ...
    %         'start_offset_x', 'start_offset_z', ...
    %         '-v7.3');
    % 
    %     fprintf('Saved snapshot at t = %.3f s\n', snapshot_time);
    % end

end
total_time = toc(t_total);


%% Timing Summary
fprintf('\n=========== Timing Summary ===========\n');
fprintf('Total wall time: %.6f s\n\n', total_time);

names = { ...
    'V_step01', ...
    'V_step02', ...
    'V_step03', ...
    'V_step04', ...
    'T_step01', ...
    'T_step02', ...
    'T_step03', ...
    'T_step04', ...
    'record receivers'};

vals = [ ...
    time_stat.V_step01, ...
    time_stat.V_step02, ...
    time_stat.V_step03, ...
    time_stat.V_step04, ...
    time_stat.T_step01, ...
    time_stat.T_step02, ...
    time_stat.T_step03, ...
    time_stat.T_step04, ...
    time_stat.record];

for k = 1:numel(vals)
    fprintf('%-32s : %12.6f s   (%6.2f%%)   avg/step = %.6e s\n', ...
        names{k}, vals(k), 100*vals(k)/total_time, vals(k)/Nt);
end
fprintf('======================================\n');
