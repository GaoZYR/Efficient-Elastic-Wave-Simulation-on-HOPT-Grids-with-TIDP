function [Vx_patches, Vz_patches, Txx_patches, Tzz_patches, Txz_patches, ...
    Vx_gather, Vz_gather, Txx_gather, Tzz_gather, Txz_gather, ...
    copy_plan, interp_plan, restriction_plan, Nx, Nz] = ...
    Initialize_wavefield_and_plan_new(patch_catalog, levels, patches, patch_offset, N, nx00, nz00)

    Vx_gather = nan(2^(levels - 1)*nz00, 2^(levels - 1)*nx00);
    Vz_gather = nan(2^(levels - 1)*nz00, 2^(levels - 1)*nx00);
    Txx_gather = nan(2^(levels - 1)*nz00, 2^(levels - 1)*nx00);
    Tzz_gather = nan(2^(levels - 1)*nz00, 2^(levels - 1)*nx00);
    Txz_gather = nan(2^(levels - 1)*nz00, 2^(levels - 1)*nx00);

    copy_plan = cell(levels, 1);
    interp_plan = cell(levels, 1);
    restriction_plan = cell(levels, 1);
    ghost_info_all = cell(levels, 1);

    for il = 1:levels

        Nz(il) = 2^(il - 1)*nz00;
        Nx(il) = 2^(il - 1)*nx00;

        Vx_patches{il} = cell(patches(il), 1);
        Vz_patches{il} = cell(patches(il), 1);
        Txx_patches{il} = cell(patches(il), 1);
        Tzz_patches{il} = cell(patches(il), 1);
        Txz_patches{il} = cell(patches(il), 1);

        copy_plan{il} = cell(patches(il), 1);
        interp_plan{il} = cell(patches(il), 1);
        restriction_plan{il} = cell(patches(il), 1);
        ghost_info_all{il} = cell(patches(il), 1);

        for ip = 1:patches(il)
            
            % Patch Parameters
            patch_par = patch_catalog(ip + patch_offset(il));
            ix_start = patch_par.i_lo + 1; ix_end = patch_par.i_hi + 1;
            iz_start = patch_par.j_lo + 1; iz_end = patch_par.j_hi + 1;
            nx = patch_par.nx; nz = patch_par.nz;

            % Initialize wavefield patches
            Vx_patches{il}{ip} = zeros(nz + 2*N, nx + 2*N);
            Vz_patches{il}{ip} = zeros(nz + 2*N, nx + 2*N);
            Txx_patches{il}{ip} = zeros(nz + 2*N, nx + 2*N);
            Tzz_patches{il}{ip} = zeros(nz + 2*N, nx + 2*N);
            Txz_patches{il}{ip} = zeros(nz + 2*N, nx + 2*N);
            

            % ============================================================
            % Build copy, interp, restriction plan 
            % type:
            %   0  : 真点(internal)
            %   1  : 同层 patch
            %   2  : 更细层 patch
            %   3  : 更粗层 patch
            %  -1  : 不落在任何 patch 中（物理边界/外部）
            % ============================================================

            nz_loc = nz + 2*N;
            nx_loc = nx + 2*N;

            ghost_type = -ones(nz_loc, nx_loc);
            src_level  = -ones(nz_loc, nx_loc); 
            src_patch  = -ones(nz_loc, nx_loc);  
            x_phys     = nan(nz_loc, nx_loc);
            z_phys     = nan(nz_loc, nx_loc);

            ghost_mask = false(nz_loc, nx_loc);
            ghost_mask(1:N, N+1:nx+N) = true;
            ghost_mask(nz+N+1:end, N+1:nx+N) = true;
            ghost_mask(N+1:nz+N, 1:N) = true;
            ghost_mask(N+1:nz+N, nx+N+1:end) = true;

            ghost_type(N+1:nz+N, N+1:nx+N) = 0;

            for jloc = 1:nz_loc
                for iloc = 1:nx_loc

                    x = patch_par.x_lo + (iloc - 1 - N) * patch_par.dx + 0.5 * patch_par.dx;
                    z = patch_par.z_lo + (jloc - 1 - N) * patch_par.dz + 0.5 * patch_par.dz;

                    x_phys(jloc, iloc) = x;
                    z_phys(jloc, iloc) = z;

                    if ~ghost_mask(jloc, iloc)
                        continue;
                    end

                    [ptype, slevel, spatch] = classify_ghost_point( ...
                        x, z, il, ip, patch_catalog, levels, patches, patch_offset);

                    ghost_type(jloc, iloc) = ptype;
                    src_level(jloc, iloc) = slevel;
                    src_patch(jloc, iloc) = spatch;
                end
            end
            
            ghost_info = struct();
            ghost_info.ghost_mask = ghost_mask;
            ghost_info.ghost_type = ghost_type;
            ghost_info.src_level = src_level;
            ghost_info.src_patch = src_patch;
            ghost_info.x_phys = x_phys;
            ghost_info.z_phys = z_phys;
            ghost_info.type_note = '0=interior, 1=same_level, 2=finer_level, 3=coarser_level, -1=outside';

            ghost_info_all{il}{ip} = ghost_info;
            
            copy_plan{il}{ip} = build_copy_plan_one_patch(ghost_info, patch_par, patch_catalog, patch_offset, il, N);
            % interp_plan{il}{ip} = [];
            % restriction_plan{il}{ip} = [];
        end
    end
    % Vx_new_patches = Vx_patches;
    % Vz_new_patches = Vz_patches;
    % Txx_new_patches = Txx_patches;
    % Tzz_new_patches = Tzz_patches;
    % Txz_new_patches = Txz_patches;

    for il = 1:levels
        for ip = 1:patches(il)
            patch_par = patch_catalog(ip + patch_offset(il));
    
            restriction_plan{il}{ip} = build_restriction_plan_one_patch( ...
                ghost_info_all{il}{ip}, patch_par, patch_catalog, patch_offset, ...
                patches, ghost_info_all, il, levels, N);
    
            interp_plan{il}{ip} = build_interp_plan_one_patch( ...
                ghost_info_all{il}{ip}, patch_par, patch_catalog, patch_offset, ...
                patches, ghost_info_all, il, levels, Nx, Nz, N);
        end
    end
end
%% 辅助函数 1：判断 ghost 点性质
function [ptype, slevel, spatch] = classify_ghost_point( ...
    x, z, il, ip, patch_catalog, levels, patches, patch_offset)

    % 默认：外部/物理边界
    ptype = -1;
    slevel = -1;
    spatch = -1;

    % --------------------------------------------------
    % 1) 同层优先
    % --------------------------------------------------
    for ip2 = 1:patches(il)
        if ip2 == ip
            continue;
        end

        patch2 = patch_catalog(ip2 + patch_offset(il));
        if point_in_patch(x, z, patch2)
            ptype = 1;
            slevel = patch2.level;
            spatch = patch2.patch_id;
            return;
        end
    end

    % --------------------------------------------------
    % 2) 更细层
    % --------------------------------------------------
    for il2 = il+1 : levels
        for ip2 = 1:patches(il2)
            patch2 = patch_catalog(ip2 + patch_offset(il2));
            if point_in_patch(x, z, patch2)
                ptype = 2;
                slevel = patch2.level;
                spatch = patch2.patch_id;
                return;
            end
        end
    end

    % --------------------------------------------------
    % 3) 更粗层
    % --------------------------------------------------
    for il2 = il-1 : -1 : 1
        for ip2 = 1:patches(il2)
            patch2 = patch_catalog(ip2 + patch_offset(il2));
            if point_in_patch(x, z, patch2)
                ptype = 3;
                slevel = patch2.level;
                spatch = patch2.patch_id;
                return;
            end
        end
    end
end
%% 辅助函数 2：判断一个物理点是否落在某个 patch 内
function inside = point_in_patch(x, z, patch_par)

    tol = 1e-10;

    % x_min = patch_par.x_lo + 0.5 * patch_par.dx - tol;
    % x_max = patch_par.x_hi - 0.5 * patch_par.dx + tol;
    % z_min = patch_par.z_lo + 0.5 * patch_par.dz - tol;
    % z_max = patch_par.z_hi - 0.5 * patch_par.dz + tol;
    x_min = patch_par.x_lo - tol;
    x_max = patch_par.x_hi + tol;
    z_min = patch_par.z_lo - tol;
    z_max = patch_par.z_hi + tol;

    inside = (x >= x_min) && (x <= x_max) && (z >= z_min) && (z <= z_max);
end
%% 辅助函数 3：build copy plan one patch
function plan = build_copy_plan_one_patch(ghost_info, dst_patch_par, patch_catalog, patch_offset, il, N)

    % 找到本 patch 里所有需要 same-level copy 的 ghost 点
    mask = (ghost_info.ghost_type == 1);
    [jloc_all, iloc_all] = find(mask);

    if isempty(iloc_all)
        plan.nblock = 0;
        plan.src_patch = [];
        plan.dst_ind = {};
        plan.src_ind = {};
        return;
    end

    % 目标点在当前 level 上的“全局网格索引”(1-based)
    ix_global = dst_patch_par.i_lo + (iloc_all - N);
    iz_global = dst_patch_par.j_lo + (jloc_all - N);

    % 这些 ghost 各自对应的 donor patch（你现在存的是 0-based patch_id）
    src_patch0_all = ghost_info.src_patch(mask);   % 0-based

    % 按 donor patch 分组，减少时间步里的循环次数
    [src_patch0_unique, ~, grp] = unique(src_patch0_all, 'stable');

    plan.nblock = numel(src_patch0_unique);
    plan.src_patch = src_patch0_unique(:);
    plan.dst_ind = cell(plan.nblock, 1);
    plan.src_ind = cell(plan.nblock, 1);

    for k = 1:plan.nblock
        m = (grp == k);

        sp0 = src_patch0_unique(k);   % 0-based
        src_patch_par = patch_catalog(sp0 + 1 + patch_offset(il));

        % donor patch 内局部索引
        iloc_src = ix_global(m) - src_patch_par.i_lo + N;
        jloc_src = iz_global(m) - src_patch_par.j_lo + N;

        % 目标 patch 内线性索引
        plan.dst_ind{k} = sub2ind([dst_patch_par.nz + 2*N, dst_patch_par.nx + 2*N], ...
                                  jloc_all(m), iloc_all(m));

        % donor patch 内线性索引
        plan.src_ind{k} = sub2ind([src_patch_par.nz + 2*N, src_patch_par.nx + 2*N], ...
                                  jloc_src, iloc_src);
    end
end
%% 辅助函数 4：build restriction plan one patch
function plan = build_restriction_plan_one_patch( ...
    ghost_info, dst_patch_par, patch_catalog, patch_offset, patches, ghost_info_all, il, levels, N)

    plan = empty_restriction_plan();

    % 最细层没有更细 donor
    if il == levels
        return;
    end

    % ------------------------------------------------------------
    % restriction 的目标不再只是 ghost_type == 2
    % 还包括：当前 coarse patch 内部、但被更细层 patch 覆盖的 coarse 点
    %
    % fill_mask = finer-ghost  +  finer-covered interior
    % ------------------------------------------------------------
    fill_mask = build_restriction_fill_mask_one_patch( ...
        ghost_info, dst_patch_par, patch_catalog, patch_offset, patches, il, levels, N);

    [jloc_all, iloc_all] = find(fill_mask);

    if isempty(iloc_all)
        return;
    end

    % 五个变量分别建 composite restriction plan
    plan.Vx  = build_restriction_var_plan('Vx',  jloc_all, iloc_all, dst_patch_par, ...
                                          patch_catalog, patch_offset, patches, ghost_info_all, il, levels, N);

    plan.Vz  = build_restriction_var_plan('Vz',  jloc_all, iloc_all, dst_patch_par, ...
                                          patch_catalog, patch_offset, patches, ghost_info_all, il, levels, N);

    plan.Txx = build_restriction_var_plan('Txx', jloc_all, iloc_all, dst_patch_par, ...
                                          patch_catalog, patch_offset, patches, ghost_info_all, il, levels, N);

    plan.Tzz = build_restriction_var_plan('Tzz', jloc_all, iloc_all, dst_patch_par, ...
                                          patch_catalog, patch_offset, patches, ghost_info_all, il, levels, N);

    plan.Txz = build_restriction_var_plan('Txz', jloc_all, iloc_all, dst_patch_par, ...
                                          patch_catalog, patch_offset, patches, ghost_info_all, il, levels, N);
end

function fill_mask = build_restriction_fill_mask_one_patch( ...
    ghost_info, dst_patch_par, patch_catalog, patch_offset, patches, il, levels, N)

    nz_loc = dst_patch_par.nz + 2*N;
    nx_loc = dst_patch_par.nx + 2*N;

    fill_mask = false(nz_loc, nx_loc);

    % ------------------------------------------------------------
    % A. 原来的 finer-ghost
    % ------------------------------------------------------------
    fill_mask = fill_mask | (ghost_info.ghost_type == 2);

    % ------------------------------------------------------------
    % B. patch 内部真点，但几何上被更细层 patch 覆盖
    % ------------------------------------------------------------
    interior_mask = false(nz_loc, nx_loc);
    interior_mask(N+1:N+dst_patch_par.nz, N+1:N+dst_patch_par.nx) = true;

    for jloc = N+1:N+dst_patch_par.nz
        for iloc = N+1:N+dst_patch_par.nx

            x = ghost_info.x_phys(jloc, iloc);
            z = ghost_info.z_phys(jloc, iloc);

            if point_in_any_finer_patch(x, z, il, patch_catalog, patch_offset, patches, levels)
                fill_mask(jloc, iloc) = true;
            end
        end
    end
end

function inside = point_in_any_finer_patch(x, z, il, patch_catalog, patch_offset, patches, levels)

    inside = false;

    for il2 = il+1:levels
        for ip2 = 1:patches(il2)
            patch2 = patch_catalog(ip2 + patch_offset(il2));
            if point_in_patch(x, z, patch2)
                inside = true;
                return;
            end
        end
    end
end

function vplan = build_restriction_var_plan(var_type, jloc_all, iloc_all, dst_patch_par, ...
                                            patch_catalog, patch_offset, patches, ghost_info_all, il, levels, N)

    % 临时 block 容器
    blocks = struct('src_level', {}, 'src_patch', {}, ...
                    'dst_list', {}, 'src_rows', {}, 'w_rows', {});

    for m = 1:numel(jloc_all)

        jloc = jloc_all(m);
        iloc = iloc_all(m);

        % 当前 coarse 点在本层的全局网格索引
        IC = dst_patch_par.i_lo + (iloc - N);
        JC = dst_patch_par.j_lo + (jloc - N);

        % 目标 patch 局部线性索引
        dst_ind = sub2ind([dst_patch_par.nz + 2*N, dst_patch_par.nx + 2*N], jloc, iloc);

        % 递归展开成最终 donor atoms
        atoms = expand_restriction_target(var_type, il, JC, IC, ...
                                          patch_catalog, patch_offset, patches, ghost_info_all, levels, N);

        if isempty(atoms.level0)
            error(['Failed to expand restriction target for var=%s, ' ...
                   'coarse patch(level=%d, patch=%d), local ghost(%d,%d).'], ...
                  var_type, dst_patch_par.level, dst_patch_par.patch_id, jloc, iloc);
        end

        % 先按 (src_level, src_patch) 分组，再把这一行塞进对应 block
        key = [atoms.level0(:), atoms.patch0(:)];
        [key_u, ~, grp] = unique(key, 'rows', 'stable');

        for g = 1:size(key_u, 1)
            src_level0 = key_u(g, 1);
            src_patch0 = key_u(g, 2);

            src_row = atoms.src_ind(grp == g).';
            w_row   = atoms.w(grp == g).';

            blocks = append_restriction_block_row(blocks, src_level0, src_patch0, ...
                                                  dst_ind, src_row, w_row);
        end
    end

    vplan = finalize_restriction_blocks(blocks);
end

function atoms = expand_restriction_target(var_type, il, JC, IC, ...
                                           patch_catalog, patch_offset, patches, ghost_info_all, levels, N)

    atoms = empty_atoms();

    % 先生成相邻 finer level 上的 restriction children
    [childJ, childI, childW] = get_restriction_children(var_type, JC, IC);

    fine_il = il + 1;
    for q = 1:numel(childW)
        atoms_q = resolve_point_on_finer_hierarchy(var_type, fine_il, ...
                                                   childJ(q), childI(q), childW(q), ...
                                                   patch_catalog, patch_offset, patches, ...
                                                   ghost_info_all, levels, N);
        atoms = concat_atoms(atoms, atoms_q);
    end

    atoms = merge_atoms(atoms);
end

function atoms = resolve_point_on_finer_hierarchy(var_type, lev, J, I, base_w, ...
                                                  patch_catalog, patch_offset, patches, ...
                                                  ghost_info_all, levels, N)

    atoms = empty_atoms();

    [found, patch0, src_ind] = find_holder_patch_on_level( ...
        lev, J, I, patch_catalog, patch_offset, patches, ghost_info_all, N);

    if found
        atoms.level0 = lev - 1;   % 0-based
        atoms.patch0 = patch0;    % 0-based
        atoms.src_ind = src_ind;
        atoms.w = base_w;
        return;
    end

    % 到最细层还找不到，就失败
    if lev == levels
        return;
    end

    % 当前层找不到 holder，就继续往更细层展开
    [childJ, childI, childW] = get_restriction_children(var_type, J, I);

    for q = 1:numel(childW)
        atoms_q = resolve_point_on_finer_hierarchy(var_type, lev + 1, ...
                                                   childJ(q), childI(q), base_w * childW(q), ...
                                                   patch_catalog, patch_offset, patches, ...
                                                   ghost_info_all, levels, N);
        atoms = concat_atoms(atoms, atoms_q);
    end

    atoms = merge_atoms(atoms);
end

function [found, patch0, src_ind] = find_holder_patch_on_level( ...
    lev, J, I, patch_catalog, patch_offset, patches, ghost_info_all, N)

    % 接受的 holder 类型：
    %   0 interior
    %   1 same-level copy ghost
    %   2 finer-donor ghost（因为运行时 finest->coarsest，这类点在此之前已会被填好）
    %
    % 不接受：
    %  -1 outside
    %   3 coarser-donor ghost（依赖 interp，restriction 时还不可用）

    found = false;
    patch0 = -1;
    src_ind = [];

    best_rank = inf;

    for ip2 = 1:patches(lev)

        patch_par = patch_catalog(ip2 + patch_offset(lev));

        iloc = I - patch_par.i_lo + N;
        jloc = J - patch_par.j_lo + N;

        % 不在该 patch 的扩展数组范围里
        if iloc < 1 || iloc > patch_par.nx + 2*N || ...
           jloc < 1 || jloc > patch_par.nz + 2*N
            continue;
        end

        % interior
        if iloc >= N+1 && iloc <= N+patch_par.nx && ...
           jloc >= N+1 && jloc <= N+patch_par.nz
            gtype = 0;
        else
            gtype = ghost_info_all{lev}{ip2}.ghost_type(jloc, iloc);
        end

        rank = holder_rank(gtype);

        if rank < best_rank
            best_rank = rank;
            found = true;
            patch0 = patch_par.patch_id;   % 0-based
            src_ind = sub2ind([patch_par.nz + 2*N, patch_par.nx + 2*N], jloc, iloc);

            % interior 已经是最优，直接停
            if rank == 1
                return;
            end
        end
    end
end

function rank = holder_rank(gtype)
    switch gtype
        case 0
            rank = 1;   % interior best
        case 2
            rank = 2;   % only if this point is already restriction-filled from finer levels
        otherwise
            rank = inf; % -1 / 1 / 3 都不接受
    end
end

function [childJ, childI, childW] = get_restriction_children(var_type, J, I)

    switch var_type
        case 'Vx'
            childJ = [2*J - 1];
            childI = [2*I - 1];
            childW = [1.0];

        case {'Txx', 'Tzz'}
            childJ = [2*J - 1, 2*J - 1];
            childI = [2*I - 1, 2*I];
            childW = [0.5, 0.5];

        case 'Txz'
            childJ = [2*J - 1, 2*J];
            childI = [2*I - 1, 2*I - 1];
            childW = [0.5, 0.5];

        case 'Vz'
            childJ = [2*J - 1, 2*J - 1, 2*J, 2*J];
            childI = [2*I - 1, 2*I, 2*I - 1, 2*I];
            childW = [0.25, 0.25, 0.25, 0.25];

        otherwise
            error('Unknown var_type: %s', var_type);
    end
end



function plan = empty_restriction_plan()
    plan.Vx  = empty_restriction_var_plan();
    plan.Vz  = empty_restriction_var_plan();
    plan.Txx = empty_restriction_var_plan();
    plan.Tzz = empty_restriction_var_plan();
    plan.Txz = empty_restriction_var_plan();
end

function vplan = empty_restriction_var_plan()
    vplan.nblock = 0;
    vplan.src_level = [];   % 0-based donor level
    vplan.src_patch = [];   % 0-based donor patch
    vplan.dst_ind = {};
    vplan.src_ind = {};
    vplan.w = {};
    vplan.nsrc = {};
    vplan.dst_union = [];  
end

function blocks = append_restriction_block_row(blocks, src_level0, src_patch0, dst_ind, src_row, w_row)

    hit = 0;
    for k = 1:numel(blocks)
        if blocks(k).src_level == src_level0 && blocks(k).src_patch == src_patch0
            hit = k;
            break;
        end
    end

    if hit == 0
        hit = numel(blocks) + 1;
        blocks(hit).src_level = src_level0;
        blocks(hit).src_patch = src_patch0;
        blocks(hit).dst_list = [];
        blocks(hit).src_rows = {};
        blocks(hit).w_rows = {};
    end

    blocks(hit).dst_list(end+1, 1) = dst_ind;
    blocks(hit).src_rows{end+1, 1} = src_row(:).';
    blocks(hit).w_rows{end+1, 1} = w_row(:).';
end

function vplan = finalize_restriction_blocks(blocks)

    if isempty(blocks)
        vplan = empty_restriction_var_plan();
        return;
    end

    nblock = numel(blocks);

    vplan.nblock = nblock;
    vplan.src_level = zeros(nblock, 1);
    vplan.src_patch = zeros(nblock, 1);
    vplan.dst_ind = cell(nblock, 1);
    vplan.src_ind = cell(nblock, 1);
    vplan.w = cell(nblock, 1);
    vplan.nsrc = cell(nblock, 1);
    vplan.dst_union = [];

    for k = 1:nblock
        vplan.src_level(k) = blocks(k).src_level;
        vplan.src_patch(k) = blocks(k).src_patch;

        M = numel(blocks(k).dst_list);
        nsrc_row = zeros(M, 1);
        Smax = 0;

        for r = 1:M
            nsrc_row(r) = numel(blocks(k).src_rows{r});
            if nsrc_row(r) > Smax
                Smax = nsrc_row(r);
            end
        end

        src_mat = ones(M, Smax);   % 注意这里改成 ones
        w_mat = zeros(M, Smax);

        for r = 1:M
            s = nsrc_row(r);
            src_mat(r, 1:s) = blocks(k).src_rows{r};
            w_mat(r, 1:s) = blocks(k).w_rows{r};
        end

        vplan.dst_ind{k} = blocks(k).dst_list;
        vplan.src_ind{k} = src_mat;
        vplan.w{k} = w_mat;
        vplan.nsrc{k} = nsrc_row;
    end

    all_dst = [];
    for k = 1:nblock
        all_dst = [all_dst; vplan.dst_ind{k}(:)];
    end
    vplan.dst_union = unique(all_dst);
end

function atoms = empty_atoms()
    atoms.level0 = [];
    atoms.patch0 = [];
    atoms.src_ind = [];
    atoms.w = [];
end

function atoms = concat_atoms(a, b)
    atoms.level0 = [a.level0(:); b.level0(:)];
    atoms.patch0 = [a.patch0(:); b.patch0(:)];
    atoms.src_ind = [a.src_ind(:); b.src_ind(:)];
    atoms.w = [a.w(:); b.w(:)];
end

function atoms = merge_atoms(atoms)

    if isempty(atoms.level0)
        return;
    end

    key = [atoms.level0(:), atoms.patch0(:), atoms.src_ind(:)];
    [key_u, ~, grp] = unique(key, 'rows', 'stable');

    w_sum = accumarray(grp, atoms.w(:));

    atoms.level0 = key_u(:, 1);
    atoms.patch0 = key_u(:, 2);
    atoms.src_ind = key_u(:, 3);
    atoms.w = w_sum;
end

%% 辅助函数 5：build interp plan one patch
function plan = build_interp_plan_one_patch( ...
    ghost_info, dst_patch_par, patch_catalog, patch_offset, ...
    patches, ghost_info_all, il, levels, Nx, Nz, N)

    % coarsest level 没有更粗 donor
    plan = empty_restriction_plan();   % 直接复用 restriction 的 plan 结构

    if il == 1
        return;
    end

    % interp 只处理 ghost_type == 3
    mask = (ghost_info.ghost_type == 3);
    [jloc_all, iloc_all] = find(mask);

    if isempty(iloc_all)
        return;
    end

    plan.Vx  = build_interp_var_plan('Vx',  jloc_all, iloc_all, dst_patch_par, ...
                                     patch_catalog, patch_offset, patches, ...
                                     ghost_info_all, il, levels, Nx, Nz, N);

    plan.Vz  = build_interp_var_plan('Vz',  jloc_all, iloc_all, dst_patch_par, ...
                                     patch_catalog, patch_offset, patches, ...
                                     ghost_info_all, il, levels, Nx, Nz, N);

    plan.Txx = build_interp_var_plan('Txx', jloc_all, iloc_all, dst_patch_par, ...
                                     patch_catalog, patch_offset, patches, ...
                                     ghost_info_all, il, levels, Nx, Nz, N);

    plan.Tzz = build_interp_var_plan('Tzz', jloc_all, iloc_all, dst_patch_par, ...
                                     patch_catalog, patch_offset, patches, ...
                                     ghost_info_all, il, levels, Nx, Nz, N);

    plan.Txz = build_interp_var_plan('Txz', jloc_all, iloc_all, dst_patch_par, ...
                                     patch_catalog, patch_offset, patches, ...
                                     ghost_info_all, il, levels, Nx, Nz, N);
end

function vplan = build_interp_var_plan(var_type, jloc_all, iloc_all, dst_patch_par, ...
                                       patch_catalog, patch_offset, patches, ...
                                       ghost_info_all, il, levels, Nx, Nz, N)

    blocks = struct('src_level', {}, 'src_patch', {}, ...
                    'dst_list', {}, 'src_rows', {}, 'w_rows', {});

    for m = 1:numel(jloc_all)

        jloc = jloc_all(m);
        iloc = iloc_all(m);

        % fine 点在本层的全局网格索引
        I = dst_patch_par.i_lo + (iloc - N);
        J = dst_patch_par.j_lo + (jloc - N);

        % 目标 patch 局部线性索引
        dst_ind = sub2ind([dst_patch_par.nz + 2*N, dst_patch_par.nx + 2*N], jloc, iloc);

        % 递归展开成最终 donor atoms
        atoms = expand_interp_target(var_type, il, J, I, ...
                                     patch_catalog, patch_offset, patches, ...
                                     ghost_info_all, levels, Nx, Nz, N);

        if isempty(atoms.level0)
            error(['Failed to expand interp target for var=%s, ' ...
                   'fine patch(level=%d, patch=%d), local ghost(%d,%d).'], ...
                  var_type, dst_patch_par.level, dst_patch_par.patch_id, jloc, iloc);
        end

        % 按 (src_level, src_patch) 分组，再把这一行塞到对应 block
        key = [atoms.level0(:), atoms.patch0(:)];
        [key_u, ~, grp] = unique(key, 'rows', 'stable');

        for g = 1:size(key_u, 1)
            src_level0 = key_u(g, 1);
            src_patch0 = key_u(g, 2);

            src_row = atoms.src_ind(grp == g).';
            w_row   = atoms.w(grp == g).';

            blocks = append_restriction_block_row(blocks, src_level0, src_patch0, ...
                                                  dst_ind, src_row, w_row);
        end
    end

    vplan = finalize_restriction_blocks(blocks);
end

function atoms = expand_interp_target(var_type, il, J, I, ...
                                      patch_catalog, patch_offset, patches, ...
                                      ghost_info_all, levels, Nx, Nz, N)

    atoms = empty_atoms();

    % 先生成相邻 coarse level 上的理想 stencil
    [childJ, childI, childW] = get_interp_children(var_type, il, J, I, Nx, Nz);

    coarse_il = il - 1;
    for q = 1:numel(childW)
        atoms_q = resolve_point_on_coarser_hierarchy(var_type, coarse_il, ...
                                                     childJ(q), childI(q), childW(q), ...
                                                     patch_catalog, patch_offset, patches, ...
                                                     ghost_info_all, levels, Nx, Nz, N);
        atoms = concat_atoms(atoms, atoms_q);
    end

    atoms = merge_atoms(atoms);
end

function atoms = resolve_point_on_coarser_hierarchy(var_type, lev, J, I, base_w, ...
                                                    patch_catalog, patch_offset, patches, ...
                                                    ghost_info_all, levels, Nx, Nz, N)

    atoms = empty_atoms();

    [found, patch0, src_ind] = find_holder_patch_on_level_interp( ...
        lev, J, I, patch_catalog, patch_offset, patches, ghost_info_all, N);

    if found
        atoms.level0 = lev - 1;   % 0-based donor level
        atoms.patch0 = patch0;    % 0-based donor patch
        atoms.src_ind = src_ind;
        atoms.w = base_w;
        return;
    end

    % 到最粗层还找不到，就失败
    if lev == 1
        return;
    end

    % 当前层找不到 holder，就继续往更粗层展开
    [childJ, childI, childW] = get_interp_children(var_type, lev, J, I, Nx, Nz);

    for q = 1:numel(childW)
        atoms_q = resolve_point_on_coarser_hierarchy(var_type, lev - 1, ...
                                                     childJ(q), childI(q), base_w * childW(q), ...
                                                     patch_catalog, patch_offset, patches, ...
                                                     ghost_info_all, levels, Nx, Nz, N);
        atoms = concat_atoms(atoms, atoms_q);
    end

    atoms = merge_atoms(atoms);
end

function [found, patch0, src_ind] = find_holder_patch_on_level_interp( ...
    lev, J, I, patch_catalog, patch_offset, patches, ghost_info_all, N)

    % interp 阶段的 holder 更宽松：
    %   0 interior
    %   1 same-level copy ghost
    %   2 restriction-filled ghost
    %   3 even-coarser-donor ghost
    %
    % 因为运行顺序应为：
    %   copy -> restriction(finest->coarsest) -> interp(coarsest->finest)
    %
    % 所以在为更细层做 interp 时，coarser level 上这些点都应已可用。

    found = false;
    patch0 = -1;
    src_ind = [];

    best_rank = inf;

    for ip2 = 1:patches(lev)

        patch_par = patch_catalog(ip2 + patch_offset(lev));

        iloc = I - patch_par.i_lo + N;
        jloc = J - patch_par.j_lo + N;

        % 不在该 patch 的扩展数组范围里
        if iloc < 1 || iloc > patch_par.nx + 2*N || ...
           jloc < 1 || jloc > patch_par.nz + 2*N
            continue;
        end

        % interior
        if iloc >= N+1 && iloc <= N+patch_par.nx && ...
           jloc >= N+1 && jloc <= N+patch_par.nz
            gtype = 0;
        else
            gtype = ghost_info_all{lev}{ip2}.ghost_type(jloc, iloc);
        end

        rank = holder_rank_interp(gtype);

        if rank < best_rank
            best_rank = rank;
            found = true;
            patch0 = patch_par.patch_id;   % 0-based
            src_ind = sub2ind([patch_par.nz + 2*N, patch_par.nx + 2*N], jloc, iloc);

            if rank == 1
                return;
            end
        end
    end
end

function rank = holder_rank_interp(gtype)
    switch gtype
        case 0
            rank = 1;   % interior best
        case 1
            rank = 2;   % same-level copy ghost
        case 2
            rank = 3;   % restriction-filled ghost
        case 3
            rank = 4;   % even-coarser-donor ghost
        otherwise
            rank = inf; % -1 not accepted
    end
end

function [childJ, childI, childW] = get_interp_children(var_type, il, J, I, Nx, Nz)

    % il: 当前“细层”的 MATLAB level 编号（1-based）
    % 该点位于 level il，要向 level il-1 展开
    if il <= 1
        childJ = [];
        childI = [];
        childW = [];
        return;
    end

    fine_sz   = [Nz(il),   Nx(il)];
    coarse_sz = [Nz(il-1), Nx(il-1)];

    switch var_type
        case 'Vx'
            p = build_plan_Vx(J, I, fine_sz, coarse_sz);
        case 'Vz'
            p = build_plan_Vz(J, I, fine_sz, coarse_sz);
        case 'Txx'
            p = build_plan_Txx(J, I, fine_sz, coarse_sz);
        case 'Tzz'
            p = build_plan_Tzz(J, I, fine_sz, coarse_sz);
        case 'Txz'
            p = build_plan_Txz(J, I, fine_sz, coarse_sz);
        otherwise
            error('Unknown var_type: %s', var_type);
    end

    src = [p.src1(:), p.src2(:), p.src3(:), p.src4(:)];
    w   = [p.w1(:),   p.w2(:),   p.w3(:),   p.w4(:)];

    src = src(1, :);
    w   = w(1, :);

    valid = (w ~= 0);

    src = src(valid);
    childW = w(valid);

    [childJ, childI] = ind2sub(coarse_sz, src);
    childJ = childJ(:).';
    childI = childI(:).';
    childW = childW(:).';
end

function [nodes, weights] = interp_1d_nodes(s, nmax)

    tol = 1e-12;

    nodes = [];
    weights = [];

    % 如果刚好落在整数 coarse index 上，只取一个点
    nearest = round(s);
    if abs(s - nearest) < tol
        if nearest >= 1 && nearest <= nmax
            nodes = nearest;
            weights = 1.0;
        end
        return;
    end

    i0 = floor(s);
    i1 = i0 + 1;

    w1 = s - i0;
    w0 = 1.0 - w1;

    % 正常情况：两个 coarse 点都在 level 内
    if i0 >= 1 && i1 <= nmax
        nodes = [i0, i1];
        weights = [w0, w1];
        return;
    end

    % 边界情况：只保留 level 内的点
    % 这不是最理想的高阶边界插值，但可以避免物理边界附近越界。
    if i0 >= 1 && i0 <= nmax
        nodes = i0;
        weights = 1.0;
        return;
    end

    if i1 >= 1 && i1 <= nmax
        nodes = i1;
        weights = 1.0;
        return;
    end
end

function plan = build_plan_Vx(iz, ix, fine_sz, coarse_sz)

    nzc = coarse_sz(1); 
    nxc = coarse_sz(2);
    n   = numel(ix);

    plan.dst_ind = sub2ind(fine_sz, iz, ix);
    plan.src1 = ones(n,1);
    plan.src2 = ones(n,1);
    plan.src3 = ones(n,1);
    plan.src4 = ones(n,1);
    plan.w1 = zeros(n,1);
    plan.w2 = zeros(n,1);
    plan.w3 = zeros(n,1);
    plan.w4 = zeros(n,1);

    oddx = mod(ix,2)==1; evenx = ~oddx;
    oddz = mod(iz,2)==1; evenz = ~oddz;

    % ---- odd, odd ----
    m = oddx & oddz;
    if any(m)
        jx = (ix(m)+1)/2;
        jz = (iz(m)+1)/2;
        plan.src1(m) = sub2ind(coarse_sz, jz, jx);
        plan.w1(m) = 1.0;
    end

    % ---- even, odd ----
    m = evenx & oddz;
    if any(m)
        jx = ix(m)/2;
        jz = (iz(m)+1)/2;

        has_r = (jx+1 <= nxc);

        idx = find(m);
        ii = idx(has_r);
        jj = idx(~has_r);

        if ~isempty(ii)
            jx2 = ix(ii)/2; jz2 = (iz(ii)+1)/2;
            plan.src1(ii) = sub2ind(coarse_sz, jz2, jx2);
            plan.src2(ii) = sub2ind(coarse_sz, jz2, jx2+1);
            plan.w1(ii) = 0.5; plan.w2(ii) = 0.5;
        end
        if ~isempty(jj)
            jx2 = ix(jj)/2; jz2 = (iz(jj)+1)/2;
            plan.src1(jj) = sub2ind(coarse_sz, jz2, jx2);
            plan.w1(jj) = 1.0;
        end
    end

    % ---- odd, even ----
    m = oddx & evenz;
    if any(m)
        jx = (ix(m)+1)/2;
        jz = iz(m)/2;

        has_d = (jz+1 <= nzc);

        idx = find(m);
        ii = idx(has_d);
        jj = idx(~has_d);

        if ~isempty(ii)
            jx2 = (ix(ii)+1)/2; jz2 = iz(ii)/2;
            plan.src1(ii) = sub2ind(coarse_sz, jz2,   jx2);
            plan.src2(ii) = sub2ind(coarse_sz, jz2+1, jx2);
            plan.w1(ii) = 0.5; plan.w2(ii) = 0.5;
        end
        if ~isempty(jj)
            jx2 = (ix(jj)+1)/2; jz2 = iz(jj)/2;
            plan.src1(jj) = sub2ind(coarse_sz, jz2, jx2);
            plan.w1(jj) = 1.0;
        end
    end

    % ---- even, even ----
    m = evenx & evenz;
    if any(m)
        idx = find(m);
        jx = ix(idx)/2;
        jz = iz(idx)/2;

        has_r = (jx+1 <= nxc);
        has_d = (jz+1 <= nzc);

        ii = idx(has_r & has_d);
        if ~isempty(ii)
            jx2 = ix(ii)/2; jz2 = iz(ii)/2;
            plan.src1(ii) = sub2ind(coarse_sz, jz2,   jx2);
            plan.src2(ii) = sub2ind(coarse_sz, jz2,   jx2+1);
            plan.src3(ii) = sub2ind(coarse_sz, jz2+1, jx2);
            plan.src4(ii) = sub2ind(coarse_sz, jz2+1, jx2+1);
            plan.w1(ii) = 0.25; plan.w2(ii) = 0.25;
            plan.w3(ii) = 0.25; plan.w4(ii) = 0.25;
        end

        ii = idx(~has_r & has_d);
        if ~isempty(ii)
            jx2 = ix(ii)/2; jz2 = iz(ii)/2;
            plan.src1(ii) = sub2ind(coarse_sz, jz2,   jx2);
            plan.src2(ii) = sub2ind(coarse_sz, jz2+1, jx2);
            plan.w1(ii) = 0.5; plan.w2(ii) = 0.5;
        end

        ii = idx(has_r & ~has_d);
        if ~isempty(ii)
            jx2 = ix(ii)/2; jz2 = iz(ii)/2;
            plan.src1(ii) = sub2ind(coarse_sz, jz2, jx2);
            plan.src2(ii) = sub2ind(coarse_sz, jz2, jx2+1);
            plan.w1(ii) = 0.5; plan.w2(ii) = 0.5;
        end

        ii = idx(~has_r & ~has_d);
        if ~isempty(ii)
            jx2 = ix(ii)/2; jz2 = iz(ii)/2;
            plan.src1(ii) = sub2ind(coarse_sz, jz2, jx2);
            plan.w1(ii) = 1.0;
        end
    end
end

function plan = build_plan_Vz(iz, ix, fine_sz, coarse_sz)

    nzc = coarse_sz(1); 
    nxc = coarse_sz(2);
    n   = numel(ix);

    plan.dst_ind = sub2ind(fine_sz, iz, ix);
    plan.src1 = ones(n,1);
    plan.src2 = ones(n,1);
    plan.src3 = ones(n,1);
    plan.src4 = ones(n,1);
    plan.w1 = zeros(n,1);
    plan.w2 = zeros(n,1);
    plan.w3 = zeros(n,1);
    plan.w4 = zeros(n,1);

    oddx = mod(ix,2)==1; evenx = ~oddx;
    oddz = mod(iz,2)==1; evenz = ~oddz;

    m = oddx & oddz;
    if any(m)
        idx = find(m);
        jx = (ix(idx)+1)/2;
        jz = (iz(idx)+1)/2;

        has_l = (jx-1 >= 1);
        has_u = (jz-1 >= 1);

        ii = idx(has_l & has_u);
        if ~isempty(ii)
            jx2 = (ix(ii)+1)/2; jz2 = (iz(ii)+1)/2;
            plan.src1(ii) = sub2ind(coarse_sz, jz2-1, jx2-1);
            plan.src2(ii) = sub2ind(coarse_sz, jz2-1, jx2);
            plan.src3(ii) = sub2ind(coarse_sz, jz2,   jx2-1);
            plan.src4(ii) = sub2ind(coarse_sz, jz2,   jx2);
            plan.w1(ii)=0.0625; plan.w2(ii)=0.1875; plan.w3(ii)=0.1875; plan.w4(ii)=0.5625;
        end

        ii = idx(has_l & ~has_u);
        if ~isempty(ii)
            jx2 = (ix(ii)+1)/2; jz2 = (iz(ii)+1)/2;
            plan.src1(ii) = sub2ind(coarse_sz, jz2, jx2-1);
            plan.src2(ii) = sub2ind(coarse_sz, jz2, jx2);
            plan.w1(ii)=0.25; plan.w2(ii)=0.75;
        end

        ii = idx(~has_l & has_u);
        if ~isempty(ii)
            jx2 = (ix(ii)+1)/2; jz2 = (iz(ii)+1)/2;
            plan.src1(ii) = sub2ind(coarse_sz, jz2-1, jx2);
            plan.src2(ii) = sub2ind(coarse_sz, jz2,   jx2);
            plan.w1(ii)=0.25; plan.w2(ii)=0.75;
        end

        ii = idx(~has_l & ~has_u);
        if ~isempty(ii)
            jx2 = (ix(ii)+1)/2; jz2 = (iz(ii)+1)/2;
            plan.src1(ii) = sub2ind(coarse_sz, jz2, jx2);
            plan.w1(ii)=1.0;
        end
    end

    m = evenx & oddz;
    if any(m)
        idx = find(m);
        jx = ix(idx)/2;
        jz = (iz(idx)+1)/2;

        has_r = (jx+1 <= nxc);
        has_u = (jz-1 >= 1);

        ii = idx(has_r & has_u);
        if ~isempty(ii)
            jx2 = ix(ii)/2; jz2 = (iz(ii)+1)/2;
            plan.src1(ii) = sub2ind(coarse_sz, jz2-1, jx2);
            plan.src2(ii) = sub2ind(coarse_sz, jz2-1, jx2+1);
            plan.src3(ii) = sub2ind(coarse_sz, jz2,   jx2);
            plan.src4(ii) = sub2ind(coarse_sz, jz2,   jx2+1);
            plan.w1(ii)=0.1875; plan.w2(ii)=0.0625; plan.w3(ii)=0.5625; plan.w4(ii)=0.1875;
        end

        ii = idx(has_r & ~has_u);
        if ~isempty(ii)
            jx2 = ix(ii)/2; jz2 = (iz(ii)+1)/2;
            plan.src1(ii) = sub2ind(coarse_sz, jz2, jx2);
            plan.src2(ii) = sub2ind(coarse_sz, jz2, jx2+1);
            plan.w1(ii)=0.75; plan.w2(ii)=0.25;
        end

        ii = idx(~has_r & has_u);
        if ~isempty(ii)
            jx2 = ix(ii)/2; jz2 = (iz(ii)+1)/2;
            plan.src1(ii) = sub2ind(coarse_sz, jz2-1, jx2);
            plan.src2(ii) = sub2ind(coarse_sz, jz2,   jx2);
            plan.w1(ii)=0.25; plan.w2(ii)=0.75;
        end

        ii = idx(~has_r & ~has_u);
        if ~isempty(ii)
            jx2 = ix(ii)/2; jz2 = (iz(ii)+1)/2;
            plan.src1(ii) = sub2ind(coarse_sz, jz2, jx2);
            plan.w1(ii)=1.0;
        end
    end

    m = oddx & evenz;
    if any(m)
        idx = find(m);
        jx = (ix(idx)+1)/2;
        jz = iz(idx)/2;

        has_l = (jx-1 >= 1);
        has_d = (jz+1 <= nzc);

        ii = idx(has_l & has_d);
        if ~isempty(ii)
            jx2 = (ix(ii)+1)/2; jz2 = iz(ii)/2;
            plan.src1(ii) = sub2ind(coarse_sz, jz2,   jx2-1);
            plan.src2(ii) = sub2ind(coarse_sz, jz2,   jx2);
            plan.src3(ii) = sub2ind(coarse_sz, jz2+1, jx2-1);
            plan.src4(ii) = sub2ind(coarse_sz, jz2+1, jx2);
            plan.w1(ii)=0.1875; plan.w2(ii)=0.5625; plan.w3(ii)=0.0625; plan.w4(ii)=0.1875;
        end

        ii = idx(has_l & ~has_d);
        if ~isempty(ii)
            jx2 = (ix(ii)+1)/2; jz2 = iz(ii)/2;
            plan.src1(ii) = sub2ind(coarse_sz, jz2, jx2-1);
            plan.src2(ii) = sub2ind(coarse_sz, jz2, jx2);
            plan.w1(ii)=0.25; plan.w2(ii)=0.75;
        end

        ii = idx(~has_l & has_d);
        if ~isempty(ii)
            jx2 = (ix(ii)+1)/2; jz2 = iz(ii)/2;
            plan.src1(ii) = sub2ind(coarse_sz, jz2,   jx2);
            plan.src2(ii) = sub2ind(coarse_sz, jz2+1, jx2);
            plan.w1(ii)=0.75; plan.w2(ii)=0.25;
        end

        ii = idx(~has_l & ~has_d);
        if ~isempty(ii)
            jx2 = (ix(ii)+1)/2; jz2 = iz(ii)/2;
            plan.src1(ii) = sub2ind(coarse_sz, jz2, jx2);
            plan.w1(ii)=1.0;
        end
    end

    m = evenx & evenz;
    if any(m)
        idx = find(m);
        jx = ix(idx)/2;
        jz = iz(idx)/2;

        has_r = (jx+1 <= nxc);
        has_d = (jz+1 <= nzc);

        ii = idx(has_r & has_d);
        if ~isempty(ii)
            jx2 = ix(ii)/2; jz2 = iz(ii)/2;
            plan.src1(ii) = sub2ind(coarse_sz, jz2,   jx2);
            plan.src2(ii) = sub2ind(coarse_sz, jz2,   jx2+1);
            plan.src3(ii) = sub2ind(coarse_sz, jz2+1, jx2);
            plan.src4(ii) = sub2ind(coarse_sz, jz2+1, jx2+1);
            plan.w1(ii)=0.5625; plan.w2(ii)=0.1875; plan.w3(ii)=0.1875; plan.w4(ii)=0.0625;
        end

        ii = idx(has_r & ~has_d);
        if ~isempty(ii)
            jx2 = ix(ii)/2; jz2 = iz(ii)/2;
            plan.src1(ii) = sub2ind(coarse_sz, jz2, jx2);
            plan.src2(ii) = sub2ind(coarse_sz, jz2, jx2+1);
            plan.w1(ii)=0.75; plan.w2(ii)=0.25;
        end

        ii = idx(~has_r & has_d);
        if ~isempty(ii)
            jx2 = ix(ii)/2; jz2 = iz(ii)/2;
            plan.src1(ii) = sub2ind(coarse_sz, jz2,   jx2);
            plan.src2(ii) = sub2ind(coarse_sz, jz2+1, jx2);
            plan.w1(ii)=0.75; plan.w2(ii)=0.25;
        end

        ii = idx(~has_r & ~has_d);
        if ~isempty(ii)
            jx2 = ix(ii)/2; jz2 = iz(ii)/2;
            plan.src1(ii) = sub2ind(coarse_sz, jz2, jx2);
            plan.w1(ii)=1.0;
        end
    end
end

function plan = build_plan_Txx(iz, ix, fine_sz, coarse_sz)

    nzc = coarse_sz(1); 
    nxc = coarse_sz(2);
    n   = numel(ix);

    plan.dst_ind = sub2ind(fine_sz, iz, ix);
    plan.src1 = ones(n,1);
    plan.src2 = ones(n,1);
    plan.src3 = ones(n,1);
    plan.src4 = ones(n,1);
    plan.w1 = zeros(n,1);
    plan.w2 = zeros(n,1);
    plan.w3 = zeros(n,1);
    plan.w4 = zeros(n,1);

    oddx = mod(ix,2)==1; evenx = ~oddx;
    oddz = mod(iz,2)==1; evenz = ~oddz;

    m = oddx & oddz;
    if any(m)
        idx = find(m);
        jx = (ix(idx)+1)/2;
        jz = (iz(idx)+1)/2;
        has_l = (jx-1 >= 1);

        ii = idx(has_l);
        if ~isempty(ii)
            jx2 = (ix(ii)+1)/2; jz2 = (iz(ii)+1)/2;
            plan.src1(ii) = sub2ind(coarse_sz, jz2, jx2-1);
            plan.src2(ii) = sub2ind(coarse_sz, jz2, jx2);
            plan.w1(ii)=0.25; plan.w2(ii)=0.75;
        end

        ii = idx(~has_l);
        if ~isempty(ii)
            jx2 = (ix(ii)+1)/2; jz2 = (iz(ii)+1)/2;
            plan.src1(ii) = sub2ind(coarse_sz, jz2, jx2);
            plan.w1(ii)=1.0;
        end
    end

    m = evenx & oddz;
    if any(m)
        idx = find(m);
        jx = ix(idx)/2;
        jz = (iz(idx)+1)/2;
        has_r = (jx+1 <= nxc);

        ii = idx(has_r);
        if ~isempty(ii)
            jx2 = ix(ii)/2; jz2 = (iz(ii)+1)/2;
            plan.src1(ii) = sub2ind(coarse_sz, jz2, jx2);
            plan.src2(ii) = sub2ind(coarse_sz, jz2, jx2+1);
            plan.w1(ii)=0.75; plan.w2(ii)=0.25;
        end

        ii = idx(~has_r);
        if ~isempty(ii)
            jx2 = ix(ii)/2; jz2 = (iz(ii)+1)/2;
            plan.src1(ii) = sub2ind(coarse_sz, jz2, jx2);
            plan.w1(ii)=1.0;
        end
    end

    m = oddx & evenz;
    if any(m)
        idx = find(m);
        jx = (ix(idx)+1)/2;
        jz = iz(idx)/2;
        has_l = (jx-1 >= 1);
        has_d = (jz+1 <= nzc);

        ii = idx(has_l & has_d);
        if ~isempty(ii)
            jx2 = (ix(ii)+1)/2; jz2 = iz(ii)/2;
            plan.src1(ii) = sub2ind(coarse_sz, jz2,   jx2);
            plan.src2(ii) = sub2ind(coarse_sz, jz2,   jx2-1);
            plan.src3(ii) = sub2ind(coarse_sz, jz2+1, jx2);
            plan.src4(ii) = sub2ind(coarse_sz, jz2+1, jx2-1);
            plan.w1(ii)=0.375; plan.w2(ii)=0.125; plan.w3(ii)=0.375; plan.w4(ii)=0.125;
        end

        ii = idx(has_l & ~has_d);
        if ~isempty(ii)
            jx2 = (ix(ii)+1)/2; jz2 = iz(ii)/2;
            plan.src1(ii) = sub2ind(coarse_sz, jz2, jx2);
            plan.src2(ii) = sub2ind(coarse_sz, jz2, jx2-1);
            plan.w1(ii)=0.75; plan.w2(ii)=0.25;
        end

        ii = idx(~has_l & has_d);
        if ~isempty(ii)
            jx2 = (ix(ii)+1)/2; jz2 = iz(ii)/2;
            plan.src1(ii) = sub2ind(coarse_sz, jz2,   jx2);
            plan.src2(ii) = sub2ind(coarse_sz, jz2+1, jx2);
            plan.w1(ii)=0.5; plan.w2(ii)=0.5;
        end

        ii = idx(~has_l & ~has_d);
        if ~isempty(ii)
            jx2 = (ix(ii)+1)/2; jz2 = iz(ii)/2;
            plan.src1(ii) = sub2ind(coarse_sz, jz2, jx2);
            plan.w1(ii)=1.0;
        end
    end

    m = evenx & evenz;
    if any(m)
        idx = find(m);
        jx = ix(idx)/2;
        jz = iz(idx)/2;
        has_r = (jx+1 <= nxc);
        has_d = (jz+1 <= nzc);

        ii = idx(has_r & has_d);
        if ~isempty(ii)
            jx2 = ix(ii)/2; jz2 = iz(ii)/2;
            plan.src1(ii) = sub2ind(coarse_sz, jz2,   jx2);
            plan.src2(ii) = sub2ind(coarse_sz, jz2,   jx2+1);
            plan.src3(ii) = sub2ind(coarse_sz, jz2+1, jx2);
            plan.src4(ii) = sub2ind(coarse_sz, jz2+1, jx2+1);
            plan.w1(ii)=0.375; plan.w2(ii)=0.125; plan.w3(ii)=0.375; plan.w4(ii)=0.125;
        end

        ii = idx(has_r & ~has_d);
        if ~isempty(ii)
            jx2 = ix(ii)/2; jz2 = iz(ii)/2;
            plan.src1(ii) = sub2ind(coarse_sz, jz2, jx2);
            plan.src2(ii) = sub2ind(coarse_sz, jz2, jx2+1);
            plan.w1(ii)=0.75; plan.w2(ii)=0.25;
        end

        ii = idx(~has_r & has_d);
        if ~isempty(ii)
            jx2 = ix(ii)/2; jz2 = iz(ii)/2;
            plan.src1(ii) = sub2ind(coarse_sz, jz2,   jx2);
            plan.src2(ii) = sub2ind(coarse_sz, jz2+1, jx2);
            plan.w1(ii)=0.5; plan.w2(ii)=0.5;
        end

        ii = idx(~has_r & ~has_d);
        if ~isempty(ii)
            jx2 = ix(ii)/2; jz2 = iz(ii)/2;
            plan.src1(ii) = sub2ind(coarse_sz, jz2, jx2);
            plan.w1(ii)=1.0;
        end
    end
end

function plan = build_plan_Tzz(iz, ix, fine_sz, coarse_sz)
    plan = build_plan_Txx(iz, ix, fine_sz, coarse_sz);
end

function plan = build_plan_Txz(iz, ix, fine_sz, coarse_sz)

    nzc = coarse_sz(1); 
    nxc = coarse_sz(2);
    n   = numel(ix);

    plan.dst_ind = sub2ind(fine_sz, iz, ix);
    plan.src1 = ones(n,1);
    plan.src2 = ones(n,1);
    plan.src3 = ones(n,1);
    plan.src4 = ones(n,1);
    plan.w1 = zeros(n,1);
    plan.w2 = zeros(n,1);
    plan.w3 = zeros(n,1);
    plan.w4 = zeros(n,1);

    oddx = mod(ix,2)==1; evenx = ~oddx;
    oddz = mod(iz,2)==1; evenz = ~oddz;

    m = oddx & oddz;
    if any(m)
        idx = find(m);
        jx = (ix(idx)+1)/2;
        jz = (iz(idx)+1)/2;
        has_u = (jz-1 >= 1);

        ii = idx(has_u);
        if ~isempty(ii)
            jx2 = (ix(ii)+1)/2; jz2 = (iz(ii)+1)/2;
            plan.src1(ii) = sub2ind(coarse_sz, jz2-1, jx2);
            plan.src2(ii) = sub2ind(coarse_sz, jz2,   jx2);
            plan.w1(ii)=0.25; plan.w2(ii)=0.75;
        end

        ii = idx(~has_u);
        if ~isempty(ii)
            jx2 = (ix(ii)+1)/2; jz2 = (iz(ii)+1)/2;
            plan.src1(ii) = sub2ind(coarse_sz, jz2, jx2);
            plan.w1(ii)=1.0;
        end
    end

    m = evenx & oddz;
    if any(m)
        idx = find(m);
        jx = ix(idx)/2;
        jz = (iz(idx)+1)/2;
        has_r = (jx+1 <= nxc);
        has_u = (jz-1 >= 1);

        ii = idx(has_r & has_u);
        if ~isempty(ii)
            jx2 = ix(ii)/2; jz2 = (iz(ii)+1)/2;
            plan.src1(ii) = sub2ind(coarse_sz, jz2-1, jx2);
            plan.src2(ii) = sub2ind(coarse_sz, jz2-1, jx2+1);
            plan.src3(ii) = sub2ind(coarse_sz, jz2,   jx2);
            plan.src4(ii) = sub2ind(coarse_sz, jz2,   jx2+1);
            plan.w1(ii)=0.125; plan.w2(ii)=0.125; plan.w3(ii)=0.375; plan.w4(ii)=0.375;
        end

        ii = idx(has_r & ~has_u);
        if ~isempty(ii)
            jx2 = ix(ii)/2; jz2 = (iz(ii)+1)/2;
            plan.src1(ii) = sub2ind(coarse_sz, jz2, jx2);
            plan.src2(ii) = sub2ind(coarse_sz, jz2, jx2+1);
            plan.w1(ii)=0.5; plan.w2(ii)=0.5;
        end

        ii = idx(~has_r & has_u);
        if ~isempty(ii)
            jx2 = ix(ii)/2; jz2 = (iz(ii)+1)/2;
            plan.src1(ii) = sub2ind(coarse_sz, jz2-1, jx2);
            plan.src2(ii) = sub2ind(coarse_sz, jz2,   jx2);
            plan.w1(ii)=0.25; plan.w2(ii)=0.75;
        end

        ii = idx(~has_r & ~has_u);
        if ~isempty(ii)
            jx2 = ix(ii)/2; jz2 = (iz(ii)+1)/2;
            plan.src1(ii) = sub2ind(coarse_sz, jz2, jx2);
            plan.w1(ii)=1.0;
        end
    end

    m = oddx & evenz;
    if any(m)
        idx = find(m);
        jx = (ix(idx)+1)/2;
        jz = iz(idx)/2;
        has_d = (jz+1 <= nzc);

        ii = idx(has_d);
        if ~isempty(ii)
            jx2 = (ix(ii)+1)/2; jz2 = iz(ii)/2;
            plan.src1(ii) = sub2ind(coarse_sz, jz2,   jx2);
            plan.src2(ii) = sub2ind(coarse_sz, jz2+1, jx2);
            plan.w1(ii)=0.75; plan.w2(ii)=0.25;
        end

        ii = idx(~has_d);
        if ~isempty(ii)
            jx2 = (ix(ii)+1)/2; jz2 = iz(ii)/2;
            plan.src1(ii) = sub2ind(coarse_sz, jz2, jx2);
            plan.w1(ii)=1.0;
        end
    end

    m = evenx & evenz;
    if any(m)
        idx = find(m);
        jx = ix(idx)/2;
        jz = iz(idx)/2;
        has_r = (jx+1 <= nxc);
        has_d = (jz+1 <= nzc);

        ii = idx(has_r & has_d);
        if ~isempty(ii)
            jx2 = ix(ii)/2; jz2 = iz(ii)/2;
            plan.src1(ii) = sub2ind(coarse_sz, jz2,   jx2);
            plan.src2(ii) = sub2ind(coarse_sz, jz2,   jx2+1);
            plan.src3(ii) = sub2ind(coarse_sz, jz2+1, jx2);
            plan.src4(ii) = sub2ind(coarse_sz, jz2+1, jx2+1);
            plan.w1(ii)=0.375; plan.w2(ii)=0.375; plan.w3(ii)=0.125; plan.w4(ii)=0.125;
        end

        ii = idx(has_r & ~has_d);
        if ~isempty(ii)
            jx2 = ix(ii)/2; jz2 = iz(ii)/2;
            plan.src1(ii) = sub2ind(coarse_sz, jz2, jx2);
            plan.src2(ii) = sub2ind(coarse_sz, jz2, jx2+1);
            plan.w1(ii)=0.5; plan.w2(ii)=0.5;
        end

        ii = idx(~has_r & has_d);
        if ~isempty(ii)
            jx2 = ix(ii)/2; jz2 = iz(ii)/2;
            plan.src1(ii) = sub2ind(coarse_sz, jz2,   jx2);
            plan.src2(ii) = sub2ind(coarse_sz, jz2+1, jx2);
            plan.w1(ii)=0.75; plan.w2(ii)=0.25;
        end

        ii = idx(~has_r & ~has_d);
        if ~isempty(ii)
            jx2 = ix(ii)/2; jz2 = iz(ii)/2;
            plan.src1(ii) = sub2ind(coarse_sz, jz2, jx2);
            plan.w1(ii)=1.0;
        end
    end
end