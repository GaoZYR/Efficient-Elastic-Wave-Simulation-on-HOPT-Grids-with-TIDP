%%
function [source_patch, receiver_patch, rec_bilinear_plan, s, ricker, ricker_diff] = ...
    build_source_receivers(patch_catalog, src_x, src_z, rec_x, rec_z, fm, Nt, dt)


    % Find the patch that source in
    source_patch = locate_point_in_highest_patch(patch_catalog, src_x, src_z);
    assert(~isempty(source_patch), 'Source is not inside any patch.');
    source_patch.x_src = src_x; source_patch.z_src = src_z;
    % Find the patch that receivers in
    [receiver_patch, rec_bilinear_plan] = locate_receivers_in_highest_patches(patch_catalog, rec_x, rec_z);
    for ir = 1:numel(receiver_patch)
        assert(~isempty(receiver_patch(ir).level), 'Receiver %d is not inside any patch.', ir);
    end
    
    % Calculate Source
    [ricker, ricker_diff] = Ricker(fm, Nt, dt);
    % ricker = FTDT(ricker);
    for iz = 1:source_patch.nz
        zc = source_patch.z_lo + 0.5 * source_patch.dz + (iz - 1) * source_patch.dz;
        for ix = 1:source_patch.nx
            xc = source_patch.x_lo + 0.5 * source_patch.dx + (ix - 1) * source_patch.dx;
            s(iz, ix) = exp( -((xc - src_x)^2)/(source_patch.dx^2) - ((zc - src_z)^2)/(source_patch.dz^2) );
        end
    end
    s_sum = sum(sum(s));
    s = s./s_sum;
end
%% This function is used to find the patches that receivers in
function [receiver_info, rec_bilinear_plan] = locate_receivers_in_highest_patches(patch_catalog, rec_x, rec_z)

    nrec = numel(rec_x);

    % Build receivers patches
    receiver_info = repmat(struct('level', [], 'patch_id', [], ...
                                  'x_lo', [], 'x_hi', [], 'z_lo', [], 'z_hi', [], ...
                                  'nx', [], 'nz', [], 'dx', [], 'dz', [], ...
                                  'i_lo', [], 'i_hi', [], 'j_lo', [], 'j_hi', []), 1, nrec);

    for ir = 1:nrec
        receiver_info(ir) = locate_point_in_highest_patch(patch_catalog, rec_x(ir), rec_z(ir));
        receiver_patch = receiver_info(ir);

        x_lo = receiver_patch.x_lo; x_hi = receiver_patch.x_hi;
        z_lo = receiver_patch.z_lo; z_hi = receiver_patch.z_hi;

        % Judge the nature of receivers
        dx = receiver_patch.dx; dz = receiver_patch.dz;
        gx = rec_x(ir); gz = rec_z(ir);

        ngx = (gx - x_lo - 0.5*dx)/dx + 1;
        ngz = (gz - z_lo - 0.5*dz)/dz + 1;

        ngx_up = ceil(ngx); ngx_low = floor(ngx);
        ngz_up = ceil(ngz); ngz_low = floor(ngz);

        rec_bilinear_plan(ir).w(1) = 0; rec_bilinear_plan(ir).x(1) = 0; rec_bilinear_plan(ir).z(1) = 0;
        rec_bilinear_plan(ir).w(2) = 0; rec_bilinear_plan(ir).x(2) = 0; rec_bilinear_plan(ir).z(2) = 0;
        rec_bilinear_plan(ir).w(3) = 0; rec_bilinear_plan(ir).x(3) = 0; rec_bilinear_plan(ir).z(3) = 0;
        rec_bilinear_plan(ir).w(4) = 0; rec_bilinear_plan(ir).x(4) = 0; rec_bilinear_plan(ir).z(4) = 0;
        rec_bilinear_plan(ir).w(5) = 0; rec_bilinear_plan(ir).x(5) = 0; rec_bilinear_plan(ir).z(5) = 0;

        tol = 1e-10;
        if (abs(ngx_up - ngx_low) < tol) && (abs(ngz_up - ngz_low) < tol)
            rec_bilinear_plan(ir).w(1) = 1;
            rec_bilinear_plan(ir).x(1) = ngx_low; rec_bilinear_plan(ir).z(1) = ngz_low;

        elseif (abs(ngx_up - ngx_low) < tol) && (abs(ngz_up - ngz_low) >= tol)
            % rz(2) = (gz - z_lo)/(z_hi - z_lo); rz(3) = (z_hi - gz)/(z_hi - z_lo);
            rz(2) = ngz - ngz_low; rz(3) = ngz_up - ngz;

            rec_bilinear_plan(ir).w(2) = rz(3);
            rec_bilinear_plan(ir).x(2) = ngx_low; rec_bilinear_plan(ir).z(2) = ngz_low;

            rec_bilinear_plan(ir).w(3) = rz(2);
            rec_bilinear_plan(ir).x(3) = ngx_low; rec_bilinear_plan(ir).z(3) = ngz_up;

        elseif (abs(ngx_up - ngx_low) >= tol) && (abs(ngz_up - ngz_low) < tol)
            % rx(2) = (gx - x_lo)/(x_hi - x_lo); rx(3) = (x_hi - gx)/(x_hi - x_lo);
            rx(2) = ngx - ngx_low; rx(3) = ngx_up - ngx;

            rec_bilinear_plan(ir).w(2) = rx(3);
            rec_bilinear_plan(ir).x(2) = ngx_low; rec_bilinear_plan(ir).z(2) = ngz_low;

            rec_bilinear_plan(ir).w(4) = rx(2);
            rec_bilinear_plan(ir).x(4) = ngx_up; rec_bilinear_plan(ir).z(4) = ngz_low;

        elseif (abs(ngx_up - ngx_low) >= tol) && (abs(ngz_up - ngz_low) >= tol)
            % rx(2) = (gx - x_lo)/(x_hi - x_lo); rx(3) = (x_hi - gx)/(x_hi - x_lo);
            % rz(2) = (gz - z_lo)/(z_hi - z_lo); rz(3) = (z_hi - gz)/(z_hi - z_lo);
            rx(2) = ngx - ngx_low; rx(3) = ngx_up - ngx;
            rz(2) = ngz - ngz_low; rz(3) = ngz_up - ngz;

            rec_bilinear_plan(ir).w(2) = rx(3)*rz(3);
            rec_bilinear_plan(ir).x(2) = ngx_low; rec_bilinear_plan(ir).z(2) = ngz_low;

            rec_bilinear_plan(ir).w(3) = rx(3)*rz(2);
            rec_bilinear_plan(ir).x(3) = ngx_low; rec_bilinear_plan(ir).z(3) = ngz_up;

            rec_bilinear_plan(ir).w(4) = rx(2)*rz(3);
            rec_bilinear_plan(ir).x(4) = ngx_up; rec_bilinear_plan(ir).z(4) = ngz_low;

            rec_bilinear_plan(ir).w(5) = rx(2)*rz(2);
            rec_bilinear_plan(ir).x(5) = ngx_up; rec_bilinear_plan(ir).z(5) = ngz_up;
        end
    end
end

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