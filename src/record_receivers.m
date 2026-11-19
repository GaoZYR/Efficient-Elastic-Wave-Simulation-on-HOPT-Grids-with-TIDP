function [record_Vx] = record_receivers(Vx_patches, record_Vx, receiver_patch, rec_bilinear_plan, nrec, N, it)

    for ir = 1:nrec

        rec = receiver_patch(ir);
        rec_level = rec.level;
        rec_patch = rec.patch_id;

        % Vx0 = Vx_patches{rec_level + 1}{rec_patch + 1};
        % Vx = Vx0(N + 1:rec.nz + N, N + 1:rec.nx + N);

        Vx = Vx_patches{rec_level + 1}{rec_patch + 1};

        w = rec_bilinear_plan(ir).w; x = rec_bilinear_plan(ir).x + N; z = rec_bilinear_plan(ir).z + N;

        for ip = 1:5
            if x(ip) > 0 && z(ip) > 0 && w(ip) > 0
                record_Vx(it, ir) = record_Vx(it, ir) + w(ip)*Vx(z(ip), x(ip));
            end
        end

    end
end