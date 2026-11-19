function [Txx_patches, Tzz_patches, Txz_patches] = apply_interp_plan_T(interp_plan, Txx_patches, Tzz_patches, Txz_patches)

    levels = numel(interp_plan);

    % interp must coarsest -> finest
    for il = 2:levels
        for ip = 1:numel(interp_plan{il})

            plan = interp_plan{il}{ip};
            Txx_patches{il}{ip} = apply_interp_var_plan(plan.Txx, Txx_patches{il}{ip}, Txx_patches);
            Tzz_patches{il}{ip} = apply_interp_var_plan(plan.Tzz, Tzz_patches{il}{ip}, Tzz_patches);
            Txz_patches{il}{ip} = apply_interp_var_plan(plan.Txz, Txz_patches{il}{ip}, Txz_patches);
        end
    end
end

function A = apply_interp_var_plan(vplan, A, field_patches)

    if vplan.nblock == 0
        return;
    end

    if ~isempty(vplan.dst_union)
        A(vplan.dst_union) = 0;
    end

    for k = 1:vplan.nblock

        sl = vplan.src_level(k) + 1;   % 0-based -> MATLAB 1-based
        sp = vplan.src_patch(k) + 1;   % 0-based -> MATLAB 1-based

        B = field_patches{sl}{sp};

        dst = vplan.dst_ind{k};      % [M x 1]
        src_mat = vplan.src_ind{k};  % [M x Smax]
        w_mat   = vplan.w{k};        % [M x Smax]

        vals = B(src_mat);

        vals(w_mat == 0) = 0;

        contrib = sum(vals .* w_mat, 2);

        A(dst) = A(dst) + contrib;
    end
end