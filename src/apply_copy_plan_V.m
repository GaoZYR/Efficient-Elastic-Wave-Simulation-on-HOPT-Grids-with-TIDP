function [Vx_patches, Vz_patches] = apply_copy_plan_V(copy_plan, Vx_patches, Vz_patches)

    levels = numel(copy_plan);

    for il = 1:levels
        for ip = 1:numel(copy_plan{il})

            plan = copy_plan{il}{ip};
            if plan.nblock == 0
                continue;
            end

            for k = 1:plan.nblock
                sp0 = plan.src_patch(k);  
                sp = sp0 + 1;             

                dst_ind = plan.dst_ind{k};
                src_ind = plan.src_ind{k};

                A = Vx_patches{il}{ip};
                B = Vx_patches{il}{sp};
                A(dst_ind) = B(src_ind);
                Vx_patches{il}{ip} = A;

                A = Vz_patches{il}{ip};
                B = Vz_patches{il}{sp};
                A(dst_ind) = B(src_ind);
                Vz_patches{il}{ip} = A;
            end
        end
    end
end