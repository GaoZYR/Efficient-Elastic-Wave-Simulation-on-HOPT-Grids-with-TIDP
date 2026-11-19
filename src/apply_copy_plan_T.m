function [Txx_patches, Tzz_patches, Txz_patches] = apply_copy_plan_T(copy_plan, Txx_patches, Tzz_patches, Txz_patches)

    levels = numel(copy_plan);

    for il = 1:levels
        for ip = 1:numel(copy_plan{il})

            plan = copy_plan{il}{ip};
            if plan.nblock == 0
                continue;
            end

            for k = 1:plan.nblock
                sp0 = plan.src_patch(k);   % 0-based
                sp = sp0 + 1;            

                dst_ind = plan.dst_ind{k};
                src_ind = plan.src_ind{k};

                A = Txx_patches{il}{ip};
                B = Txx_patches{il}{sp};
                A(dst_ind) = B(src_ind);
                Txx_patches{il}{ip} = A;

                A = Tzz_patches{il}{ip};
                B = Tzz_patches{il}{sp};
                A(dst_ind) = B(src_ind);
                Tzz_patches{il}{ip} = A;

                A = Txz_patches{il}{ip};
                B = Txz_patches{il}{sp};
                A(dst_ind) = B(src_ind);
                Txz_patches{il}{ip} = A;
            end
        end
    end
end