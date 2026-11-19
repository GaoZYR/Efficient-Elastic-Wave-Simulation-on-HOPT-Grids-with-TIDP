function [rec_index_x, rec_index_z] = record_location_index(receiver_patch, nrec, rec_x, rec_z, N)
    for ig = 1:nrec
        patch = receiver_patch(ig);
        x_lo = patch.x_lo;
        z_lo = patch.z_lo;
        dx = patch.dx;
        dz = patch.dz;
        

        gx = rec_x(ig);
        gz = rec_z(ig);

        rec_index_x(ig) = round((gx - x_lo)/dx + N);
        rec_index_z(ig) = round((gz - z_lo)/dz + N);


    end
