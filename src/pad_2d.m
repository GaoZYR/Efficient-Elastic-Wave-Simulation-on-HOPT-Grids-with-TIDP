function [v]=pad_2d(npd_top,npd_blr,v0)
    [Nz,Nx]=size(v0);
    v=zeros(Nz+npd_top+npd_blr,Nx+2*npd_blr);
    v(npd_top+1:npd_top+Nz,npd_blr+1:npd_blr+Nx)=v0;
    v(1:npd_top,1:npd_blr)=v0(1,1);
    v(1:npd_top,npd_blr+Nx+1:end)=v0(1,end);
    v(npd_top+Nz+1:end,1:npd_blr)=v0(end,1);
    v(npd_top+Nz+1:end,npd_blr+Nx+1:end)=v0(end,end);
    for iz=1:npd_top
        v(iz,npd_blr+1:npd_blr+Nx)=v0(1,:);
    end
    for iz=npd_top+Nz+1:Nz+npd_top+npd_blr
        v(iz,npd_blr+1:npd_blr+Nx)=v0(end,:);
    end
    for ix=1:npd_blr
        v(npd_top+1:npd_top+Nz,ix)=v0(:,1);
    end
    for ix=npd_blr+Nx+1:Nx+npd_top+npd_blr
        v(npd_top+1:npd_top+Nz,ix)=v0(:,end);
    end
end