function [M1,M2,M3,M4,D1,D2,D3,D4] = add_mpml( ...
    M1,M2,M3,M4,D1,D2,D3,D4,ix,iz,nx,nz,N,index, ...
    ax,bx,az,bz,kx,kz,axh,bxh,azh,bzh,kxh,kzh,dt,eta,kind)
% Cross-direction damping extension of the original CPML recurrence.
% eta=0 reproduces add_cpml. Separate stress/velocity stagger positions matter.
    if eta==0
        [M1,M2,M3,M4,D1,D2,D3,D4]=add_cpml( ...
            M1,M2,M3,M4,D1,D2,D3,D4,ix,iz,nx,nz,N,index, ...
            ax,bx,az,bz,kx,kz,axh,bxh,azh,bzh,kxh,kzh);
        return
    end
    xx=ix:ix+nx+2*N-1; zz=iz:iz+nz+2*N-1;
    [dx,alx]=recover(ax(xx),bx(xx),kx(xx),dt);
    [dxh,alxh]=recover(axh(xx),bxh(xx),kxh(xx),dt);
    [dz,alz]=recover(az(zz)',bz(zz)',kz(zz)',dt);
    [dzh,alzh]=recover(azh(zz)',bzh(zz)',kzh(zz)',dt);
    if strcmp(kind,'stress')
        [a1,b1]=coeff(dxh+eta*dz,alxh,kxh(xx),dt);
        [a2,b2]=coeff(dz+eta*dxh,alz,kz(zz)',dt);
        [a3,b3]=coeff(dx+eta*dzh,alx,kx(xx),dt);
        [a4,b4]=coeff(dzh+eta*dx,alzh,kzh(zz)',dt);
    elseif strcmp(kind,'velocity')
        [a1,b1]=coeff(dxh+eta*dzh,alxh,kxh(xx),dt);
        [a2,b2]=coeff(dz+eta*dx,alz,kz(zz)',dt);
        [a3,b3]=coeff(dx+eta*dz,alx,kx(xx),dt);
        [a4,b4]=coeff(dzh+eta*dxh,alzh,kzh(zz)',dt);
    else
        error('kind must be stress or velocity');
    end
    M1{index}=b1.*M1{index}+a1.*D1;
    M2{index}=b2.*M2{index}+a2.*D2;
    M3{index}=b3.*M3{index}+a3.*D3;
    M4{index}=b4.*M4{index}+a4.*D4;
    D1=D1./kxh(xx)+M1{index};
    D2=D2./kz(zz)'+M2{index};
    D3=D3./kx(xx)+M3{index};
    D4=D4./kzh(zz)'+M4{index};
end

function [d,alpha]=recover(a,b,k,dt)
    rate=-log(b)/dt;
    d=zeros(size(a)); mask=(b~=1);
    d(mask)=a(mask).*k(mask).^2.*rate(mask)./(b(mask)-1);
    alpha=max(0,rate-d./k);
end

function [a,b]=coeff(d,alpha,k,dt)
    b=exp(-(d./k+alpha)*dt);
    a=zeros(size(d)); den=k.*(d+k.*alpha);mask=(d>1e-6);
    a(mask)=d(mask).*(b(mask)-1)./den(mask);
end
