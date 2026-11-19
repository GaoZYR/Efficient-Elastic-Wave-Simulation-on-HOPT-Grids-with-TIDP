function [ddx, ddx_half, dkx, dkx_half, dalphax, dalphax_half, ...
    ddz, ddz_half, dkz, dkz_half, dalphaz, dalphaz_half, ...
    a_x, a_x_half, b_x, b_x_half, a_z, a_z_half, b_z, b_z_half] = cpml_coeff(nLx, nLz, fm, vmax, Nx, Nz, dx, dt, free_surface)
    R = 0.00000001;
    n1 = 2; n2 = 2; n3 = 1;

    kmax = 1;
    alphamax = pi*fm;
    % vmax = max(max(vp));
    dmax_x = - (n1 + 1)*vmax*log(R)/(2*nLx*dx);
    dmax_z = - (n1 + 1)*vmax*log(R)/(2*nLz*dx);
    
    ddx = zeros(1, Nx); dkx = zeros(1, Nx) + 1; dalphax = zeros(1, Nx);
    ddx_half = zeros(1, Nx);dkx_half = zeros(1, Nx) + 1; dalphax_half = zeros(1, Nx);
    ddz = zeros(1, Nz); dkz = zeros(1, Nz) + 1; dalphaz = zeros(1, Nz);
    ddz_half = zeros(1, Nz); dkz_half = zeros(1, Nz) + 1; dalphaz_half = zeros(1, Nz);
    
    a_x = zeros(1, Nx); a_x_half = zeros(1, Nx);
    b_x = zeros(1, Nx); b_x_half = zeros(1, Nx);
    a_z = zeros(1, Nz); a_z_half = zeros(1, Nz);
    b_z = zeros(1, Nz); b_z_half = zeros(1, Nz);

    %% x-direction
    % full grid
    for i = 1:Nx
        if i - 1 <= nLx
            ddx(i) = dmax_x*((nLx - i + 1)/nLx)^n1;
            dkx(i) = 1 + (kmax - 1)*((nLx - i + 1)/nLx)^n2;
            dalphax(i) = alphamax*(1 - (nLx - i + 1)/nLx)^n3;      
        elseif i - 1 > nLx && i - 1 <= fix(Nx/2)
            ddx(i) = 0;
            dkx(i) = 1; 
            dalphax(i) = 0;
        else
            ddx(i) = ddx(Nx - i + 1);
            dkx(i) = dkx(Nx - i + 1);
            dalphax(i) = dalphax(Nx - i + 1);
        end
        if dalphax(i) < 0
            dalphax(i) = 0;
        end
        b_x(i) = exp( - (ddx(i)/dkx(i) + dalphax(i))*dt);
        if abs(ddx(i)) > 1e-6
            a_x(i) = ddx(i)*(b_x(i) - 1)/(dkx(i)*(ddx(i) + dkx(i)*dalphax(i)));
        end
    end
    % half grid
    for i = 1:Nx
        if i - 1 - 0.5 <= nLx
            ddx_half(i) = dmax_x*((nLx - i + 1 + 0.5)/nLx)^n1;
            dkx_half(i) = 1 + (kmax - 1)*((nLx - i + 1 + 0.5)/nLx)^n2;
            dalphax_half(i) = alphamax*(1 - (nLx - i + 1 + 0.5)/nLx)^n3;      
        elseif i - 1 - 0.5 > nLx && i - 1 - 0.5 <= fix(Nx/2)
            ddx_half(i) = 0;
            dkx_half(i) = 1; 
            dalphax_half(i) = 0;
        else
            ddx_half(i) = ddx_half(Nx - i + 1);
            dkx_half(i) = dkx_half(Nx - i + 1);
            dalphax_half(i) = dalphax_half(Nx - i + 1);
        end
        if dalphax_half(i) < 0
            dalphax_half(i) = 0;
        end
        b_x_half(i) = exp( - (ddx_half(i)/dkx_half(i) + dalphax_half(i))*dt);
        if abs(ddx_half(i)) > 1e-6
            a_x_half(i) = ddx_half(i)*(b_x_half(i) - 1)/(dkx_half(i)*(ddx_half(i) + dkx_half(i)*dalphax_half(i)));
        end
    end

    %% z-direction
    % full grid
    for i = 1:Nz
        if i - 1 <= nLz
            ddz(i) = dmax_z*((nLz - i + 1)/nLz)^n1;
            dkz(i) = 1 + (kmax - 1)*((nLz - i + 1)/nLz)^n2;
            dalphaz(i) = alphamax*(1 - (nLz - i + 1)/nLz)^n3;
        elseif i - 1 > nLz && i - 1 <= fix(Nz/2)
            ddz(i) = 0;
            dkz(i) = 1; 
            dalphaz(i) = 0;
        else
            ddz(i) = ddz(Nz - i + 1);
            dkz(i) = dkz(Nz - i + 1);
            dalphaz(i) = dalphaz(Nz - i + 1);
        end
        if dalphaz(i) < 0
            dalphaz(i) = 0;
        end
        b_z(i) = exp( - (ddz(i)/dkz(i) + dalphaz(i))*dt);
        if abs(ddz(i)) > 1e-6
            a_z(i) = ddz(i)*(b_z(i) - 1)/(dkz(i)*(ddz(i) + dkz(i)*dalphaz(i)));
        end
    end
    % half grid
    for i = 1:Nz
        if i - 1 - 0.5 <= nLz
            ddz_half(i) = dmax_z*((nLz - i + 1 + 0.5)/nLz)^n1;
            dkz_half(i) = 1 + (kmax - 1)*((nLz - i + 1 + 0.5)/nLz)^n2;
            dalphaz_half(i) = alphamax*(1 - (nLz - i + 1 + 0.5)/nLz)^n3;
        elseif i - 1 - 0.5 > nLz && i - 1 - 0.5 <= fix(Nz/2)
            ddz_half(i) = 0;
            dkz_half(i) = 1; 
            dalphaz_half(i) = 0;
        else
            ddz_half(i) = ddz_half(Nz - i + 1);
            dkz_half(i) = dkz_half(Nz - i + 1);
            dalphaz_half(i) = dalphaz_half(Nz - i + 1);
        end
        if dalphaz_half(i) < 0
            dalphaz_half(i) = 0;
        end
        b_z_half(i) = exp( - (ddz_half(i)/dkz_half(i) + dalphaz_half(i))*dt);
        if abs(ddz_half(i)) > 1e-6
            a_z_half(i) = ddz_half(i)*(b_z_half(i) - 1)/(dkz_half(i)*(ddz_half(i) + dkz_half(i)*dalphaz_half(i)));
        end
    end
%% For source on the top
if free_surface == 1
    for i = 1:Nz
        if i <= nLz
            ddz(i) = 0;
            ddz_half(i) = 0;
            dkz(i) = 1;
            dkz_half(i) = 1;
            dalphaz(i) = 0;
            dalphaz_half(i) = 0;
            a_z(i) = 0;
            a_z_half(i) = 0;
            b_z(i) = 1;
            b_z_half(i) = 1;
        end
    end
end