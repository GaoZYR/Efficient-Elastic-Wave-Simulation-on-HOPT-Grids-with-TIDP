function [ricker,ricker_diff]=Ricker(fm,Nt,dt)
    t=0:dt:(Nt-1)*dt;
    tt=t-fix((1/fm)/dt)*dt;
    ricker=(1-2*pi*pi*fm*fm*tt.^2).*exp(-pi*pi*fm*fm*tt.^2);
    ricker_diff=diff(ricker);