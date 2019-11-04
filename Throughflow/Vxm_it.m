function Vxm = Vxm_it(r,dVx2,Vt,H0,mdot,T0,P0)
%Vxm_it Sign change iteration for Vxm

global g R T01rm cp
a = sqrt(abs(min(dVx2)));
b = sqrt(g*R*T01rm); % max speed is ~sonic
e = 1;
i = 0;
while abs(e) > 0.001
    Vxa = (a^2) + dVx2;
    Ta = ((H0/cp) - (((Vt.^2)+Vxa)./(2*cp)));
    ra = P0*(T0^(-g/(g-1)))*(Ta.^(1/(g-1)))/R;
    Vxaa = sqrt(Vxa);
    Vxaa = real(Vxaa) - imag(Vxaa);
    ma = trapz(r,r.*ra.*Vxaa) - mdot;
    c = (a+b)/2;
    Vxc = (c^2) + dVx2;
    Tc = ((H0/cp) - (((Vt.^2)+Vxc)./(2*cp)));
    rc = P0*(T0^(-g/(g-1)))*(Tc.^(1/(g-1)))/R;
    Vxcc = sqrt(Vxc);
    Vxcc = real(Vxcc) - imag(Vxcc);
    e = trapz(r,r.*rc.*Vxcc) - mdot;
    if sign(ma) == sign(e)
        a = c;
    else
        b = c;
    end
    if i == 1000
        a = sqrt(abs(min(dVx2)));
        Vxa = (a^2) + dVx2;
        Ta = ((H0/cp) - (((Vt.^2)+Vxa)./(2*cp)));
        ra = P0*(T0^(-g/(g-1)))*(Ta.^(1/(g-1)))/R;
        Vxaa = sqrt(Vxa);
        Vxaa = real(Vxaa) - imag(Vxaa);
        ma = trapz(r,r.*ra.*Vxaa) - mdot
        if ma < 0.001
            Vxm = a;
            return
        else
            error('not converged')
        end
    end
    i = i+1; 
end
Vxm = c;
end

