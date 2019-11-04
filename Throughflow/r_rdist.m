function r_new = r_rdist(r, Vx2, Vt, H0, mdot, T0, P0)

global cp g R
n = length(r);
n_i = 10001;
n_step = round(n_i/n);
rh = min(r);
rt = max(r);
dr = (rt-rh)/(n_i-1);
r_i = linspace(rh,rt,n_i);
Vx_i = interp1(r,sqrt(Vx2),r_i);
Vt_i = interp1(r,Vt,r_i);
H0_i = interp1(r,H0,r_i);
dmdot = mdot/(n-1);
T = ((H0_i/cp) - (((Vt_i.^2)+(Vx_i.^2))./(2*cp)));
rho = P0*(T0^(-g/(g-1)))*(T.^(1/(g-1)))/R;
fnc = r_i.*Vx_i.*rho;
r_new = zeros(1,n);
r_new(1) = rh;
k = 1;
for j = 1:n-2
    a = k;
    b = round(k+(3*n_step));
    if b>n_i
        b = n_i;
    end
    c = floor((a+b)/2);
    while c ~= a && c ~= b
        if k==a
            ma = -dmdot;
        else
            %ma = trapz(r_i(k:a),fnc(k:a))-dmdot;
            ma = dr*sum(fnc(k:a)) - 0.5*dr*(fnc(k)+fnc(a)) - dmdot;
        end
        %mc = trapz(r_i(k:c),fnc(k:c))-dmdot;
        mc = dr*sum(fnc(k:c)) - 0.5*dr*(fnc(k)+fnc(c)) - dmdot;
        if sign(ma) == sign(mc)
            a = c;
        else
            b = c;
        end
        c = round((a+b)/2);
    end
    k = c;
    r_new(j+1) = r_i(k);
end
r_new(end) = rt;
r_new = (0.9*r) + (0.1*r_new);
end

