function [r,Vx,Vt,H0] = Through_Flow_2(dr)
%THROUGH_FLOW_2 Through flow solver for a set of anulus areas
%   See Graham's notes for details 

global phi rm Kt w T01rm P01rm p q cp g R CRT
n_ST = 101;
if mod(n_ST,2) == 0
    nm = n_ST/2;
else
    nm = (n_ST+1)/2;
end
nTF_it = 30;
n_H01_it = 1;
if p == 0
    Vtm = Kt(1);
    Vtmr = Vtm - (w(1)*rm);
else
    Vtm = (Kt(1)*(rm^p))-(Kt(2)/rm);
    Vtmr = Vtm - (w(1)*rm);
end
T01 = T01rm + (((Vtm^2)-(Vtmr^2))/(2*cp));
P01 = P01rm*((T01/T01rm)^(g/g-1));
r1 = linspace(rm-dr(1),rm+dr(1),n_ST);
r2 = linspace(rm-dr(2),rm+dr(2),n_ST);
r3 = linspace(rm-dr(3),rm+dr(3),n_ST);
Vx1m = (phi*rm*(w(2)-w(1)))^2;
X = [r1(1)^2 r1(1) 1; rm^2 rm 1; r1(end)^2 r1(end) 1;];
Y = [(1+q)*cp*T01; cp*T01; (1-q)*cp*T01;];
A = X\Y;
H01 = ((A(1)*(r1.^2))+(A(2)*r1)+A(3));
if p == 0
    Vt1 = Kt(1)*(r1./r1);
%elseif p == -1 && q == 0
%    Vt1 = (Kt(1)*(r1.^p))-(Kt(2)./r1)+((w(2)*(rm^2))./r1);
else
    Vt1 = (Kt(1)*(r1.^p))-(Kt(2)./r1);%+(w(2)*rm);
end
figure(50)
hold on
plot(Vt1,r1)
for k = 1:n_H01_it  
    %disp('bing')
    r1 = linspace(rm-dr(1),rm+dr(1),n_ST);
    r2 = linspace(rm-dr(2),rm+dr(2),n_ST);
    r3 = linspace(rm-dr(3),rm+dr(3),n_ST);
    for i = 1:nTF_it
        r1avg2 = ((r1 + r1(nm))/2).^2;
        if k == 1
            if p == 0
                Vt1 = Kt(1)*(r1./r1);
            elseif p == -1 && q == 0
                Vt1 = (Kt(1)*(r1.^p))-(Kt(2)./r1);
                figure(10)
                plot(Vt1)
            else
                Vt1 = (Kt(1)*(r1.^p))-(Kt(2)./r1);
            end
            H01 = ((A(1)*(r1.^2))+(A(2)*r1)+A(3));
            Vx12 = Vx1m + (2*((H01-H01(nm)) - ((((r1.*Vt1).^2)-((r1(nm)*Vt1(nm))^2))./(2*r1avg2))));
            for j = 1:n_ST
                if sign(Vx12(j))==-1
                    disp(Vx12)
                    error('negative velocity')
                end
            end
            T = ((H01/cp) - (((Vt1.^2)+Vx12)./(2*cp)));
            rho = P01*(T01^(-g/(g-1)))*(T.^(1/(g-1)))/R;
            mdot = trapz(r1,r1.*rho.*sqrt(Vx12));
        elseif k~=1 && i == 1
            H01m = cp*T01;
            dH01 = 0.5*(((Vx1m.^2)-(Vx1m(nm)^2)) + ((((r1.*Vt1).^2)-((r1(nm)*Vt1(nm))^2))./(r1avg2)));
            H01 = H01m + dH01;
            T = ((H01/cp) - (((Vt1.^2)+(Vx1m.^2))./(2*cp)));
            rho = P01*(T01^(-g/(g-1)))*(T.^(1/(g-1)))/R;
            mdot = trapz(r1,r1.*rho.*Vx1m);
            Vx12 = Vx1m.^2;
        else
            T = ((H01/cp) - (((Vt1.^2)+Vx12)./(2*cp)));
            rho = P01*(T01^(-g/(g-1)))*(T.^(1/(g-1)))/R;
            mdot = trapz(r1,r1.*rho.*sqrt(Vx12));
        end
        %Inlet plane 
        
        %March along streamlines
        if p == 0
            Vt2 = Kt(2)*(r2./r2);
        elseif p == -1 && q == 0
            Vt2 = (Kt(1)*(r2.^p))+(Kt(2)./r2);
        else
            Vt2 = (Kt(1)*(r2.^p))+(Kt(2)./r2);
        end
        H02 = (w(1)*((r2.*Vt2)-(r1.*Vt1))) + H01;  
    
        if p == 0
            Vt3 = Kt(2)*(r3./r3);
        elseif p == -1 && q == 0
            Vt3 = (Kt(1)*(r3.^p))-(Kt(2)./r3);
        else
            Vt3 = (Kt(1)*(r3.^p))-(Kt(2)./r3);
        end
        H03 = (w(2)*((r3.*Vt3)-(r2.*Vt2))) + H02;
        %Mid stage plane
        r2avg2 = ((r2 + r2(nm))/2).^2;
        dVx22 = (2*((H02-H02(nm)) - ((((r2.*Vt2).^2)-((r2(nm)*Vt2(nm))^2))./(2*r2avg2))));
        Vxm2 = Vxm_it(r2,dVx22,Vt2,H02,mdot,T01,P01);
        Vx22 = dVx22 + (Vxm2^2);
        %End stage plane
        r3avg2 = ((r3 + r3(nm))/2).^2;
        dVx32 = (2*((H03-H03(nm)) - ((((r3.*Vt3).^2)-((r3(nm)*Vt3(nm))^2))./(2*r3avg2))));
        Vxm3 = Vxm_it(r3,dVx32,Vt3,H03,mdot,T01,P01);
        Vx32 = dVx32 + (Vxm3^2);

        %Redistribute streamlines
        r1_new = r_rdist(r1, Vx12, Vt1, H01, mdot, T01, P01);
        r2 = r_rdist(r2, Vx22, Vt2, H02, mdot, T01, P01);
        r3 = r_rdist(r3, Vx32, Vt3, H03, mdot, T01, P01); 
        H01 = interp1(r1,H01,r1_new);
        Vt1 = interp1(r1,Vt1,r1_new);
        r1 = r1_new;
    end
    %[Vt1,Vx1m] = rpt_stg_2(r3,Vt3,sqrt(Vx32),r1,Vt1,sqrt(Vx12),n_ST);
    %figure(100+CRT)
    %hold on
    %plot(r1,H01/cp)
end
r = [r1;r2;r3];
Vx = [sqrt(Vx12);sqrt(Vx22);sqrt(Vx32);];
Vt = [Vt1;Vt2;Vt3;];
H0 = [H01;H02;H03];
end

