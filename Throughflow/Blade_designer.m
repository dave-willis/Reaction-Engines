clear all
close all

%Ideal midspan conditions
global phi psi L 
phi = 0.3;
psi = 0.8;
L = 0.45;

%Gad constants
global g R cp
g = 1.6625;
R = 2067;
cp = 5187;

%Design conditions
global T01rm P01rm M3r
T01rm = 950;
P01rm = 14500000;
M3r = 0.3;

%Geometry
global rm H_init Cx
Cx = 0.0073;
AR = 1.7;
HTR = 0.966;
H_init = AR*Cx;
rm = H_init*(1+HTR)/(2*(1-HTR));
%H_init = 2*rm*(1-HTR)/(1+HTR);


%Vt distribution
global p q
p = -1;
q = 0;
% Vt = A*(r^p) -+ B/r (if p = 0, B = 0 for constant tangental velocity)

%Generate areas, rotational speeds and vt constants
global w Kt CRT
CRT = 0;
[w, Kt, drn] = Flow_guess();
%Begin throughflow iteration
nA_it = 1;
dr_histn = zeros(3,nA_it);
Vxm_des = phi*(w(2)-w(1))*rm;
for i = 1:nA_it
    %Loop for area iteration
    disp(i)
    dr_histn(:,i) = drn; %Update dr history to observe convergence
    [rn,Vxn,Vtn,H0n] = Through_Flow_2(drn);  %Do through flow solution for a given dr
    %Update dr
    rn_n = zeros(size(rn));
    Vxn_n = zeros(size(rn));
    Vtn_n = zeros(size(rn));
    H0n_n = zeros(size(rn));
    if i == nA_it
        
        for j = 1:3
            rn_n(j,:) = linspace(rm-drn(j),rm+drn(j),length(rn_n(j,:)));
            Vxn_n(j,:) = interp1(rn(j,:),Vxn(j,:),rn_n(j,:));
            Vtn_n(j,:) = interp1(rn(j,:),Vtn(j,:),rn_n(j,:));
            H0n_n(j,:) = interp1(rn(j,:),H0n(j,:),rn_n(j,:));
        end
    else
        for j = 2:3
            Vxm = interp1(rn(j,:),Vxn(j,:),rm);
            drn(j) = (0.9*drn(j)) + (0.1*drn(j)*Vxm/Vxm_des);
        end
    end
end
Kt_array = Kt;
w_array = w;
CRT_array = CRT;
CRT  = 1;
[w, Kt, drc] = Flow_guess();
%Begin throughflow iteration
dr_histc = zeros(3,nA_it);
Vxm_des = phi*(w(2)-w(1))*rm;
for i = 1:nA_it
    %Loop for area iteration
    disp(i)
    dr_histc(:,i) = drc; %Update dr history to observe convergence
    [rc,Vxc,Vtc,H0c] = Through_Flow_2(drc);  %Do through flow solution for a given dr
    %Update dr
    rc_n = zeros(size(rc));
    Vxc_n = zeros(size(rn));
    Vtc_n = zeros(size(rn));
    H0c_n = zeros(size(rn));
    if i == nA_it
        
        for j = 1:3
            rc_n(j,:) = linspace(rm-drc(j),rm+drc(j),length(rc_n(j,:)));
            Vxc_n(j,:) = interp1(rc(j,:),Vxc(j,:),rc_n(j,:));
            Vtc_n(j,:) = interp1(rc(j,:),Vtc(j,:),rc_n(j,:));
            H0c_n(j,:) = interp1(rc(j,:),H0c(j,:),rc_n(j,:));
        end
    else
        for j = 2:3
            Vxm = interp1(rc(j,:),Vxc(j,:),rm);
            drc(j) = (0.9*drc(j)) + (0.1*drc(j)*Vxm/Vxm_des);
        end
    end
end
w_array = [w; w_array];
Kt_array = [Kt; Kt_array];
CRT_aaray = [CRT; CRT_array];
Vx = [Vxc_n; Vxn_n];
Vt = [Vtc_n; Vtn_n];
H0 = [H0c_n; H0n_n];
r = [rc_n; rn_n];
Vx2 = [Vxc; Vxn];
Vt2 = [Vtc; Vtn];
H02 = [H0c; H0n];
r2 = [rc; rn];
X = Flow_Angles(r2,Vx2,Vt2,w_array);
h = H0 - 0.5*((Vx.^2)+(Vt.^2));
Lam = zeros(2,length(h(1,:)));
Lam(1,:) = (h(2,:)-h(3,:))./(h(1,:)-h(3,:));
Lam(2,:) = (h(5,:)-h(6,:))./(h(4,:)-h(6,:));
%Graphic_6_2(r,Vx,'Vx')
%Graphic_6_2(r,Vt,'Vt')
%Graphic_6_2(r,H0,'H0')
%Graphic_8_2(r,X,'Blade angle')
global Z_ideal 
Z_ideal = 0.8;
%[Zx,pcrx,nbx,Zm,pcrm,nbm] = PitchCalc_2(X,r);
% Graphic_8_2(r,Zx,'Zweifel coefficient')
% Graphic_8_2(r,Zm,'Zweifel coefficient')
% figure()
% hold on 
% plot(drc,'b')
% plot(drn,'r')
% legend('CRT','NRT')
%figure()
%hold on 
%plot(Lam(1,:),'b')
%plot(Lam(2,:),'r')
%legend('CRT','NRT')
%Graphic_mdot(r,r2,Vx2,Vt2,H02)