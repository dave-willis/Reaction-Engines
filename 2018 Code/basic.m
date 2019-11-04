clear all 
%initial code used to get Smith Charts, possibly outdated if used for
%losses
%% Physical Parameters (given)
rho=7.6;
Wout=17*10^6; %W
Poin=145*10^5; %Pa
Toin=950; %K
mdot=16; %kg/s
Omega=6782*2*pi/60; %rad/s
cp=5190 ;
gam=1.667;
R=2080;
hoin=cp*Toin;

%% Geometric Parameters/Coefficients (set by me)

stages=12;
psi= [1.0:0.05:2.5];    %=((Wout/mdot)./stages)/((r.*Omega).^2);

phi=[0.2:0.01:0.5];

Lambda=0.5;                      %reaction

[PSI,PHI] = meshgrid(psi,phi);

%% Output Equations

r=(Wout./(mdot.*stages.*PSI.*Omega.^2)).^(0.5);
Vx = Omega.*r.*PHI;
Hin=mdot./(rho.*2*pi.*r.*Vx);

U=Omega.*r;
%% Flow Angles

alpha1=atand((1-0.5*PSI-Lambda)./PHI);
alpha3=alpha1;
alpha3rel=atand(tand(alpha1)-(1./PHI));    
alpha2=-alpha3rel;

alpha2rel= atand(tand(alpha2)-(1./PHI));
alpha1rel= alpha3rel;

%Every stage entry will have angle alpha1 (or alpha3), every stator exit flow angle is
%alpha2
%--------------------------------------------------------------------------
%% Changes in quantities
%flow velocities along a stage
V1 = Vx./cosd(alpha1);
V2 = Vx./cosd(alpha2);         %Stator exit velocity
V3 = Vx./cosd(alpha3);
%--------------------------------------------------------------------------
deltaho = -Wout/mdot; 
deltahostage = deltaho/stages; %unsure

%Inlet and Outlet Conditions
Tin=Toin-0.5*(V1.^2)./cp; %31x31
Toout = Toin + deltaho./cp;    %scalar
Toutisen = Toout - (V3.^2)./(2*cp); %31x31

Pin = Poin -0.5.*rho.*(V1.^2);  %31x31
Poutisen = Pin.*((Toutisen./Tin).^(gam./(gam-1)));   %31x31

%--------------------------------------------------------------------------
% T2 = To2-0.5.*Vx.^2./cp;
% M2 = V2./sqrt(gam*R.*T1);  %stator exit Mach no
%-------------------------------------------
%% Try Varying Annulus with const radius

Pstageisen = zeros(31, 31, stages);     %measured at stage exit
Tstageisen = zeros(31, 31, stages);   
hostage = zeros(31,31,stages);
hstage = zeros(31,31,stages);
Hstage = zeros(31,31,stages);
rhostageisen = zeros(31,31,stages);
M3stage = zeros(31,31,stages);

Tstageisen(:,:, 1)= Tin +((deltahostage-0.5.*(V3.^2))./cp);
Pstageisen(:,:, 1)= Pin.*((Tstageisen(:,:, 1)./Tin).^(gam./(gam-1)));
hostage(:,:,1)=hoin;
Hstage(:,:,1) = Hin;
rhostageisen(:,:,1) = Pin./(Tin*R);


for i = 1:stages-1
    
    hostage(:,:,i+1)= hostage(:,:,i)+deltahostage;
    hstage(:,:,i) = hostage(:,:,i)-0.5.*(V3.^2);
   
    Tstageisen(:,:, i+1) = Tstageisen(:,:,i) + ((deltahostage-0.5.*(V3.^2))./cp);
    Pstageisen(:,:,i+1) = Pstageisen(:,:,i).*((Tstageisen(:,:,i+1)./Tstageisen(:,:,i)).^(gam./(gam-1)));
    rhostageisen(:,:,i+1) = Pstageisen(:,:,i+1)./(R*Tstageisen(:,:,i+1));
    
    Hstage(:,:,i+1)=rhostageisen(:,:,i).*Hstage(:,:,i)./rhostageisen(:,:,i+1);
    
    M3stage(:,:,i+1) = V3./sqrt(gam*R*Tstageisen(:,:,i+1));
   
    
    
    
end
 
 
 %% Loss Calculations
%first round fix AR and cycle through stages
%second round iterate through AR & fix stage
%USING STAGE 5 AS FIXED POINT

AR = 0.2:0.2:5; 
% Cx = zeros(31,31,stages);
Cpb = -0.15;
tte = 0.001;
Cx = zeros(31,31,length(AR));
%set up chord lengths at each stage
for j =1:length(AR)
    Cx(:,:,j) = Hstage(:,:,5)./AR(j);
end

Cd = 0.002;
Z = 0.8;
p=zeros(31,31, length(AR)); %pitch
w=zeros(31,31, length(AR)); %throat area
Vbar = (V2+V1)./2;
deltaV = (V2-V1)./2;
%--------------------------------------------------------------------------
%Profile Loss Calc
zetaprofile=zeros(31,31, length(AR));
deltasprofile=zeros(31,31,length(AR));

for i = 1:length(AR)
    p(:,:,i) = (Z*Cx(:,:,i)./2)./(cosd(alpha2rel).^2.*(-tand(alpha1rel)+tand(alpha2rel)));
    w(:,:,i)=p(:,:,i).*cosd(alpha2);
    
    zetaprofile(:,:,i) = Cd.*(2.*sqrt(3) + 6.*1/sqrt(3)).*(tand(alpha2)-tand(alpha1));
    deltasprofile(:,:,i)=zetaprofile(:,:,i).*0.5.*Vbar.^2./Tstageisen(:,:,5);
    %^Profile Loss
end
%--------------------------------------------------------------------------
%Trailing Edge Loss Calc
deltasTE=zeros(31,31,length(AR));
zetaTE = -Cpb.*tte./w; %Trailing Edge Loss per stage

for j =1:length(AR)
    deltasTE(:,:,j)=zetaTE(:,:,j).*(hostage(:,:,5)-hstage(:,:,5))./(Tstageisen(:,:,5));
end
%--------------------------------------------------------------------------
%endwall scrubbing
Cdwall=0.0014;
%CAN STILL USE 3D BEC T FALLS LINEARLY W STAGES, similar ish for rho
endwall1=zeros(31,31,length(AR));
endwall2=zeros(31,31,length(AR));
endwallmid=zeros(31,31,length(AR));

for i=1:length(AR)
endwall1(:,:,i) = 0.25.*Cdwall.*(V1.^3).*rhostageisen(:,:,5).*p(:,:,i).*Cx(:,:,i)./Tstageisen(:,:,5);
endwall2(:,:,i) = 0.25.*Cdwall.*(V2.^3).*rhostageisen(:,:,5).*p(:,:,i).*Cx(:,:,i)./Tstageisen(:,:,5);
endwallmid(:,:,i) = Cdwall.*((V1.^3+V2.^3)./2).*rhostageisen(:,:,5).*p(:,:,i).*Cx(:,:,i)./Tstageisen(:,:,5);
end

deltasendwall = 2.*(endwall1+endwall2+endwallmid)./mdot;

%--------------------------------------------------------------------------
%secondary flow calc
alpham=0;
Ysec=zeros(31,31,length(AR));
deltas2nd=zeros(31,31,length(AR));
for i = 1:length(AR)
   Ysec(:,:,i)=(0.375.*0.1336).*Cx(:,:,i).*(cosd(alpha2).^3).*((tand(alpha1)-tand(alpha2)).^2)./...
       (Hin.*sqrt(cosd(alpha1)).*cosd(alpham));
   
   deltas2nd(:,:,i)=Ysec(:,:,i).*(hostage(:,:,5)-hstage(:,:,5))./(Tstageisen(:,:,5));
end
%--------------------------------------------------------------------------

stot = deltasendwall+deltasTE+deltasprofile+deltas2nd;

%--------------------------------------------------------------------------

%Efficiency Calc
%for stage 5
hloss=zeros(31,31,length(AR));
holoss=zeros(31,31,length(AR));
effstatic=zeros(31,31,length(AR));
efftotal=zeros(31,31,length(AR));

for i= 1:length(AR)
 holoss(:,:,i)=stot(:,:,i).*(Tstageisen(:,:,5)+0.5.*(Vbar.^2)./cp);
 hloss(:,:,i) = holoss(:,:,i);
 effstatic(:,:,i)=(-deltahostage-holoss(:,:,i))./(-deltahostage+0.5.*(V2.^2));
 efftotal(:,:,i)=(-deltahostage)./(-deltahostage+holoss(:,:,i));
 
end

%% Loss Plots
% figure(7)
%  surf(PSI, PHI, deltasTE(:,:,5), 'EdgeColor','r')
%  hold on
%  surf(PSI, PHI, deltasTE(:,:,14), 'EdgeColor','g')
%  hold on
%  surf(PSI, PHI, deltasTE(:,:,22), 'EdgeColor','b ')
%  hold off
%  title('TE \Delta S - AR = 1, 2.8, 4.4')
%  xlabel('\psi')
%  ylabel('\phi')
%  
%  figure(8)
%  surf(PSI, PHI, deltasprofile(:,:,5), 'DisplayName',num2str(1))
%  hold on
%  surf(PSI, PHI, deltasprofile(:,:,14), 'DisplayName',num2str(5))
%  hold on
%  surf(PSI, PHI, deltasprofile(:,:,22), 'DisplayName',num2str(9))
%  hold off
% title('\Delta S_{profile} (independent of AR) ')
%  xlabel('\psi')
%  ylabel('\phi') 
 
%  figure(13)
%  surf(PSI, PHI, deltasendwall(:,:,5), 'EdgeColor','r')
%  hold on
%  surf(PSI, PHI, deltasendwall(:,:,14), 'EdgeColor','g')
%  hold on
%  surf(PSI, PHI, deltasendwall(:,:,22), 'EdgeColor','b')
%  hold off
% title('\Delta S_{endwall}- - AR = 1, 2.8, 4.4' )
%  xlabel('\psi')
%  ylabel('\phi') 
%  
%  figure(6)
%  mesh(PSI, PHI, stot(:,:,5), 'EdgeColor','r')
%   hold on
%  surf(PSI, PHI, stot(:,:,14), 'EdgeColor','g')
%   hold on
%  surf(PSI, PHI, stot(:,:,22), 'EdgeColor','b')
%  title('Total Entropy Rise')
%  
% figure(15)
%  surf(PSI, PHI, effstatic(:,:,5), 'EdgeColor','r')
%  hold on
%  surf(PSI, PHI, effstatic(:,:,14), 'EdgeColor','g')
%  hold on
%  surf(PSI, PHI, effstatic(:,:,22), 'EdgeColor','b ')
%  hold off
%  title('Total to Static Efficiency AR = 1, 2.8, 4.4')
%  xlabel('\psi')
%  ylabel('\phi')
%  
% figure(16)
%  surf(PSI, PHI, efftotal(:,:,5), 'EdgeColor','r')
%  hold on
%  surf(PSI, PHI, efftotal(:,:,14), 'EdgeColor','g')
%  hold on
%  surf(PSI, PHI, efftotal(:,:,22), 'EdgeColor','b ')
%  hold off
%  title('Total to Total Efficiency AR = 1, 2.8, 4.4')
%  xlabel('\psi')
%  ylabel('\phi')
%  %% Test Value Plots (i.e at a particular PHI-PSI)
% 
% testp=squeeze(Pstageisen(28,25,:));
% testt=squeeze(Tstageisen(28,25,:));
% testH=squeeze(Hstage(28,25,:));
% testrho = squeeze(rhostageisen(28,25,:));
% testdeltasTE = squeeze(deltasTE(28,25, :));
% testdeltasprofile = squeeze(deltasprofile(28,25, :));
% testchord = squeeze(Cx(28,25, :));
% testpitch = squeeze(p(28,25, :));
% testendwall = squeeze(deltasendwall(28,25,:));
% testendwall1 = squeeze(endwall1(28,25,:));
% testendwall2 = squeeze(endwall2(28,25,:));
% testendwallmid = squeeze(endwallmid(28,25,:));
% teststot =  squeeze(stot(28,25,:));
% testdeltas2nd = squeeze(deltas2nd(28,25,:));
% testhloss = squeeze(holoss(28,25,:));
% % figure(9)
% % plot(AR,testpitch)
% % hold on
% % plot(AR,testchord)
% % hold off
% % xlabel('AR')
% % ylabel('pitch - Stage 5')
% 
% figure(10)
% plot(AR,testdeltasTE, 'DisplayName', 'TE')
% hold on
% plot(AR,testdeltasprofile, 'DisplayName', 'Profile')
% hold off
% xlabel('AR')
% ylabel('\Delta S_{TE} & \Delta S_{profile} at Stage 5')
% 
% figure(11)
% plot(AR,testendwall, 'DisplayName', '\Delta S_{endwall} ')
% xlabel('AR')
% ylabel('Endwall Loss at Stage 5')
% 
% figure(12)
% plot(AR, teststot)
% xlabel('AR')
% ylabel('Total Specific Entropy Rise')
% 
% figure(14)
% plot(AR, testdeltas2nd)
% xlabel('AR')
% ylabel('Secondary Flow Entropy Rise')
% 
% figure(17)
% plot(AR, testhloss)
% xlabel('AR')
% ylabel('Static Enthalpy Loss at stage 5')
%% operating point plots
figure(1)
contour(PHI, PSI, r,'ShowText', 'on')

hold on
contour(PHI,PSI,Hin,'ShowText','on')
hold off
xlabel('\phi')
ylabel('\psi')
title('radius & H ')
figure(2)
contour(PHI, PSI, Vx,'ShowText', 'on')
xlabel('\phi')
ylabel('\psi')
title('Vx')

figure(3)
contour(PHI, PSI, U,'ShowText', 'on')
xlabel('\phi')
ylabel('\psi')
title('Blade Speeds')


% figure(4)
% plot(1:stages, testH)
% xlabel('Stages')
% ylabel('H')
% title('Stage Exit H')

figure(5)
contour(PHI, PSI, alpha2,'ShowText', 'on')
xlabel('\phi')
ylabel('\psi')
title('Stator Exit angle')




