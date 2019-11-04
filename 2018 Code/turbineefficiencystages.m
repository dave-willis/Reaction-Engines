function out= turbineefficiencystages(phi, psi, stages, AR, S,g,G)
%finds losses for fixed phi and psi, loop in handoptimiser for multiple
%values for stages. Preferably input an array for AR rather than a single
%value
%S is a binary value for shrouded (1) or Unshrouded (0)
%g is the tip clearance
%G is the gear ratio
%% Physical Parameters (given)
rho=7.28; %assumed, at inlet
Wout=17*10^6; %W
Poin=145*10^5; %Pa
Toin=950; %K
mdot=16; %kg/s
Omega=G*6782*2*pi/60; %rad/s
cp=5190 ;
gam=1.667;
R=2080;
hoin=cp*Toin;
mu=31*10^-6; 
%% Geometric Parameters/Coefficients (set by me)
Lambda=0.5;                      %reaction
%% Output Equations

r=(Wout./(mdot.*stages.*psi.*Omega.^2)).^(0.5);
Vx = Omega.*r.*phi;
Hin=mdot./(rho.*2*pi.*r.*Vx);

U=Omega.*r; %dim(psi)
%% Flow Angles

alpha1=atand((1-0.5*psi-Lambda)./phi);
alpha3=alpha1;
alpha3rel=atand(tand(alpha1)-(1./phi));    
alpha2=-alpha3rel;

alpha2rel= atand(tand(alpha2)-(1./phi));
alpha1rel= alpha3rel;

%Every stage entry will have angle alpha1 (or alpha3), every stator exit flow angle is
%alpha2
%-------------------------------------------------------------------------
%% Changes in quantities
V1 = Vx./cosd(alpha1);
V2 = Vx./cosd(alpha2);         %Stator exit velocity
V3 = Vx./cosd(alpha3);
V2rel = Vx./cosd(alpha2rel);
V3rel = Vx./cosd(alpha3rel);
deltaho = -Wout/mdot;          %scalar
Toout = Toin + deltaho./cp;    %scalar
Toutisen = Toout - (V3.^2)./(2*cp); %dim(phi./psi)

deltahostage = linspace(deltaho./stages,deltaho,stages); %Cumulative, defined AT STAGE EXIT

Tin=Toin-0.5*(V1.^2)./cp; 
%--------------------------------------------------------------------------
%% Try Varying Annulus with const radius (with iteration through different max stages)
%rows = what stage we're on

Pin = Poin -0.5.*rho.*(V1.^2);  
%Poutisen = Pin.*((Toutisen./Tin).^(gam./(gam-1)));   

%measured at stage exit:
hostage=hoin + deltahostage;
hstage = hostage-0.5.*(V2.^2);   
Tstageisen = Tin +((deltahostage-0.5.*(V3.^2))./cp);
Pstageisen= Pin.*((Tstageisen./Tin).^(gam./(gam-1)));


%rhoin = Pin./(Tin*R);
rhostageisen = Pstageisen./(R*Tstageisen);
Hstage=zeros(1,stages);
Hstage(1)=rho.*Hin./rhostageisen(1);
for i=2:stages
Hstage(i)= rhostageisen(i-1).*Hstage(i-1)./rhostageisen(i);
end
M3relstage = (Vx./cosd(alpha3rel))./sqrt(gam*R.*Tstageisen);
M3stage = V3./sqrt(gam*R.*Tstageisen);    

      
%% Loss Calcs
Cpb = -0.15;
tte = 0.001; %TRAILING EDGE THICKNESS!
%set up chord lengths at each stage
Cx=zeros(stages, length(AR));

for i =1:length(AR)
Cx(:,i) = Hstage./AR(i);
end
Z = 0.8;
p=zeros(stages, length(AR)); %pitch
w=zeros(stages, length(AR)); %throat width
Vbar = (V2+V1)./2;
%--------------------------------------------------------------------------
%Profile Loss Calc
zetaprofile  =zeros(stages, length(AR));
deltasprofile=zeros(stages,length(AR));
Restages = zeros(stages,length(AR));
Cd = zeros(stages,length(AR));
for j = 1:stages
     Restages(j,:)=rhostageisen(j).*V2.*Cx(j,:)./mu;
    Cd(j,:)=0.002.*(Restages(j,:)./500000).^-0.2;
end

for i = 1:length(AR)
    p(:,i) = (Z*Cx(:,i)./2)./(cosd(alpha2).^2.*abs(tand(alpha1)-tand(alpha2)));
    w(:,i)=p(:,i).*cosd(alpha2);
    
    zetaprofile(:,i) = Cd(:,i).*(2.*sqrt(3) + 6.*1/sqrt(3)).*((tand(alpha2)-tand(alpha1))+abs(tand(alpha3rel)-tand(alpha2rel)));
    deltasprofile(1,i)=zetaprofile(1,i).*0.5.*Vbar.^2./Tstageisen(1);
for j = 2:stages
    deltasprofile(j,i)=zetaprofile(j,i).*0.5.*Vbar.^2./Tstageisen(j)+deltasprofile(j-1,i);
    %^Profile Loss
end
end
%--------------------------------------------------------------------------
%Trailing Edge Loss Calc
deltasTE=zeros(stages,length(AR));
zetaTE = 2*(-Cpb.*tte./w + (tte./w +1.4.*zetaprofile./2).^2); %Trailing Edge Loss per stage

for i =1:length(AR)
    deltasTE=zetaTE(1,:).*(hostage(1)-hstage(1))./(Tstageisen(1));
    for j =2:stages
        deltasTE(j,:) = deltasTE(j-1,:) + zetaTE(j,:).*(hostage(j)-hstage(j))./(Tstageisen(j));
    end
end
%--------------------------------------------------------------------------
%endwall scrubbing
Cdwall=0.002*ones(stages, length(AR)); %used a matrix of ones in case you want to replace 0.002 with a function

endwall1=zeros(stages,length(AR));
endwall3=zeros(stages,length(AR));
endwallmid=zeros(stages,length(AR));

for i=1:length(AR)
    endwall1(1,i) = 0.25.*Cdwall(1,i).*(V1.^3).*rhostageisen(1).*2.*pi.*(2.*r).*Cx(1,i)./Tstageisen(1);
    endwall3(1,i) = 0.25.*Cdwall(1,i).*(V3rel.^3).*rhostageisen(1).*2.*pi.*(2.*r).*Cx(1,i)./Tstageisen(1);
    endwallmid(1,i) = 0.5.*Cdwall(1,i).*((V2rel.^3+V2.^3)./2).*rhostageisen(1).*2.*pi.*(2.*r).*Cx(1,i)./Tstageisen(1);
    for j=2:stages
        endwall1(j,i) = endwall1(j-1,i)+0.25.*Cdwall(j,i).*(V1.^3).*rhostageisen(j).*2.*pi.*(2.*r).*Cx(j,i)./Tstageisen(j);
        endwall3(j,i) =endwall3(j-1,i)+ 0.25.*Cdwall(j,i).*(V3rel^3).*rhostageisen(j).*2.*pi.*(2.*r).*Cx(j,i)./Tstageisen(j);
        endwallmid(j,i) =endwallmid(j-1,i)+ 0.5.*Cdwall(j,i).*((V2rel.^3+V2.^3)./2).*rhostageisen(j).*2.*pi.*(2.*r).*Cx(j,i)./Tstageisen(j);
    end
end

deltasendwall = (endwall1+endwall3+endwallmid)./mdot;

%--------------------------------------------------------------------------
%secondary flow calc
alpham=-acotd(0.5.*(cotd(alpha1)+cotd(alpha2)));

Ysec=zeros(stages,length(AR));
deltas2nd=zeros(stages ,length(AR));

for i = 1:length(AR)
    
    Ysec(1,i)=2.*(0.375.*0.1336).*Cx(1,i).*(cosd(alpha2).^3).*((tand(alpha1)-tand(alpha2)).^2)./...
        (Hstage(1).*sqrt(cosd(alpha1)).*cosd(alpham));
    
    deltas2nd(1,i)=Ysec(1,i).*(hostage(1)-hstage(1))./(Tstageisen(1));
    for j=2:stages
        Ysec(j,i)=2.*(0.375.*0.1336).*Cx(j,i).*(cosd(alpha2).^3).*((tand(alpha1)-tand(alpha2)).^2)./...
            (Hstage(j).*sqrt(cosd(alpha1)).*cosd(alpham));
        
        deltas2nd(j,i)=Ysec(j,i).*(hostage(j)-hstage(j))./(Tstageisen(j)) + deltas2nd(j-1,i);
    end
end

%--------------------------------------------------------------------------
%Tip Clearance Loss

if S > 0
       Cc=0.6;
    'shrouded'
     
    fracleakagerot = zeros(stages, length(AR));
    fracleakagestat = zeros(stages, length(AR));
    for i = 1:length(AR)
    fracleakagerot(:,i)=g.*Cc.*sqrt(-((secd(alpha3).^2)-(tand(alpha2).^2)))./Hstage;
    fracleakagerot(:,i)=g.*Cc.*sqrt((secd(alpha3rel).^2)-(tand(alpha2rel).^2))./Hstage;
    fracleakagestat(:,i)= fracleakagerot(:,i);
    end
   
    deltastip(1,:)=(((Vx./cosd(alpha3rel)).^2).*fracleakagerot(1,:).*(1-tand(alpha2rel).*(sind(alpha3rel).^2)./tand(alpha3rel))...
        + (V2.^2).*fracleakagestat(1,:).*(1-tand(alpha1).*(sind(alpha2).^2)./tand(alpha2)))./Tstageisen(1);
    for j=2:stages 
    deltastip(j,:)=((Vx./cosd(alpha3rel).^2).*fracleakagerot(j,:).*(1-tand(alpha2rel).*(sind(alpha3rel).^2)./tand(alpha3rel))...
        + (V2.^2).*fracleakagestat(j,:).*(1-tand(alpha1).*(sind(alpha2).^2)./tand(alpha2)))./Tstageisen(j)...
        +deltastip(j-1,:);
    end
else
     'unshrouded'
     Cdtip=0.8;
    deltastip=zeros(stages, length(AR));
    deltastip(1,:)=(Cdtip.*g.*Cx(1,:).*0.65).*((Vx./cosd(alpha3rel)).^2)./(Hstage(1).*p(1,:).*cosd(alpha3rel).*Tstageisen(1))...
        +(Cdtip.*g.*Cx(1,:).*0.65).*(V2.^2)./(Hstage(1).*p(1,:).*cosd(alpha2).*Tstageisen(1));
    for j=2:stages
    deltastip(j,:)=deltastip(j-1,:)+Cdtip.*g.*Cx(j,:).*0.65.*(((Vx./cosd(alpha3rel)).^2)./Hstage(j).*p(j,:).*cosd(alpha3rel).*Tstageisen(j))...
        +Cdtip.*g.*Cx(j,:).*0.65.*(V2.^2)./(Hstage(j).*p(j,:).*cosd(alpha2).*Tstageisen(j));
        
    end

end
%--------------------------------------------------------------------------

sovr = deltasendwall+deltasTE+deltasprofile+deltas2nd+deltastip;
stot = sovr(end,:);
%--------------------------------------------------------------------------

%Efficiency Calc

 holosstot = stot.*(Toout);
 
%  effstatic =(-deltaho)./(-deltaho+0.5.*(V3.^2)+holosstot);
 efftotal =(-deltaho)./(-deltaho+holosstot);

%Lost Efficiencies
losteffTE = Toout.*deltasTE(end,:)./(-deltaho);
losteffprofile = Toout.*deltasprofile(end,:)./(-deltaho );
losteffendwall = Toout.*deltasendwall(end,:)./(-deltaho);
losteff2nd = Toout.*deltas2nd(end,:)./(-deltaho);
lostefftip = Toout.*deltastip(end,:)./(-deltaho);
lostefftotal = losteffTE+losteffprofile+losteffendwall+losteff2nd+lostefftip;

%% Outputs

stagewidth=3.*Cx;
volstage=zeros(stages,length(AR));
for j=1:stages
volstage(j,:) = pi*stagewidth(j,:).*(r+Hstage(j)./2).^2;
end
vol=sum(volstage);

if length(AR)>1
figure(1)
scatter(stages,losteffTE(7),'b')
hold on
scatter(stages,losteffprofile(7),'r')
hold on
scatter(stages,losteffendwall(7),'g')
hold on
scatter(stages,losteff2nd(7),'c')
hold on
scatter(stages,lostefftip(7),'m')
legend('TE','Profile','Endwall','Secondary','Tip Clearance')
hold on
% hold off
xlabel('No. of Stages')
ylabel('T_0\Delta s/\Delta h_0 ')
title('Loss Breakdown(AR = 1.4)')

figure(2)
plot(AR, efftotal, 'DisplayName', num2str(stages))
hold on
% plot(AR, effstatic,'DisplayName','Total-to-Static')
% hold off
xlabel('Aspect Ratio')
ylabel('Efficiency')
title('\eta_{tt} variation across AR, for different n_{stages}')

figure(3)
plot(AR, vol/(pi.*r.^3),'DisplayName', num2str(stages))
hold on
xlabel('AR')
ylabel('L/r_{mean}')
title('Turbine Length to Mean Radius')


figure(6)
plot(AR, vol,'DisplayName', num2str(stages))
xlabel('Aspect Ratio')
ylabel('Volume (m^3)')
title('Total Volume vs Aspect Ratio')
hold on
else
    
    
    efftotal
    
    vol
end 
figure(4)
plot(1:stages, Hstage, 'DisplayName',num2str(stages))
hold on
% line([0,1],[r-Hin,r-Hstage(1)])
% hold on
% line([0,stages],[r,r])
grid on
grid minor
xlabel('Stages through Turbine')
ylabel('Span (m)')
%ylim([0.4 1])
%axis([0 inf 0 1])
title('Span Variation')
% figure(5)
% plot(1:stages, M3relstage, 'DisplayName', num2str(stages))
% hold on
% % plot(1:stages, M3stage, 'DisplayName', 'Absolute')
% % hold off
% xlabel('Stages')
% ylabel('M')
% title('Stage Exit Relative Mach Numbers')

figure(7)
scatter(stages,max(efftotal),'DisplayName', num2str(stages))
xlabel('n_{stages}')
ylabel('\eta')
title('Maximum Efficiency at Different n_{stages} (Optimal AR)')
hold on


figure(8)
plot(AR, losteffendwall, 'DisplayName', num2str(stages))
ylim([0 inf])
xlabel('AR')
ylabel('T_0\Delta s_{endwall}/\Delta h_0')
title('Endwall Scrubbing Loss')
hold on

figure(18)
plot(AR, losteff2nd, 'DisplayName', num2str(stages))
ylim([0 inf])
xlabel('AR')
ylabel('T_0\Delta s_{Sec}/\Delta h_0')
title('Secondary Flow Loss')
hold on

figure(17)
plot(AR, lostefftip, 'DisplayName', num2str(stages))
ylim([0 inf])
xlabel('AR')
ylabel('T_0\Delta s_{tip}/\Delta h_0')
title('Tip Leakage Loss')
hold on

figure(16)
plot(AR, losteffTE, 'DisplayName', num2str(stages))
ylim([0 inf])
xlabel('AR')
ylabel('T_0\Delta s_{TE}/\Delta h_0')
title('Trailing Edge Loss')
hold on

figure(15)
plot(AR, losteffprofile, 'DisplayName', num2str(stages))
ylim([0 inf])
xlabel('AR')
ylabel('T_0\Delta s_{profile}/\Delta h_0')
title('Profile Loss')
hold on

out(1:stages, 1)=r;
out(1:stages, 2)=alpha1;
out(1:stages, 3)=alpha2;
out(1:stages, 4)=Vx;
out(1:stages, 5)=Tstageisen';
out(1:stages, 6)=Hstage';
end


