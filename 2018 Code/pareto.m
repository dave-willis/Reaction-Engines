%% Description
% Runs a Pareto Optimum Analysis for designs with constant mean
% radius. Plots Pareto Curve, stores ordered Pareto Set and does an Axial
% Force and Vessel thickness analysis
clear 
clc
%hold off
%% Physical Parameters (given)
rho=7.28; %assumed at inlet
Wout=17*10^6; %W
Poin=145*10^5; %Pa
Toin=950; %K
mdot=16; %kg/s
G=1; %Gear Ratio
Omega=G*6782*2*pi/60; %rad/s
cp=5190 ;
gam=1.667;
R=2080;
hoin=cp*Toin;
mu=31*10^-6; 
%% Geometric Parameters/Coefficients (set by me)
for z=40 %set as range of values to plot multiple curves
    for g =  0.0005;%0.0002:0.0002:0.001

stages=z;
psi= [0.8:0.05:2.4];    %psi and phi must be same size

phi=[0.1:0.025:0.9];

n=length(phi);
Lambda=0.5;                      %reaction

[phi,psi] = meshgrid(phi,psi);

AR=0.6:0.1:3;
%% Output Equations

deltaho = -Wout/mdot;          %scalar
deltahostage = linspace(0,deltaho,stages+1); %AT STAGE EXIT

r=(Wout./(mdot.*stages.*psi.*Omega.^2)).^(0.5);
Vx = Omega.*r.*phi;
Hin=mdot./(rho.*2*pi.*r.*Vx);

U=Omega.*r; %dim(psi)
%% Flow Angles

alpha1=atand((1-(0.5.*psi)-Lambda)./phi);
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

Toout = Toin + deltaho./cp;    %scalar
Toutisen = Toout - (V3.^2)./(2*cp); %dim(phi./psi)


Tin=Toin-0.5*(V1.^2)./cp; 
%--------------------------------------------------------------------------
%% Stage By Stage Quantities


Pin = Poin -0.5.*rho.*(V1.^2);  
%Poutisen = Pin.*((Toutisen./Tin).^(gam./(gam-1)));   

%measured at stage exit:
hostage=hoin + deltahostage(2:end);
hstage =zeros(n,n,stages); 
Tstageisen=zeros(n,n,stages); 
Pstageisen=zeros(n,n,stages); 
rhostageisen=zeros(n,n,stages); 
Hstage=zeros(n,n,stages); 


M3stage=zeros(n,n,stages); 
M3relstage=zeros(n,n,stages); 
Tostageisen = zeros(n,n,stages);
for j=1:stages
hstage(:,:,j) = hostage(j)-0.5.*(V2.^2);   
Tstageisen(:,:,j) = Tin +((deltahostage(j+1)-0.5.*(V3.^2))./cp);
Pstageisen(:,:,j)= Pin.*((Tstageisen(:,:,j)./Tin).^(gam./(gam-1)));

rhostageisen = Pstageisen./(R*Tstageisen);
Hstage(:,:,1)=rho.*Hin./rhostageisen(:,:,1);  
Tostageisen(:,:,j)=Toin + deltahostage(j+1)./cp;%Tstageisen(:,:,j)+V3.^2./cp;
if j>1
Hstage(:,:,j)= rhostageisen(:,:,j-1).*Hstage(:,:,j-1)./rhostageisen(:,:,j);
end
M3relstage = (Vx./cosd(alpha3rel))./sqrt(gam*R.*Tstageisen(:,:,j));
M3stage = V3./sqrt(gam*R.*Tstageisen(:,:,j));    
end
      
%% Loss Calcs
m=length(AR);
Cpb = -0.15;
tte = 0.0005; %trailing edge thickness!
%set up chord lengths at each stage
Cx=zeros(n,n,stages, m);

for i =1:m
Cx(:,:,:,i)  = Hstage./AR(i);
end
Z = 0.8;
p=zeros(n,n,stages, m); %pitch
w=zeros(n,n,stages, m); %throat width
Vbar = (V2+V1)./2;
%--------------------------------------------------------------------------
%Profile Loss Calc
zetaprofile=zeros(n,n,stages, m);
deltasprofile=zeros(n,n,stages, m);
Restages = zeros(n,n,stages, m);
Cd = zeros(n,n,stages, m);
for i = 1:m
    for j = 1:stages
        Restages(:,:,j,i)=rhostageisen(:,:,j).*V2.*Cx(:,:,j,i)./mu;
        Cd(:,:,j,i)=0.002.*(Restages(:,:,j,i)./500000).^-0.2;
    end
end

for j = 1:stages
    for i = 1:m
        p(:,:,j,i) = (Z*Cx(:,:,j,i)./2)./(cosd(alpha2).^2.*abs(tand(alpha1)-tand(alpha2)));
        w(:,:,j,i)=p(:,:,j,i).*cosd(alpha2);
        
        zetaprofile(:,:,j,i)= Cd(:,:,j,i).*(2.*sqrt(3) + 6.*1/sqrt(3)).*((tand(alpha2)-tand(alpha1))+abs(tand(alpha3rel)-tand(alpha2rel)));
        deltasprofile(:,:,1,i)=zetaprofile(:,:,1,i).*0.5.*Vbar.^2./Tstageisen(:,:,1);
        if j>1
            deltasprofile(:,:,j,i)=zetaprofile(:,:,j,i).*0.5.*Vbar.^2./Tstageisen(:,:,j)+deltasprofile(:,:,j-1,i);
            %^Profile Loss
        end
    end
end
%--------------------------------------------------------------------------
%Trailing Edge Loss Calc
deltasTE=zeros(n,n,stages,m);
zetaTE = 2*(-Cpb.*tte./w + (tte./w +1.4.*zetaprofile./2).^2); %Trailing Edge Loss per stage

for i =1:m
    deltasTE(:,:,1,i)=zetaTE(:,:,1,i).*(hostage(1)-hstage(:,:,1))./(Tstageisen(:,:,1));
    for j =2:stages
        deltasTE(:,:,j,i) = deltasTE(:,:,j-1,i) + zetaTE(:,:,j,i).*(hostage(j)-hstage(:,:,j))./(Tstageisen(:,:,j));
    end
end
%--------------------------------------------------------------------------
%endwall scrubbing
Cdwall=0.002*ones(n,n,stages, length(AR));
 
endwall1=zeros(n,n,stages,length(AR));
endwall3=zeros(n,n,stages,length(AR));
endwallmid=zeros(n,n,stages,length(AR));

for i=1:length(AR)
    endwall1(:,:,1,i) = 0.25.*Cdwall(:,:,1,i).*(V1.^3).*rhostageisen(:,:,1).*2.*pi.*(2.*r).*Cx(1,i)./Tstageisen(:,:,1);
    endwall3(:,:,1,i) = 0.25.*Cdwall(:,:,1,i).*(V3rel.^3).*rhostageisen(:,:,1).*2.*pi.*(2.*r).*Cx(:,:,1,i)./Tstageisen(:,:,1);
    endwallmid(:,:,1,i) = 0.5*Cdwall(:,:,1,i).*((V2rel.^3+V2.^3)./2).*rhostageisen(:,:,1).*2.*pi.*(2.*r).*Cx(:,:,1,i)./Tstageisen(:,:,1);
    for j=2:stages
        endwall1(:,:,j,i) = endwall1(:,:,j-1,i)+0.25.*Cdwall(:,:,j,i).*(V1.^3).*rhostageisen(:,:,j).*2.*pi.*(2.*r).*Cx(:,:,j,i)./Tstageisen(:,:,j);
        endwall3(:,:,j,i) =endwall3(:,:,j-1,i)+ 0.25.*Cdwall(:,:,j,i).*(V3rel.^3).*rhostageisen(:,:,j).*2.*pi.*(2.*r).*Cx(:,:,j,i)./Tstageisen(:,:,j);
        endwallmid(:,:,j,i) =endwallmid(:,:,j-1,i)+ 0.5*Cdwall(:,:,j,i).*((V2rel.^3+V2.^3)./2).*rhostageisen(:,:,j).*2.*pi.*(2.*r).*Cx(:,:,j,i)./Tstageisen(:,:,j);
    end
end

deltasendwall = (endwall1+endwall3+endwallmid)./mdot;

%--------------------------------------------------------------------------
%secondary flow calc
alpham=-acotd(0.5.*(cotd(alpha1)+cotd(alpha2)));%nxn

Ysec=zeros(n,n,stages,length(AR));
deltas2nd=zeros(n,n,stages ,length(AR));

for i = 1:length(AR)
    
    Ysec(:,:,1,i)=2*(0.375.*0.1336).*Cx(:,:,1,i).*(cosd(alpha2).^3).*((tand(alpha1)-tand(alpha2)).^2)./...
        (Hstage(:,:,1).*sqrt(cosd(alpha1)).*cosd(alpham));
    
    deltas2nd(:,:,1,i)=Ysec(:,:,1,i).*(hostage(1)-hstage(:,:,1))./(Tstageisen(:,:,1));
    for j=2:stages
        Ysec(:,:,j,i)=2*(0.375.*0.1336).*Cx(:,:,j,i).*(cosd(alpha2).^3).*((tand(alpha1)-tand(alpha2)).^2)./...
            (Hstage(:,:,j).*sqrt(cosd(alpha1)).*cosd(alpham));
        
        deltas2nd(:,:,j,i)=Ysec(:,:,j,i).*(hostage(j)-hstage(:,:,j))./(Tstageisen(:,:,j)) + deltas2nd(:,:,j-1,i);
    end
end

%--------------------------------------------------------------------------
%Tip Clearance Loss
S=1;

g;
if S > 0
       Cc=0.6;
    'shrouded'
     
    fracleakagerot = zeros(n,n,stages, length(AR));
    fracleakagestat = zeros(n,n,stages, length(AR));
    deltastip = zeros(n,n,stages, length(AR));
    for i = 1:length(AR)
        for j =1:stages
           
            fracleakagerot(:,:,j,i)=g.*Cc.*sqrt((secd(alpha3rel).^2)-(tand(alpha2rel).^2))./Hstage(:,:,j);
            fracleakagestat(:,:,j,i)= fracleakagerot(:,:,j,i);
            
            deltastip(:,:,1,i)=(((Vx./cosd(alpha3rel)).^2).*fracleakagerot(:,:,1,i).*(1-tand(alpha2rel).*(sind(alpha3rel).^2)./tand(alpha3rel))...
                + (V2.^2).*fracleakagestat(:,:,1,i).*(1-tand(alpha1).*(sind(alpha2).^2)./tand(alpha2)))./Tstageisen(:,:,1);
        if j>1
            deltastip(:,:,j,i)=((Vx./cosd(alpha3rel).^2).*fracleakagerot(:,:,j,i).*(1-tand(alpha2rel).*(sind(alpha3rel).^2)./tand(alpha3rel))...
                + (V2.^2).*fracleakagestat(:,:,j,i).*(1-tand(alpha1).*(sind(alpha2).^2)./tand(alpha2)))./Tstageisen(:,:,j)...
                +deltastip(:,:,j-1,i);
            
        end
        
        end
        
    
    end
   
    
    


else
     'unshrouded'
     Cdtip=0.8;
    deltastip=zeros(n,n,stages, length(AR));
    for i=1:m
        deltastip(:,:,1,i)=(Cdtip.*g.*Cx(:,:,1,i).*0.65).*((Vx./cosd(alpha3rel)).^2)./(Hstage(:,:,1).*p(:,:,1,i).*cosd(alpha3rel).*Tstageisen(:,:,1))...
            +(Cdtip.*g.*Cx(:,:,1,i).*0.65).*(V2.^2)./(Hstage(:,:,1).*p(:,:,1,i).*cosd(alpha2).*Tstageisen(:,:,1));
        for j=2:stages
            deltastip(:,:,j,i)=deltastip(:,:,j-1,i)+Cdtip.*g.*Cx(:,:,j,i).*0.65.*(((Vx./cosd(alpha3rel)).^2)./Hstage(:,:,j).*p(:,:,j,i).*cosd(alpha3rel).*Tstageisen(:,:,j))...
                +Cdtip.*g.*Cx(:,:,j,i).*0.65.*(V2.^2)./(Hstage(:,:,j).*p(:,:,j,i).*cosd(alpha2).*Tstageisen(:,:,j));
            
        end
    end
end
%--------------------------------------------------------------------------

sovr = deltasendwall+deltasTE+deltasprofile+deltas2nd+deltastip;
stot = squeeze(sovr(:,:,end,:));
deltas2nd = squeeze(deltas2nd);
deltasTE = squeeze(deltasTE);
deltasendwall = squeeze(deltasendwall);
deltasprofile = squeeze(deltasprofile);
deltastip = squeeze(deltastip);
%--------------------------------------------------------------------------

%Efficiency Calc

 holosstot = stot.*(Toout);
 
%  effstatic =(-deltaho)./(-deltaho+0.5.*(V3.^2)+holosstot);

efftotal =(-deltaho)./(-deltaho+holosstot);

%Lost Efficiencies
losteffTE = Toout.*deltasTE(:,:,end,:)./(-deltaho);
losteffprofile = Toout.*deltasprofile(:,:,end,:)./(-deltaho );
losteffendwall = Toout.*deltasendwall(:,:,end,:)./(-deltaho);
losteff2nd = Toout.*deltas2nd(:,:,end,:)./(-deltaho);
lostefftip = Toout.*deltastip(:,:,end,:)./(-deltaho);
lostefftotal = losteffTE+losteffprofile+losteffendwall+losteff2nd+lostefftip;

%% Plotting Pareto
stagewidth=3.*Cx;
volstage= zeros(n,n,stages,m);
for j=1:stages
for i=1:m
    volstage(:,:,j,i) = pi*stagewidth(:,:,j,i).*((r+Hstage(:,:,j)./2).^2);   
end
end
machinelength=sum(stagewidth,3);
vol = sum(volstage,3);
paretovol = vol;
paretoneff = 1-efftotal;

paretovol = reshape(paretovol, 1,[]);
paretoneff = reshape(paretoneff,1,[]);
%reshape will order as: vary psi, vary ph  i, vary AR
store=zeros(length(paretoneff),5);


for i=1:n*m
    store(n*(i-1)+1:n*(i),2)=psi(:,1);
   
end
phi2=reshape(phi,1,[]);

for j=1:m
    store((n*n*(j-1))+1:n*n*j,1)=phi2;
    store(n*n*(j-1)+1:n*n*(j),3)=repelem(AR(j),n*n);
end

store(:,4)= paretoneff';
store(:,5)= paretovol';
for i = 1:length(paretoneff)%point we're comparing to
    for j = 1:length(paretoneff) %point we're looking at
         
        if  store(j,4)~=0 && paretovol(j)>paretovol(i) && paretoneff(j)>paretoneff(i)
            
            store(j,:)=0;
%             paretovol(j)=0;
%             paretoneff(j)=0;
            
        end
    end
end

x=length(nonzeros(store(:,1)));
store2=zeros(x,5);
for col=1:5
    store2(:,col)=nonzeros(store(:,col));
end

mkr=['*','+','o','.','s'];
% figure(1)
%  scatter(paretovol, paretoneff, 5,'DisplayName', num2str(z))
%  hold on
%  xlabel('Volume')
%  ylabel('1 - \eta_{tt}')

figure(6)
scatter(store2(:,5),store2(:,4),3, mkr(5),'DisplayName', num2str(stages))
hold on
 xlabel('Volume')
 ylabel('1 - \eta_{tt}')
 
%% Shape Plots

Hplot = zeros(stages,3);
rtipplot= zeros(stages,3);
rhubplot = zeros(stages,3);
stagewidthplot = zeros(stages,3);
ttetowplot = zeros(stages,3);
hubtotipplot = zeros(stages,3);
nbladeplot = zeros(stages,3);

[~,idx] = sort(store2(:,4)); % sort just the 4th column
store3 = store2(idx,:);  


medianval= round(length(store3(:,1))./2);
rplot(1) = (Wout./(mdot.*stages.*store3(1,2).*Omega.^2)).^(0.5);
rplot(2) = (Wout./(mdot.*stages.*store3(medianval,2).*Omega.^2)).^(0.5);
rplot(3) = (Wout./(mdot.*stages.*store3(end,2).*Omega.^2)).^(0.5);
 
idxphi(1) = find(phi(1,:)==store3(1,1));
idxpsi(1) = find(psi(:,1)==store3(1,2));
Hinplot(1)=Hin(idxpsi(1),idxphi(1));
Hplot(:,1)=Hstage(idxpsi(1),idxphi(1),:);
rtipplot(:,1) = rplot(1)+Hplot(:,1)./2;
alpha2plot(1)=alpha2(idxpsi(1),idxphi(1));
alpha1plot(1)=alpha1(idxpsi(1),idxphi(1));
ptoCplot(1) = (Z./2)./(cosd(alpha2plot(1)).^2.*abs(tand(alpha1plot(1))-tand(alpha2plot(1))));
idxAR(1) = find(AR == store3(1,3));
stagewidthplot(:,1) = stagewidth(idxpsi(1),idxphi(1),:,idxAR(1));
lplot(1)=sum(stagewidthplot(:,1));
rhubplot(:,1)=  rplot(1)-Hplot(:,1)./2;
hubtotipplot(:,1)=rhubplot(:,1)./rtipplot(:,1);
ttetowplot(:,1)=tte./w(idxpsi(1),idxphi(1),:,idxAR(1));


idxphi(2) = find(phi(1,:)==store3(medianval,1));
idxpsi(2) = find(psi(:,1)==store3(medianval,2));
Hinplot(2)=Hin(idxpsi(2),idxphi(2));
Hplot(:,2)=Hstage(idxpsi(2),idxphi(2),:);
rtipplot(:,2) = rplot(2)+Hplot(:,2)./2;
alpha2plot(2)=alpha2(idxpsi(2),idxphi(2));
alpha1plot(2)=alpha1(idxpsi(2),idxphi(2));
ptoCplot(2) = (Z./2)./(cosd(alpha2plot(2)).^2.*abs(tand(alpha1plot(2))-tand(alpha2plot(2))));
idxAR(2) = find(AR == store3(medianval,3));
stagewidthplot(:,2) = stagewidth(idxpsi(2),idxphi(2),:,idxAR(2));
lplot(2)=sum(stagewidthplot(:,2));
rhubplot(:,2)=  rplot(2)-Hplot(:,2)./2;
hubtotipplot(:,2)=rhubplot(:,2)./rtipplot(:,2);
ttetowplot(:,2)=tte./w(idxpsi(2),idxphi(2),:,idxAR(2));

idxphi(3) = find(phi(1,:)==store3(end,1));
idxpsi(3) = find(psi(:,1)==store3(end,2));
Hinplot(3)=Hin(idxpsi(3),idxphi(3));
Hplot(:,3)=Hstage(idxpsi(3),idxphi(3),:);
rtipplot(:,3) = rplot(3)+Hplot(:,3)./2;
alpha2plot(3)=alpha2(idxpsi(3),idxphi(3));
alpha1plot(3)=alpha1(idxpsi(3),idxphi(3));
ptoCplot(3) = (Z./2)./(cosd(alpha2plot(3)).^2.*abs(tand(alpha1plot(3))-tand(alpha2plot(3))));
idxAR(3) = find(AR == store3(end,3));
stagewidthplot(:,3) = stagewidth(idxpsi(3),idxphi(3),:,idxAR(3));
lplot(3)=sum(stagewidthplot(:,3));
rhubplot(:,3)=  rplot(3)-Hplot(:,3)./2;
hubtotipplot(:,3)=rhubplot(:,3)./rtipplot(:,3);
ttetowplot(:,3)=tte./w(idxpsi(3),idxphi(3),:,idxAR(3));

pitchplot=zeros(stages,3);
for x=1:3
    pitchplot(:,x)=p(idxpsi(x),idxphi(x),:,idxAR(x));
    nbladeplot(:,x)=2.*pi.*rplot(x)./p(idxpsi(x),idxphi(x),:,idxAR(x));
end
figure(10)
plot((1:stages).*stagewidthplot(:,1)', rtipplot(:,1), (1:stages).*stagewidthplot(:,1)', rtipplot(:,1)-Hplot(:,1),[0 stagewidthplot(1,1)],[rplot(1)-Hinplot(1)./2 rtipplot(1,1)-Hplot(1,1)],[0 stagewidthplot(1,1)],[rplot(1)+Hinplot(1)./2 rtipplot(1,1)])
hold on
plot((1:stages).*stagewidthplot(:,3)', rtipplot(:,3), (1:stages).*stagewidthplot(:,3)', rtipplot(:,3)-Hplot(:,3),[0 stagewidthplot(1,3)],[rplot(3)-Hinplot(3)./2 rtipplot(1,3)-Hplot(1,3)],[0 stagewidthplot(1,3)],[rplot(3)+Hinplot(3)./2 rtipplot(1,3)])
 

ylim([0 0.7])

figure(11)
plot(1:stages, ttetowplot(:,1), 'DisplayName', 'Constant r_{mean}')
title('t_{te}/w')
hold on

figure(12)
plot(1:stages, nbladeplot(:,2), 'DisplayName', num2str(z))
title('n_{blades}')
hold on

figure(13)
scatter(z,abs(alpha1plot(2)-alpha2plot(2)))
title('flow turning')
hold on

figure(14)
plot(1:stages, pitchplot(:,2), 'DisplayName',num2str(z))
title('Pitch')
hold on

%% Axial Loads & Pressure Vessel
%postive = force on fluid in upstream direction

Postageisen=zeros(n,n,stages);
Pstagerotor=zeros(n,n,stages);
Postageisen(:,:,1)=Poin.*((Tostageisen(:,:,1)./Toin).^(gam./(gam-1)));
Pstagerotor(:,:,1)=Poin-0.5.*rho.*V2.^2; %evaluated at stage entry
%i.e, other 'stageisen' values are evaluated at stage exit, but we need to
%use values at stage entry to get sensible static pressures on either side
%of the rotor

Fx = zeros(n,n,stages);
Fx(:,:,1)= 2.*pi.*r.*(Pstageisen(:,:,1).*Hstage(:,:,1)-(Pstagerotor(:,:,1).*Hin+(Pstagerotor(:,:,1)+Pstageisen(:,:,1))./2.*(Hstage(:,:,1)-Hin)));
for j = 2:stages
    Postageisen(:,:,j) = Postageisen(:,:,j-1).*((Tostageisen(:,:,j)./Tostageisen(:,:,j-1)).^(gam./(gam-1)));
    Pstagerotor(:,:,j) = Postageisen(:,:,j-1)-0.5.*rhostageisen(:,:,j-1).*V2.^2;
    Fx(:,:,j) =  2.*pi.*r.*(Pstageisen(:,:,j).*Hstage(:,:,j)-(Pstagerotor(:,:,j).*Hstage(:,:,j-1)+(Pstagerotor(:,:,j)+Pstageisen(:,:,j))./2.*(Hstage(:,:,j)-Hstage(:,:,j-1))));
end
Fxplotstages=zeros(3,stages);
Fxplottot = zeros(1,3);
for x=1:3
Fxplotstages(x,:)= Fx(idxpsi(x),idxphi(x),:);
Fxplottot(x) =squeeze(sum(Fxplotstages(x,:),2));
figure(15)
scatter(x,Fxplottot(x))
hold on
end

sigmayield = 370*10^6;
vesselthickness=zeros(1,3);

for x = 1:3  
vesselthickness(x)=Pin(idxpsi(x),idxphi(x)).*max(rtipplot(:,x))./sigmayield;
end

    end
end
hold off
hold off