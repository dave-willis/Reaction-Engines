clear 
clc
%Runs a Pareto Optimum Analysis for designs with constant mean radius and
%chord as a function of stage (as opposed to constant AR)
%I couldn't get anything useful out of it, maybe you can!
%% Physical Parameters (given)
rho=7.28;
Wout=17*10^6; %W
Poin=145*10^5; %Pa
Toin=950; %K
mdot=16; %kg/s
G=1;
Omega=G*6782*2*pi/60; %rad/s
cp=5190 ;
gam=1.667;
R=2080;
hoin=cp*Toin;
mu=31*10^-6; 
%% Geometric Parameters/Coefficients (set by me)

for z=6:2:18
 for g = 0.0005
stages=z;
psi= [1.0:0.05:2.5];    

phi=[0.2:0.01:0.5];

n=length(phi);
Lambda=0.5;                      %reaction

[phi,psi] = meshgrid(phi,psi);
Cx=zeros(stages,11);
for j = 1:stages
Cx(j,:)=0.0008 +0.0063.*j.*[1:0.05:1.5] -0.0001.*(j.^2).*[1:0.05:1.5];
end
figure(10)
surf(Cx)
xlabel('Stages')
ylabel('Chord Factor Sweep')
zlabel('Chord Length')
hold on

%% Output Equations

%% Output Equations

r=(Wout./(mdot.*stages.*psi.*Omega.^2)).^(0.5);
Vx = Omega.*r.*phi;
Hin=mdot./(rho.*2*pi.*r.*Vx);

U=Omega.*r; %dim(psi)
%% Flow Angles

alpha1=atand((1-0.5.*psi-Lambda)./phi);
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

deltahostage = linspace(deltaho./stages,deltaho,stages); %AT STAGE EXIT

Tin=Toin-0.5*(V1.^2)./cp; 
%--------------------------------------------------------------------------
%% Stage By Stage Quantities


Pin = Poin -0.5.*rho.*(V1.^2);  
%Poutisen = Pin.*((Toutisen./Tin).^(gam./(gam-1)));   

%measured at stage exit:
hostage=hoin + deltahostage;
hstage =zeros(n,n,stages); 
Tstageisen=zeros(n,n,stages); 
Pstageisen=zeros(n,n,stages); 
rhostageisen=zeros(n,n,stages); 
Hstage=zeros(n,n,stages); 


M3stage=zeros(n,n,stages); 
M3relstage=zeros(n,n,stages); 
for j=1:stages
hstage(:,:,j) = hostage(j)-0.5.*(V2.^2);   
Tstageisen(:,:,j) = Tin +((deltahostage(j)-0.5.*(V3.^2))./cp);
Pstageisen(:,:,j)= Pin.*((Tstageisen(:,:,j)./Tin).^(gam./(gam-1)));

rhostageisen = Pstageisen./(R*Tstageisen);
Hstage(:,:,1)=rho.*Hin./rhostageisen(1);  

if j>1
Hstage(:,:,j)= rhostageisen(:,:,j-1).*Hstage(:,:,j-1)./rhostageisen(:,:,j);
end
M3relstage = (Vx./cosd(alpha3rel))./sqrt(gam*R.*Tstageisen(:,:,j));
M3stage = V3./sqrt(gam*R.*Tstageisen(:,:,j));    
end
      
%% Loss Calcs
m=length(Cx(1,:));
Cpb = -0.15;
tte = 0.001;
%set up chord lengths at each stage
AR=zeros(n,n,stages, m);

for j = 1:stages
for i =1:m
AR(:,:,j,i)  = Hstage(:,:,j)./Cx(j,i);
end
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
        Restages(:,:,j,i)=rhostageisen(:,:,j).*V2.*Cx(j,i)./mu;
        Cd(:,:,j,i)=0.002.*(Restages(:,:,j,i)./500000).^-0.2;
    end
end

for j = 1:stages
    for i = 1:m
        p(:,:,j,i) = (Z*Cx(j,i)./2)./(cosd(alpha2).^2.*abs(tand(alpha1)-tand(alpha2)));
        w(:,:,j,i)=p(:,:,j,i).*cosd(alpha2);
        
        zetaprofile(:,:,j,i)= Cd(:,:,j,i).*(2.*sqrt(3) + 6.*1/sqrt(3)).*(tand(alpha2)-tand(alpha1));
        deltasprofile(:,:,1,i)=zetaprofile(:,:,1,i).*0.5.*Vbar.^2./Tstageisen(:,:,1);
        if j < stages
            deltasprofile(:,:,j+1,i)=zetaprofile(:,:,j+1,i).*0.5.*Vbar.^2./Tstageisen(:,:,j)+deltasprofile(:,:,j,i);
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

endwall1=zeros(n,n,stages,m);
endwall3=zeros(n,n,stages,m);
endwallmid=zeros(n,n,stages,m);

for i=1:m
    endwall1(:,:,1,i) = 0.25.*Cdwall(:,:,1,i).*(V1.^3).*rhostageisen(:,:,1).*2.*pi.*(2.*r).*Cx(1,i)./Tstageisen(:,:,1);
    endwall3(:,:,1,i) = 0.25.*Cdwall(:,:,1,i).*(V3rel.^3).*rhostageisen(:,:,1).*2.*pi.*(2.*r).*Cx(1,i)./Tstageisen(:,:,1);
    endwallmid(:,:,1,i) = 0.5*Cdwall(:,:,1,i).*((V2rel.^3+V2.^3)./2).*rhostageisen(:,:,1).*2.*pi.*(2.*r).*Cx(1,i)./Tstageisen(:,:,1);
    for j=2:stages
        endwall1(:,:,j,i) = endwall1(:,:,j-1,i)+0.25.*Cdwall(:,:,j,i).*(V1.^3).*rhostageisen(:,:,j).*2.*pi.*(2.*r).*Cx(j,i)./Tstageisen(:,:,j);
        endwall3(:,:,j,i) =endwall3(:,:,j-1,i)+ 0.25.*Cdwall(:,:,j,i).*(V3rel.^3).*rhostageisen(:,:,j).*2.*pi.*(2.*r).*Cx(j,i)./Tstageisen(:,:,j);
        endwallmid(:,:,j,i) =endwallmid(:,:,j-1,i)+ 0.5.*Cdwall(:,:,j,i).*((V2rel.^3+V2.^3)./2).*rhostageisen(:,:,j).*2.*pi.*(2.*r).*Cx(j,i)./Tstageisen(:,:,j);
    end
end

deltasendwall = (endwall1+endwall3+endwallmid)./mdot;

%--------------------------------------------------------------------------
%secondary flow calc
alpham=-acotd(0.5.*(cotd(alpha1)+cotd(alpha2)));%nxn

Ysec=zeros(n,n,stages,m);
deltas2nd=zeros(n,n,stages ,m);

for i = 1:m
    
    Ysec(:,:,1,i)=2*(0.375.*0.1336).*Cx(1,i).*(cosd(alpha2).^3).*((tand(alpha1)-tand(alpha2)).^2)./...
        (Hstage(:,:,1).*sqrt(cosd(alpha1)).*cosd(alpham));
    
    deltas2nd(:,:,1,i)=Ysec(:,:,1,i).*(hostage(1)-hstage(:,:,1))./(Tstageisen(:,:,1));
    for j=2:stages
        Ysec(:,:,j,i)=2*(0.375.*0.1336).*Cx(j,i).*(cosd(alpha2).^3).*((tand(alpha1)-tand(alpha2)).^2)./...
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
     
    fracleakagerot = zeros(n,n,stages, m);
    fracleakagestat = zeros(n,n,stages, m);
    deltastip = zeros(n,n,stages, m);
    for i = 1:m
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

volstage = zeros(n,n,stages,m);
for j=1:stages
for i=1:m
    volstage(:,:,j,i) = pi*stagewidth(j,i).*((r+Hstage(:,:,j)./2).^2);  
end
end
vol =sum(volstage,3);
paretovol = vol;
paretoneff = 1-efftotal;

paretovol = reshape(paretovol, 1,[]);
paretoneff = reshape(paretoneff,1,[]);
%reshape will order as: vary phi, vary psi, vary AR
store=zeros(length(paretoneff),5);


for i=1:n*m
    store(n*(i-1)+1:n*(i),1)=phi(1,:);
   
end
psi2=reshape(psi',1,[]);

for k=1:m
    store((n*n*(k-1))+1:n*n*k,2)=psi2;
end

for j=1:m
    store(n*n*(j-1)+1:n*n*(j),3)=repelem(Cx(1,j),n*n);
end


for i = 1:length(paretoneff)%point we're comparing to
    for j = 1:length(paretoneff) %point we're looking at
        store(j,4)= paretoneff(j);
        store(j,5)= paretovol(j);
        if  paretovol(j)~=0 && paretovol(i)~=0 && paretovol(j)>paretovol(i) && paretoneff(j)>paretoneff(i)
            
            store(j,:)=0;
            paretovol(j)=0;
            paretoneff(j)=0;
            
        end
    end
end

x=length(nonzeros(store(:,1)));
store2=zeros(x,5);
for col=1:5
    store2(:,col)=nonzeros(store(:,col));
end
figure(2)
 scatter(paretovol, paretoneff, 5,'DisplayName', num2str(z))
 hold on
 xlabel('Volume')
 ylabel('1 - \eta_{tt}')
 end
end
hold off
