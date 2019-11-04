function out = stresscalc(omega, r, H, T )
%Calculates stresses in disk
%Physical Properties: Inconel @ 1000K
rho=8440;
alpha=15*10^-6;
E=160.6*10^9;
G=60*10^9;
nu=0.34;

%--------------------------------------------------------------------------
%Centrifugal Stress, assuming untapered blade
%sigmac=((3+nu)./8)*(rho.*omega.^2.*(r-H./2).^2);

A=((3+nu)./8).*rho.*(omega.^2).*((r-H./2).^2); %'A' term in Lame's Equations

sigmarrdisk=((3+nu)./4).*rho.*(omega.^2).*(r.*H); %from Lame's Equations, I haven't used it because,
%strain formula only really needs the 'A' term

sigmarrbladeroot=rho.*omega.^2.*(r.*H); %stress from spinning mass of blade hanging off blade root

%Thermal Stress = none, assuming entire disk at same temp (unpurged disk)

figure(10)
plot(1:length(T),sigmarrbladeroot )
hold on
xlabel('Stage')
ylabel('\sigma_{rr}')

%strains
strnthetahubcentrifugal=(1/E).*(A.*(1-nu));
strnthetatemp=alpha.*T;
rt=r+H./2; %tip radius
rh=r-H./2; %hib radius

 
uhubcentrifugal=strnthetahubcentrifugal.*(rh);
utipcentrifugal=(rho.*omega.^2./(2*E)).*((2.*(rt.^3)./3) - rt.^2.*rh+(rh.^3)./3);
utemp=rt.*strnthetatemp;
utip=utemp+utipcentrifugal;

%now output whatever you want
figure(11)
plot(1:length(T),utip)
hold on
xlabel('Stage')
ylabel('Total Tip Stretching')
figure(12)
plot(1:length(T),utipcentrifugal)
out = [sigmarrbladeroot; utip];
end

