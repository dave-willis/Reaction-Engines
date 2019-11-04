clear
hold off
%calls subroutines for the efficiency and stress calculations
for z = 0:1
effcalc = turbineefficiencyphi(0.23,1,20,[0.4:0.2:3],z,0.0005,1);
r=effcalc(1,1);
T=effcalc(:,5)';
H=effcalc(:,6)';
stresscalc(6782*2*pi/60,r,H,T-400)
hold on
end
hold off

