function [w, Kt, dr] = Flow_guess()
%FLOW_GUESS Generate flow conditions and initial area guess
%   Use the midspan ideal design to find the rotational speed of the
%   turbine and constants for the tangental velocity distribution.  Then
%   finds an initial guess for the annulus area assuming constant Vx and rho.
 
global phi psi L T01rm P01rm M3r rm p g R cp H_init CRT

w = zeros(1,2);
Kt = zeros(1,2);
dr = zeros(3,1);

tm1 = (1-L-(psi/2))/phi;
tm2 = (1-L+(psi/2))/phi;
sm1 = sec(atan(tm1));
sm2 = sec(atan(tm2));
Urm = sqrt(((M3r^2)*g*R*T01rm)/(((phi*sm2)^2)+((g*R*(M3r^2)/(2*cp))*((2*psi)+((sm1*phi)^2)))));
Vx = Urm*phi;
if CRT == 0
    w(1) = 0;
    w(2) = Urm/rm;
elseif CRT == 1
    w(1) = -Urm/(2*rm);
    w(2) = Urm/(2*rm);
else
    error('CRT variable must equal 1 or 0')
end

Vtm1 = (Vx*tm1) + (rm*w(1))
Vtm2 = (Vx*tm2) + (rm*w(1))
% Kt(1) = ((Vx*tm1) + (w(1)*rm))/(rm^p); %A1
% Kt(2) = (Vx*tm2) + (w(1)*rm)/(rm^p); %A2
if p == 0
    Kt(1) = (Vx*tm1)+(w(1)*rm); %A1
    Kt(2) = (Vx*tm2)+(w(1)*rm); %A2
else
    Kt(1) = ((Vtm1)+(Vtm2))/(2*(rm^p));
    Kt(2) = (-(Vtm1)+(Vtm2))*rm/2;
%     if CRT == 1
%         Kt(2) = Kt(2) + (w(1)*(rm^2));
%     end
end
T01 = T01rm + ((((Vx*tm1) + (w(1)*rm))^2)-((Vx*tm1)^2))/(2*cp);
T02 = T01 - abs((rm*w(1)*Vx*(tm2-tm1))/cp);
T03 = T02 - abs((rm*w(2)*Vx*(tm2-tm1))/cp); 
% rho1 = (P01rm/(R*(T01rm^(g/(g-1)))))*((T01 -(((Kt(1)^2)+(Vx^2))/(2*cp)))^(1/(g-1)));
% rho2 = (P01rm/(R*(T01rm^(g/(g-1)))))*((T02 -(((Kt(2)^2)+(Vx^2))/(2*cp)))^(1/(g-1)));
% rho3 = (P01rm/(R*(T01rm^(g/(g-1)))))*((T03 -(((Kt(1)^2)+(Vx^2))/(2*cp)))^(1/(g-1)));
if p == 0
    rho1 = (P01rm/(R*(T01rm^(g/(g-1)))))*((T01 -(((Kt(1)^2)+(Vx^2))/(2*cp)))^(1/(g-1)));
    rho2 = (P01rm/(R*(T01rm^(g/(g-1)))))*((T02 -(((Kt(2)^2)+(Vx^2))/(2*cp)))^(1/(g-1)));
    rho3 = (P01rm/(R*(T01rm^(g/(g-1)))))*((T03 -(((Kt(1)^2)+(Vx^2))/(2*cp)))^(1/(g-1)));
else
    rho1 = (P01rm/(R*(T01rm^(g/(g-1)))))*((T01 -(((((Kt(1)*(rm^p))-(Kt(2)/rm)+(w(1)*rm))^2)+(Vx^2))/(2*cp)))^(1/(g-1)));
    rho2 = (P01rm/(R*(T01rm^(g/(g-1)))))*((T02 -(((((Kt(1)*(rm^p))+(Kt(2)/rm)+(w(1)*rm))^2)+(Vx^2))/(2*cp)))^(1/(g-1)));
    rho3 = (P01rm/(R*(T01rm^(g/(g-1)))))*((T03 -(((((Kt(1)*(rm^p))-(Kt(2)/rm)+(w(1)*rm))^2)+(Vx^2))/(2*cp)))^(1/(g-1)));
end
% if p == -1 || p == 0
%     Kt(3) = 0;
% else
%     Kt(3) = 1;
% end
mdot = H_init*Vx*rho1;
dr(1) = H_init/2;
dr(2) = 0.5*mdot/(Vx*rho2);
dr(3) = 0.5*mdot/(Vx*rho3);
%w = tm1;
%Kt = tm2;
end

