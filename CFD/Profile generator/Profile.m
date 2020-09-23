%Function based on profgen.f (JDD) - GENERATES SIMPLE PROFILES
%EXCUSE THE WEIRD FORMULATION, THIS IS CONVERSION FROM FORTRAN
%USER HARD CODED VALUES BELOW - COULD BE TAKEN OUTSIDE FUNCTION CJC

function [Z, RTH, te_index] = Profile(X1, X2, TKTE, Cx)

controls_cam = [0 -0.5 -0.8];
Rle = 0.05;
Tte = TKTE/Cx;
Beta_te = 4.0;
Oval_frac = 0.3;
controls_thk_x = [0.0 0.5 1.0];
controls_thk_y = [(2*Rle)^0.5*(1-Oval_frac) 0.25 tan(deg2rad(Beta_te))+Tte/2];

[XLIN, YLIN, XUIN, YUIN, thk] = F_Make(X1, X2, controls_cam, controls_thk_x,...
    controls_thk_y, Tte, Oval_frac);

XIN = (XUIN+XLIN)/2;
YIN = (YUIN+YLIN)/2;
XScale = max([max(XUIN), max(XLIN)]);
XUIN = XUIN/XScale;
XLIN = XLIN/XScale;
XIN = XIN/XScale;

yoffset = mean(YIN);
YUIN=(YUIN-yoffset)/XScale;
YLIN=(YLIN-yoffset)/XScale;
YIN = (YIN-yoffset)/XScale;

Xnew = zeros(length(XUIN)+length(XLIN)-1, 1);
Xnew(1:200) = flip(XUIN);
Xnew(200:end) = XLIN;

Ynew = zeros(length(YUIN)+length(YLIN)-1, 1);
Ynew(1:200)= flip(YUIN);
Ynew(200:end) = YLIN;

Snew = zeros(length(Ynew), 1);
for i=2:length(Xnew)
    Snew(i) = Snew(i-1)+((Xnew(i)-Xnew(i-1))^2+(Ynew(i)-Ynew(i-1))^2)^0.5;
end

for i=2:length(Xnew)
    if Xnew(i) == min(Xnew)
        Sle = Snew(i);
        Xle = Xnew(i);
        break
    end
end

points = 501;
XUIN = flip(interp1(Snew, Xnew, linspace(0, Sle, points)));
YUIN = flip(interp1(Snew, Ynew, linspace(0, Sle, points)));
XLIN = interp1(Snew, Xnew, linspace(Sle, Snew(end), points));
YLIN = interp1(Snew, Ynew, linspace(Sle, Snew(end), points));

XScale = max([max(XUIN), max(XLIN)])-Xle;
XUIN = (XUIN-Xle)/XScale;
XLIN = (XLIN-Xle)/XScale;
XIN = (XIN-Xle)/XScale;

YUIN=(YUIN-yoffset)/XScale;
YLIN=(YLIN-yoffset)/XScale;
YIN =(YIN-yoffset)/XScale;

UTE_m = (XUIN(end-1)-XUIN(end))/(YUIN(end)-YUIN(end-1));
LTE_m = (XLIN(end-1)-XLIN(end))/(YLIN(end)-YLIN(end-1));
TE_circ_cx = (YUIN(end)-YLIN(end)+LTE_m*XLIN(end)-UTE_m*XUIN(end))/(LTE_m-UTE_m);
TE_circ_cy = YUIN(end)+UTE_m*(TE_circ_cx-XUIN(end));
TE_circ_r = ((XLIN(end)-TE_circ_cx)^2+(YLIN(end)-TE_circ_cy)^2)^0.5;
U_theta = acos((XUIN(end)-TE_circ_cx)/TE_circ_r);
L_theta = -acos((XLIN(end)-TE_circ_cx)/TE_circ_r);
TE_n_points = 100;
TE_theta = linspace(0.95*U_theta, L_theta, TE_n_points);
TEx = TE_circ_cx - TE_circ_r*cos(pi-TE_theta);
TEy = TE_circ_cy + TE_circ_r*sin(pi-TE_theta);

miny = min([min(YUIN) min(YLIN) min(TEy)]);

XLIN = XLIN(1:end)*Cx;
YLIN = (YLIN(1:end)-miny)*Cx;
XUIN = XUIN(1:end)*Cx;
YUIN = (YUIN(1:end)-miny)*Cx;
TEx = TEx*Cx;
TEy = (TEy-miny)*Cx;

if X1 < X2
    Z = [XLIN flip(TEx) flip(XUIN)];
    RTH = [YLIN flip(TEy) flip(YUIN)];
else
    Z = [XUIN TEx flip(XLIN)];
    RTH = [YUIN TEy flip(YLIN)];
end
te_index = round(length(Z)/2,0);
end