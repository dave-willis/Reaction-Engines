%Generate x,y of a bezier from control points

function [xcam, ycam] = GEN_TURNING(Xi1, Xi2, controls)

spacings = linspace(0,1,length(controls)+2);
spacings = spacings(2:end-1);

CPs = zeros(length(controls)+2, 2);
CPs(1,:) = [0 0];
for i=1:length(controls)
    CPs(i+1,:) = [spacings(i)+min(spacings(i),1-spacings(i))*controls(i),...
        spacings(i)-min(spacings(i),1-spacings(i))*controls(i)];
end
CPs(end,:) = [1 1];

steps = 100;
controlPoints = CPs;

x = zeros(steps,1);
y = zeros(steps,1);
bez_range = bezier_curve_range(steps, controlPoints);
for i=1:steps;
    x(i) = bez_range(i,1);
    y(i) = bez_range(i,2);
end

xcam = zeros(steps,1);
ycam = zeros(steps,1);

DXi = Xi2-Xi1;
for i=2:steps
    Xi = Xi1+DXi*(y(i)+y(i-1))*0.5;
    xcam(i) = xcam(i-1)+(x(i)-x(i-1))*cos(Xi);
    ycam(i) = ycam(i-1)+(x(i)-x(i-1))*sin(Xi);
end

DX = xcam(end);
ycam = ycam/DX;
xcam = xcam/DX;
end