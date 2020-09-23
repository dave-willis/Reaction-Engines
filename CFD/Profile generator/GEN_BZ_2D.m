%Generate x,y of a bezier from control points

function [x, y, CPs] = GEN_BZ_2D(controls_x, controls_y)

CPs = zeros(length(controls_x),2);

for i=1:length(controls_x)
    CPs(i, :) = [controls_x(i), controls_y(i)];
end

steps = 100;
controlPoints = CPs;
x = zeros(steps,1);
y = zeros(steps,1);
bez_range = bezier_curve_range(steps, controlPoints);
for i=1:steps;
    x(i) = bez_range(i,1);
    y(i) = bez_range(i,2);
end
end
