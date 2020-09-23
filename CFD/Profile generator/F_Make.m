%Compile a blade from reduced variables

function [X_L, Z_L, X_U, Z_U, thk] = F_Make(Xi1, Xi2, controls_cam,...
    controls_thk_x, controls_thk_y, Tte, Oval_frac)

Xi1 = deg2rad(Xi1);
Xi2 = deg2rad(Xi2);
[x_cam, y_cam] = GEN_TURNING(Xi1, Xi2, controls_cam);
Gamma = rad2deg(atan(y_cam(end)/x_cam(end)));
s2 = zeros(length(x_cam), 1);
for i=2:length(x_cam)
    s2(i) = s2(i-1)+((x_cam(i)-x_cam(i-1))^2+(y_cam(i)-y_cam(i-1))^2)^0.5;	
end
psi = s2/max(s2);
psi_new = dist_vino(200, 0, 1, 1/2000, 1/500);
x_cam = interp1(psi, x_cam, psi_new);
y_cam = interp1(psi, y_cam, psi_new);
norm = calc_norm(x_cam, y_cam);
y_cam = transpose(y_cam);
x_cam = transpose(x_cam);
thk = F_TF_bez(controls_thk_x, controls_thk_y, Tte, Oval_frac);
thk = transpose(thk);
Z_U = y_cam+thk.*norm(:, 2);
Z_L = y_cam-thk.*norm(:, 2);
X_U = x_cam+thk.*norm(:, 1);
X_L = x_cam-thk.*norm(:, 1);

end 