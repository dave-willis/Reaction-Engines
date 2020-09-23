%Return thickness contributions due to parameters

function TF = F_TF_bez(controls_cam_x, controls_cam_y, TE, Oval_frac)

psi = dist_vino(200, 0, 1, 1/2000, 1/500);
[X, S, cp] = GEN_BZ_2D(controls_cam_x, controls_cam_y);
S = interp1(X, S, psi);
S = S+Oval_frac*S(1).*(1-psi).^18;
C = F_CF(psi);

TF = S.*C+psi.*TE/2;
end