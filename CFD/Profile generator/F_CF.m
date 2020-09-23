%Generate class function of psi

function C = F_CF(psi)

C = (psi.^0.5).*(1-psi);
end