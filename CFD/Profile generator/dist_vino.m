%Vinokur distribution

function xs = dist_vino(n, xst, xen, s1, s2)

s1 = s1/(xen-xst);
s2 = s2/(xen-xst);
eta = linspace(0,n-1,n);
a = sqrt(s2/s1);
b = 1/((n-1)*sqrt(s1*s2));
nm1 = n-1;
trans = @(delta) sinh(delta)/delta - b;
x0 = [0.00001 100];
delta = fzero(trans, x0);
u = 0.5*(1+tanh(delta*(eta/nm1-0.5))/tanh(delta*0.5));
s = u./(a+(1-a)*u);
xs = xst+(xen-xst)*s;
end