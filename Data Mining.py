"""Calculations for evaluating turbine performance"""

from Stage import turbine, angles, spline
from Profile import b2b_plot, annulus, b2b_variable
import numpy as np
import time
from multiprocessing import cpu_count
from multiprocessing.pool import Pool
from scipy.optimize import minimize

Po1 = 145*10**5
To1 = 950
W = 17*10**6
mdot = 16
Omega = 6782*2*np.pi/60
t = 0.0003
g = 0.0003

cp = 5187
del_ho = W/mdot
To3 = To1-del_ho/cp

##20 stages
#phi = 0.3
#psi = 0.8
#Lambda = 0.5
#AR = 1.5
#ptoC = -1.0
#n = 20
##10 stages
#phi = 0.3
#psi = 0.8
#Lambda = 0.5
#AR = 1.0
#ptoC = 1.056
#n = 10
##5 stages
#phi = 0.4
#psi = 1.1
#Lambda = 0.5
#AR = 1.0
#ptoC = 1.0
#n = 5

#10 stages
phi = [0.317, 0.335]
psi = [1.065, 1.067]
Lambda = [0.496, 0.527]
AR = [1.0, 1.0]
ptoC = [1.067, 1.058]
n = 10
dho = [1.002, 1.008]

calcs = ''
plot = 'yes'
save = ''
start_time = time.time()
result = turbine(Po1, To1, mdot, Omega, phi, psi, Lambda, AR, W, dho, n, t, g, ptoC)
print('Time:{} s'.format(time.time()-start_time))
#print('Angles [a1,a2,b2,a3,b3]=', np.round(result[10],2))
#print('Chords [Cxst,Cxro]=', np.round([[i[3],i[4]] for i in result[5]],6))
print('Work = {} W'.format(result[1]))
print('Efficiency = {}'.format(result[0]))
print('Mass = {} kg'.format(result[2]))
print('Volume = {} m^3'.format(result[3]))
print('No. Blades = {}'.format(int(result[6])))
print('')
if plot == 'yes':
    b2b_variable(result)
#    b2b_plot(result)
#    annulus(result)


if calcs == 'brute force':
    start_time = time.time()
    
#    #Create arrays of all the loading combinations 
#    phi_set = np.arange(0.25, 0.75, 0.2)
#    psi_set = np.arange(0.75, 1.75, 0.2)
#    Lambda = [0.5, 0.5, 0.5]
#    dh_set = np.arange(1, 2, 0.25)
#    
#    phis = np.empty([len(phi_set)**3,3])
#    psis = np.empty([len(psi_set)**3,3])
#    dhs = np.empty([len(dh_set)**3-len(dh_set)+1,3])
#    
#    n_phi = 0
#    for i in phi_set:
#        for j in phi_set:
#            for k in phi_set:
#                phis[n_phi] = [i, j, k]
#                n_phi += 1
#    n_psi = 0
#    for i in psi_set:
#        for j in psi_set:
#            for k in psi_set:
#                psis[n_psi] = [i, j, k]
#                n_psi += 1
#    n_dh = 1
#    dhs[0] = [1, 1, 1]
#    for i in dh_set:
#        for j in dh_set:
#            for k in dh_set:
#                if i != j or i != k or j != k:
#                    dhs[n_dh] = [i, j, k]
#                    n_dh += 1
                    
    #Create arrays of all the loading combinations 
    phi_set = np.arange(0.3, 0.9, 0.1)
    psi_set = np.arange(0.8, 2, 0.1)
    Lambda = [0.5, 0.5]
    dh_set = np.arange(1, 1.2, 0.02)
    
    phis = np.empty([len(phi_set)**2,2])
    psis = np.empty([len(psi_set)**2,2])
    dhs = np.empty([len(dh_set),2])
    
    n_phi = 0
    for i in phi_set:
        for j in phi_set:
            phis[n_phi] = [i, j]
            n_phi += 1
    n_psi = 0
    for i in psi_set:
        for j in psi_set:
            psis[n_psi] = [i, j]
            n_psi += 1
    n_dh = 0
    for i in dh_set:
        dhs[n_dh] = [1, i]
        n_dh += 1
    
    variables = []
    for n in np.arange(4, 10):
        To1 = 950-del_ho/(cp*(n+1))
        Po1 = (145*10**5)*(To1/950)**(1.6625/(0.8*0.6625))
        for phi in phis:
            for psi in psis:
                ph = spline(n, phi)
                ps = spline(n, psi)
                L = spline(n, Lambda)
                ang_1 = 0 #Initialise the inlet angle
                max_ang = False 
                for l in range(len(ph)):
                    ang_check = angles(ph[l], ps[l], L[l], ang_1)
                    ang_2 = np.degrees(ang_check[1])
                    ang_3 = np.degrees(ang_check[4])
                    ang_1 = ang_check[3]
                    if abs(ang_2) > 73 or abs(ang_3) > 73:
                        max_ang = True
                if max_ang:
                    continue
                for dho in dhs:
                    for AR in np.arange(1.0, 2.0, 0.2):
                        for a1i in np.arange(0, 20, 2):
                            variables.append([To1, Po1, n, phi, psi, dho, AR, a1i])
                        
    def turbine_calcs(var):
        To1, Po1, n, phi, psi, dho, AR, a1i = [i for i in var]
        eff = turbine(Po1, To1, mdot, Omega, phi, psi, Lambda, AR, W*n/(n+1), dho, n, t, g, ptoC, a1i)[0]
        return [eff, n, phi, psi, Lambda, dho, AR, a1i]
    
    p = Pool(processes=cpu_count()-2)
    results = p.map(turbine_calcs, variables)
    p.close()
    p.join()
    
    hours = int((time.time()-start_time)/3600)
    minutes = int((time.time()-start_time-hours*3600)/60)
    seconds = time.time()-start_time-hours*3600-minutes*60
    
    print('Turbine calculations done in', hours, 'h', minutes, 'm', seconds, 's')

    results.sort(key=lambda x: (x[1], x[0]))
    
if calcs == 'opt':
    
    Po1 = 145*10**5
    To1 = 950
    W = 17*10**6
    mdot = 16
    Omega = 6782*2*np.pi/60
    t = 0.0003
    g = 0.0003
    
    n = 5
    
    phi0 = 0.35
    psi0 = 1.1
    Lam0 = 0.5
    AR0 = 1.0
    ptoC0 = 1.0
    dho0 = 1.0
         
    phi_lim = (0.2, 1.0)
    psi_lim = (0.5, 2.5)
    Lam_lim = (0, 1)
    AR_lim = (1, 2)
    ptoC_lim = (0.7, 1.3)
    dh_lim = (1, 5)
    
    def turbine_calcs(args):
        
        ph1, ph2, ps1, ps2, L1, L2, AR1, AR2, ptoC1, ptoC2, dh1, dh2 = [i for i in args]
        phi = [ph1, ph2]
        psi = [ps1, ps2]
        Lambda = [L1, L2]
        dho = [dh1, dh2]
        AR = [AR1, AR2]
        ptoC = (ptoC1, ptoC2)
             
        eff = turbine(Po1, To1, mdot, Omega, phi, psi, Lambda, AR, W, dho, n, t, g, ptoC)[0]
        
        return -eff
    
    def constraint_a2(args):
        
        ph1, ph2, ps1, ps2, L1, L2, AR1, AR2, ptoC1, ptoC2, dh1, dh2 = [i for i in args]
        
        phi = [ph1, ph2]
        psi = [ps1, ps2]
        Lambda = [L1, L2]
        
        phi = spline(n, phi)
        psi = spline(n, psi)
        Lambda = spline(n, Lambda)
        
        a2_max = 0
        a1 = 0
        for i in range(len(phi)):
            ang_check = angles(phi[i], psi[i], Lambda[i], a1)
            a2 = np.degrees(ang_check[1])
            a1 = ang_check[3] 
            if abs(a2) > a2_max:
                a2_max = abs(a2)
        
        return 73-a2_max

    def constraint_b3(args):
        
        ph1, ph2, ps1, ps2, L1, L2, AR1, AR2, ptoC1, ptoC2, dh1, dh2 = [i for i in args]
        
        phi = [ph1, ph2]
        psi = [ps1, ps2]
        Lambda = [L1, L2]
        
        phi = spline(n, phi)
        psi = spline(n, psi)
        Lambda = spline(n, Lambda)
        
        b3_max = 0
        a1 = 0
        for i in range(len(phi)):
            ang_check = angles(phi[i], psi[i], Lambda[i], a1)
            b3 = np.degrees(ang_check[4])
            a1 = ang_check[3] 
            if  abs(b3) > b3_max:
                b3_max = abs(b3)
        
        return 73-b3_max    
        
    x0 = [phi0, phi0, psi0, psi0, Lam0, Lam0, AR0, AR0, ptoC0, ptoC0, dho0, dho0]
    bnds = (phi_lim, phi_lim, psi_lim, psi_lim, Lam_lim, Lam_lim, AR_lim, AR_lim, ptoC_lim, ptoC_lim, dh_lim, dh_lim)
    cons = cons = ({'type': 'ineq', 'fun': constraint_a2}, {'type': 'ineq', 'fun': constraint_b3})
    
    res = minimize(turbine_calcs, x0, method='SLSQP', bounds=bnds, constraints = cons)
    
    phi = [res['x'][0], res['x'][1]]
    psi = [res['x'][2], res['x'][3]]
    Lambda = [res['x'][4], res['x'][5]]
    AR = [res['x'][6], res['x'][7]]
    ptoC = [res['x'][8], res['x'][9]]
    dho = [res['x'][10], res['x'][11]]
    
    print("Optimum efficiency =", -res['fun'])
    print("phi =", np.round(phi, 4))
    print("psi =", np.round(psi, 4))
    print("Lambda =", np.round(Lambda, 4))
    print("AR =", np.round(AR, 4))
    print("ptc =", np.round(ptoC, 4))
    print("dho =", np.round(dho, 4))
    
    result = turbine(Po1, To1, mdot, Omega, phi, psi, Lambda, AR, W, dho, n, t, g, ptoC)

    if plot == 'opt':
        b2b_plot(result)
        annulus(result)
            