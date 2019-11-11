"""Calculations for evaluating turbine performance"""

from Stage import turbine, angles, spline, optimise
from Profile import b2b_plot, annulus, b2b_variable
import numpy as np
import time
from multiprocessing import cpu_count
from multiprocessing.pool import Pool
from scipy.optimize import minimize
import sys

Po1 = 145*10**5
To1 = 950
W = 17*10**6
mdot = 16
Omega = 6782*2*np.pi/60
t = 0.00015
g = 0.00015

cp = 5187
del_ho = W/mdot
To3 = To1-del_ho/cp

##20 stages
#phi = 0.3
#psi = 0.8
#Lambda = 0.5
#AR = 1.5
#ptc = -1.0
#n = 20
##10 stages
#phi = 0.3
#psi = 0.8
#Lambda = 0.5
#AR = 1.0
#ptc = 1.056
#n = 10
##5 stages
#phi = 0.4
#psi = 1.1
#Lambda = 0.5
#AR = 1.0
#ptc = 1.0
#n = 5

phi_lim = (0.2, 1.0)
psi_lim = (0.5, 2.5)
Lam_lim = (0, 1)
AR_lim = (1, 2)
ptc_lim = (0.7, 1.3)
dh_lim = (1, 5)

#10 stages
phi = [0.281, 0.307]
psi = [0.833, 0.892]
Lambda = [0.5, 0.5]
AR = [1.0, 1.0]
ptc = [1.1, 1.1]
n = 10
ain = 0
dho = [1.0, 1.1]

#phi = [0.35, 0.35]
#psi = [1,1, 1.1]
#Lambda = [0.5, 0.5]
#AR = [1.0, 1.0]
#ptc = [1.3, 1.3]
#n = 10
#dho = [1.0, 1.0]

calcs = ''
plot = 'yes'
save = ''
start_time = time.time()
result = turbine(Po1, To1, mdot, Omega, W, t, g, phi, psi, Lambda, AR, dho, n, ptc, ain)
print('Time: {} s'.format(time.time()-start_time))
#print('Angles [a1,a2,b2,a3,b3]=', np.round(result[10],2))
#print('Chords [Cxst,Cxro]=', np.round([[i[3],i[4]] for i in result[5]],6))
print('Work = {} W'.format(result[1]))
print('Efficiency = {}'.format(result[0]))
print('Mass = {} kg'.format(result[2]))
print('Volume = {} m^3'.format(result[3]))
print('No. Blades = {}'.format(int(result[6])))
print('Axial force on rotor = {} N'.format(result[13]))
print(result[12])
print('')
if plot == 'yes':
    b2b_variable(result)
#    b2b_plot(result)
#    annulus(result)


if calcs == 'brute force':
    start_time = time.time()
                    
    #Create arrays of all the loading combinations 
    phi_set = np.arange(0.3, 0.9, 0.1)
    psi_set = np.arange(0.8, 2, 0.1)
    Lambda = [0.5, 0.5]
    dh_set = np.arange(1, 1.2, 0.02)
    
    phis = np.zeros([len(phi_set)**2,2])
    psis = np.zeros([len(psi_set)**2,2])
    dhs = np.zeros([len(dh_set),2])
    
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
                        
    def turb_calcs(var):
        To1, Po1, n, phi, psi, dho, AR, a1i = [i for i in var]
        eff = turbine(Po1, To1, mdot, Omega, phi, psi, Lambda, AR, W*n/(n+1), dho, n, t, g, ptc, a1i)[0]
        return [eff, n, phi, psi, Lambda, dho, AR, a1i]
    
    p = Pool(processes=cpu_count()-2)
    results = p.map(turb_calcs, variables)
    p.close()
    p.join()
    
    hours = int((time.time()-start_time)/3600)
    minutes = int((time.time()-start_time-hours*3600)/60)
    seconds = time.time()-start_time-hours*3600-minutes*60
    
    print('Turbine calculations done in', hours, 'h', minutes, 'm', seconds, 's')

    results.sort(key=lambda x: (x[1], x[0]))    
    
if calcs == 'opt1':
    
    start = time.time()
    
    phi, psi, Lambda, AR, dho = optimise(result)
    
    print("Optimizer time: {}".format(time.time()-start))
    
    turbine_data = turbine(Po1, To1, mdot, Omega, W, t, g, phi, psi, Lambda, AR, dho, n, ptc, ain)
    
    print("Optimum efficiency = {}".format(turbine_data[0]))
    print("phi = {}".format(np.round(phi, 4)))
    print("psi = {}".format(np.round(psi, 4)))
    print("Lambda = {}".format(np.round(Lambda, 4)))
    print("AR = {}".format(np.round(AR, 4)))
    print("dho = {}".format(np.round(dho, 4)))
    print("ain = {}".format(np.round(turbine_data[11][14], 4)))
    print("Max angle = {}".format(np.round(turbine_data[12][1], 4)))
    print('Mass = {} kg'.format(result[2]))
    print('No. Blades = {}'.format(int(result[6])))
    print('Axial force on rotor = {} N'.format(result[13]))

    if plot == 'opt':
        b2b_plot(turbine_data)
        annulus(turbine_data)
        
if calcs == 'optall':
    
    start=time.time()
    #Create arrays of all the loading combinations 
    phi_set = np.arange(0.2, 0.9, 0.3)
    psi_set = np.arange(0.6, 2, 0.4)
    Lambda_set = np.arange(0.2, 0.9, 0.3)
    dh_set = np.arange(1, 1.2, 0.1)
    
    phis = np.zeros([len(phi_set)**2,2])
    psis = np.zeros([len(psi_set)**2,2])
    Lambdas = np.zeros([len(Lambda_set)**2,2])
    dhs = np.zeros([len(dh_set)**2,2])
    
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
        
    n_Lambda = 0
    for i in Lambda_set:
        for j in Lambda_set:
            Lambdas[n_Lambda] = [i, j]
            n_Lambda += 1

    n_dh = 0
    for i in dh_set:
        for j in dh_set:
            dhs[n_dh] = [j, i]
            n_dh += 1
    
    variables = []
    i = 1
    for phi in phis:
        for psi in psis:
            for Lambda in Lambdas:
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
                    for AR in np.arange(0.9, 2.0, 0.5):
                        ARs = [AR, AR]
                        number = [i, i]
                        variables.append([phi, psi, Lambda, ARs, dho, number])
                        i += 1
   
    def turbine_opts(args):
    
        ph1, ph2, ps1, ps2, L1, L2, AR1, AR2, dh1, dh2 = [i for i in args]
        phi = [ph1, ph2]
        psi = [ps1, ps2]
        Lambda = [L1, L2]
        dho = [dh1, dh2]
        AR = [AR1, AR2]
             
        eff = turbine(Po1, To1, mdot, Omega, W, t, g, phi, psi, Lambda, AR, dho, n, ptc, ain)[0]
        
        return -eff
    
    def con_a2(args):
    
        ph1, ph2, ps1, ps2, L1, L2, AR1, AR2, dh1, dh2 = [i for i in args]
        
        phi = [ph1, ph2]
        psi = [ps1, ps2]
        Lambda = [L1, L2]
        
        phi = spline(n, phi)
        psi = spline(n, psi)
        Lambda = spline(n, Lambda)
        
        a2_max = 0
        a1 = ain
        for i in range(len(phi)):
            ang_check = angles(phi[i], psi[i], Lambda[i], a1)
            a2 = np.degrees(ang_check[1])
            a1 = ang_check[3] 
            if abs(a2) > a2_max:
                a2_max = abs(a2)
        
        return 73-a2_max
    
    def con_b3(args):
        
        ph1, ph2, ps1, ps2, L1, L2, AR1, AR2, dh1, dh2 = [i for i in args]
        
        phi = [ph1, ph2]
        psi = [ps1, ps2]
        Lambda = [L1, L2]
        
        phi = spline(n, phi)
        psi = spline(n, psi)
        Lambda = spline(n, Lambda)
        
        b3_max = 0
        a1 = ain
        for i in range(len(phi)):
            ang_check = angles(phi[i], psi[i], Lambda[i], a1)
            b3 = np.degrees(ang_check[4])
            a1 = ang_check[3] 
            if  abs(b3) > b3_max:
                b3_max = abs(b3)
        
        return 73-b3_max
                         
    def optimize(vals):
        
        phi01, psi01, Lam01, AR01, dho01, num = [i[0] for i in vals]
        phi02, psi02, Lam02, AR02, dho02, num = [i[1] for i in vals]
        
        x0 = [phi01, phi02, psi01, psi02, Lam01, Lam02, AR01, AR02, dho01, dho02]
        bnds = (phi_lim, phi_lim, psi_lim, psi_lim, Lam_lim, Lam_lim, AR_lim, AR_lim, dh_lim, dh_lim)
        cons = cons = ({'type': 'ineq', 'fun': con_a2}, {'type': 'ineq', 'fun': con_b3})
        res = minimize(turbine_opts, x0, method='SLSQP', bounds=bnds, constraints=cons)
        print(-res['fun'], num)
        
        return res
    
    p = Pool(processes=cpu_count())
    try:
        results = p.map_async(optimize, variables).get(9999999)
    except KeyboardInterrupt:
        p.terminate()
        sys.exit('KeyboardInterrupt')
    p.close()
    p.join()
    
    hours = int((time.time()-start_time)/3600)
    minutes = int((time.time()-start_time-hours*3600)/60)
    seconds = time.time()-start_time-hours*3600-minutes*60
    
    print('Turbine calculations done in', hours, 'h', minutes, 'm', seconds, 's')
    
    results.sort(key=lambda x: -x['fun'])
    
    for i in results[-9:]:
        print('Efficiency = {}'.format(np.round(-i['fun'], 4)))
        print('phi = {}'.format(np.round([i['x'][0], i['x'][1]],3)))
        print('psi = {}'.format(np.round([i['x'][2], i['x'][3]],3)))
        print('Lambda = {}'.format(np.round([i['x'][4], i['x'][5]],3)))
        print('AR = {}'.format(np.round([i['x'][6], i['x'][7]],3)))
        print('dho = {}'.format(np.round([i['x'][8], i['x'][9]],3)))
        print('')
        
#if calcs == 'sensitivity':
#    
#    base_eff = turbine(Po1, To1, mdot, Omega, W, t, g, phi, psi, Lambda, AR, dho, n, ptc, ain)
    
    