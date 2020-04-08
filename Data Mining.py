"""Calculations for evaluating turbine performance"""

from turbine import turbine, angles, spline, optimise, free_vortex
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
t = 0.0003
g = 0.0003
cp = 5187
gas = 'He'
del_ho = W/mdot
To3 = To1-del_ho/cp

# # n = 5
# phi = [0.28, 0.32]
# psi = [0.83, 1.08]
# Lambda = [0.50, 0.52]
# AR = [0.48, 0.62]
# ptc = [1.1, 1.1]
# n = 5
# ain = 0
# dho = [1.00, 1.11]

# n = 10
phi = [0.27, 0.29]
psi = [0.78, 0.90]
Lambda = [0.50, 0.51]
AR = [0.64, 0.80]
ptc = [1.1, 1.1]
n = 10
ain = 0
dho = [1.0, 1.11]

# # n = 15
# phi = [0.27, 0.28]
# psi = [0.73, 0.81]
# Lambda = [0.49, 0.51]
# AR = [0.76, 0.92]
# ptc = [1.1, 1.1]
# n = 15
# ain = 0
# dho = [1.0, 1.11]

# # n = 20
# phi = [0.26, 0.27]
# psi = [0.69, 0.74]
# Lambda = [0.49, 0.51]
# AR = [0.85, 1.01]
# ptc = [1.1, 1.1]
# n = 20
# ain = 0
# dho = [1.0, 1.11]

# #Properties for representative real turbine
# Po1 = 5*10**5
# To1 = 1200
# W = 8*10**6
# mdot = 30
# Omega = 800
# t = 0.001
# g = 0.0002
# cp = 1200
# gas = 'A1'
# del_ho = W/mdot
# To3 = To1-del_ho/cp

# phi = [0.55, 0.65]
# psi = [0.9, 0.95]
# Lambda = [0.5, 0.5]
# AR = [2, 2]
# ptc = [1.1, 1.1]
# n = 1
# dho = [1.0, 1.0]

phi_lim = (0.1, 1.5)
psi_lim = (0.4, 3.0)
Lam_lim = (0, 1)
AR_lim = (0.1, 5)
ptc_lim = (0.7, 1.5)
dh_lim = (1, 5)

calcs = ''
plot = 'yes'
save = ''
save_geom = True
start_time = time.time()
result = turbine(Po1, To1, mdot, Omega, W, t, g, phi, psi, Lambda, AR, dho, n, ptc, ain, gas)
print('Time: {} s'.format(time.time()-start_time))
# print('Angles [a1,a2,b2,a3,b3]=', np.round(result[10],2))
# print('Chords [Cxst,Cxro]=', np.round([[i[4],i[5]] for i in result[5]],6))
print('Work = {} W'.format(result[1]))
print('Efficiency = {}'.format(result[0]))
print('Mass = {} kg'.format(result[2]))
print('No. Blades = {}'.format(int(result[6])))
# print('Axial force on rotor = {} N'.format(result[13]))
# print('Average Re = {}'.format(result[14]))
# print('Cold-stat, cold-rot, warm-rot, hot-rot:', result[15])
print(result[12])
print('')
if 1 == 2:
    
    dims = result[5]
    
    d_hubin = [2*i[9] for i in dims]
    d_hubout = [2*(i[0]-i[2]/2) for i in dims]
    d_casout = [2*i[8] for i in dims]
    d_casin = [2*(i[0]+i[2]/2) for i in dims]
    length = [1.5*(i[4]+i[5]) for i in dims]
    
    np.savez('Dims 10', d_hubin, d_hubout, d_casin, d_casout, length)
    
if save_geom:
    #Save geometry arrays for Autogrid
    stages_rm = [i[0] for i in result[5]]
    stages_H1 = [i[1] for i in result[5]]
    stages_H2 = [i[2] for i in result[5]]
    stages_H3 = [i[3] for i in result[5]]
    stages_Cst = [i[4] for i in result[5]]
    stages_Cro = [i[5] for i in result[5]]
    stages_phi = result[16][0]
    hub_r = [stages_rm[0]-stages_H1[0]/2]
    hub_x = [0]
    case_r = [stages_rm[0]+stages_H1[0]/2]
    case_x = [0]
    stages_Cx = []
    stages_X1_hub = []
    stages_X2_hub = []
    stages_X1_case = []
    stages_X2_case = []
    stages_X1 = []
    stages_X2 = []
    stages_blades = []
    z = 0
    n = n
    for i in range(n):
        
        a1s, a2s, b2s, b3s = free_vortex(result[10][i], [stages_rm[i], stages_H1[i], stages_H2[i]], stages_phi[i])
        hub_r.append(stages_rm[i]-stages_H2[i]/2)
        hub_x.append(z+1.5*stages_Cst[i])
        hub_r.append(stages_rm[i]-stages_H3[i]/2)
        hub_x.append(z+1.5*stages_Cst[i]+1.5*stages_Cro[i])
        case_r.append(stages_rm[i]+stages_H2[i]/2)
        case_x.append(z+1.5*stages_Cst[i])
        case_r.append(stages_rm[i]+stages_H3[i]/2)
        case_x.append(z+1.5*stages_Cst[i]+1.5*stages_Cro[i])
        stages_Cx.append(stages_Cst[i])
        stages_Cx.append(stages_Cro[i])
        stages_X1.append(a1s[0])
        stages_X1.append(a1s[1])
        stages_X1.append(a1s[2])
        stages_X2.append(a2s[0])
        stages_X2.append(a2s[1])
        stages_X2.append(a2s[2])
        stages_X1.append(b2s[0])
        stages_X1.append(b2s[1])
        stages_X1.append(b2s[2])
        stages_X2.append(b3s[0])
        stages_X2.append(b3s[1])
        stages_X2.append(b3s[2])
        stages_blades.append(result[5][i][10])
        stages_blades.append(result[5][i][11])
        
        z += 1.5*(stages_Cst[i]+stages_Cro[i])

    import csv
    
    with open('20_stage.csv', mode='w') as turbine_geom:
        turbine_writer = csv.writer(turbine_geom, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)

        turbine_writer.writerow(stages_blades)
        turbine_writer.writerow(stages_X1)
        turbine_writer.writerow(stages_X2)
        turbine_writer.writerow(stages_Cx)
        turbine_writer.writerow(hub_x)
        turbine_writer.writerow(hub_r)
        turbine_writer.writerow(case_x)
        turbine_writer.writerow(case_r)


if plot == 'yes':
    from GUI import b2b_variable, b2b_plot, annulus
    b2b_variable(result)
    # b2b_plot(result)
    # annulus(result)

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
        eff = turbine(Po1, To1, mdot, Omega, phi, psi, Lambda, AR, W*n/(n+1), dho, n, t, g, ptc, a1i, gas)[0]
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
    
    print("Optimiser time: {}".format(time.time()-start))
    
    turbine_data = turbine(Po1, To1, mdot, Omega, W, t, g, phi, psi, Lambda, AR, dho, n, ptc, ain, gas)
    
    print("Optimum efficiency = {}".format(turbine_data[0]))
    print("phi = {}".format(np.round(phi, 4)))
    print("psi = {}".format(np.round(psi, 4)))
    print("Lambda = {}".format(np.round(Lambda, 4)))
    print("AR = {}".format(np.round(AR, 4)))
    print("dho = {}".format(np.round(dho, 4)))
    print("ain = {}".format(np.round(turbine_data[11][14], 4)))
    print("Max angle = {}".format(np.round(turbine_data[12][1], 4)))
    print('Mass = {} kg'.format(turbine_data[2]))
    print('No. Blades = {}'.format(int(turbine_data[6])))
    print('Axial force on rotor = {} N'.format(turbine_data[13]))
    print('Average Re = {}'.format(turbine_data[14]))
    
    if plot == 'opt':
        from GUI import b2b_variable, b2b_plot, annulus
        # b2b_plot(turbine_data)
        annulus(turbine_data)
        
if calcs == 'optall':
    
    start=time.time()
    #Create arrays of all the loading combinations 
    phi_set = np.arange(0.2, 0.9, 0.3)
    psi_set = np.arange(0.6, 3, 0.4)
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
                    for AR in np.arange(0.9, 3.0, 0.5):
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
             
        eff = turbine(Po1, To1, mdot, Omega, W, t, g, phi, psi, Lambda, AR, dho, n, ptc, ain, gas)[0]
        
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
    
    p = Pool(processes=1)
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
        
if calcs == 'sensitivity':
    
    result = turbine(Po1, To1, mdot, Omega, W, t, g, phi, psi, Lambda, AR, dho, n, ptc, ain, gas)
    base_eff = result[0]

    Po1, To1, mdot, Omega, W, t, g, phi, psi, Lambda, AR, dho, n, ptc, ain = result[11]

    for ph in np.arange(0.9*phi[0], 1.1*phi[0], 0.01*phi[0]):
        phi_0 = [ph, phi[-1]]
        new_eff = turbine(Po1, To1, mdot, Omega, W, t, g, phi_0, psi, Lambda, AR, dho, n, ptc, ain)[0]
        if abs(base_eff-new_eff)/base_eff > 0.001:
            print('phi 0:', ph, new_eff)
        phi_0 = [phi[0], ph]
        new_eff = turbine(Po1, To1, mdot, Omega, W, t, g, phi_0, psi, Lambda, AR, dho, n, ptc, ain)[0]
        if abs(base_eff-new_eff)/base_eff > 0.001:
            print('phi 1:', ph, new_eff)

    for ps in np.arange(0.9*psi[0], 1.1*psi[0], 0.01*psi[0]):
        psi_0 = [ps, psi[-1]]
        new_eff = turbine(Po1, To1, mdot, Omega, W, t, g, phi, psi_0, Lambda, AR, dho, n, ptc, ain)[0]
        if abs(base_eff-new_eff)/base_eff > 0.001:
            print('psi 0:', ps, new_eff)
        psi_0 = [psi[0], ps]
        new_eff = turbine(Po1, To1, mdot, Omega, W, t, g, phi, psi_0, Lambda, AR, dho, n, ptc, ain)[0]
        if abs(base_eff-new_eff)/base_eff > 0.001:
            print('psi 1:', ps, new_eff)

    for la in np.arange(0.9*Lambda[0], 1.1*Lambda[0], 0.01*Lambda[0]):
        Lambda_0 = [la, Lambda[-1]]
        new_eff = turbine(Po1, To1, mdot, Omega, W, t, g, phi, psi, Lambda_0, AR, dho, n, ptc, ain)[0]
        if abs(base_eff-new_eff)/base_eff > 0.001:
            print('Lambda 0:', la, new_eff)
        phi_0 = [Lambda[0], la]
        new_eff = turbine(Po1, To1, mdot, Omega, W, t, g, phi, psi, Lambda_0, AR, dho, n, ptc, ain)[0]
        if abs(base_eff-new_eff)/base_eff > 0.001:
            print('Lambda 1:', la, new_eff)

    for dh in np.arange(0.9*dho[0], 1.1*dho[0], 0.01*dho[0]):
        dho_0 = [dh, dho[-1]]
        new_eff = turbine(Po1, To1, mdot, Omega, W, t, g, phi, psi, Lambda, AR, dho_0, n, ptc, ain)[0]
        if abs(base_eff-new_eff)/base_eff > 0.001:
            print('dho 0:', dh, new_eff)
        dho_0 = [dho[0], dh]
        new_eff = turbine(Po1, To1, mdot, Omega, W, t, g, phi, psi, Lambda, AR, dho_0, n, ptc, ain)[0]
        if abs(base_eff-new_eff)/base_eff > 0.001:
            print('dho 1:', dh, new_eff)

    for ar in np.arange(0.9*AR[0], 1.1*AR[0], 0.01*AR[0]):
        ar_0 = [ar, AR[-1]]
        new_eff = turbine(Po1, To1, mdot, Omega, W, t, g, phi, psi, Lambda, ar_0, dho, n, ptc, ain)[0]
        if abs(base_eff-new_eff)/base_eff > 0.001:
            print('AR 0:', ar, new_eff)
        ar_0 = [AR[0], ar]
        new_eff = turbine(Po1, To1, mdot, Omega, W, t, g, phi, psi, Lambda, ar_0, dho, n, ptc, ain)[0]
        if abs(base_eff-new_eff)/base_eff > 0.001:
            print('AR 1:', ar, new_eff)
