"""Calculations for evaluating turbine performance"""

from turbine import turbine, optimise
import matplotlib.pyplot as plt
import numpy as np
import time
import tikzplotlib

Po1 = 145*10**5
To1 = 950
W = 17*10**6
mdot = 16
Omega = 6782*2*np.pi/60
t = 0.0003
g = 0.0003
cp = 5187
gas = 'He'

phi = [0.272, 0.293]
psi = [0.778, 0.901]
Lambda = [0.496, 0.507]
AR = [0.643, 0.796]
ptc = [1.1, 1.1]
n = 10
ain = 0
dho = [1.0, 1.106]

# phi = [0.28, 0.32]
# psi = [0.83, 1.08]
# Lambda = [0.50, 0.52]
# AR = [0.48, 0.62]
# ptc = [1.1, 1.1]
# n = 5
# ain = 0
# dho = [1.00, 1.11]

# Properties for representative real turbine
# Po1 = 16*10**5
# To1 = 1500
# W = 15*10**6
# mdot = 32.1
# Omega = 1169
# t = 0.001
# g = 0.0002
# gas = 'A1'

# # Optimised parameters
# phi = [0.22, 0.27]
# psi = [0.55, 0.65]
# Lambda = [0.46, 0.51]
# AR = [0.61, 0.74]
# ptc = [1.1, 1.1]
# n = 2
# ain = 0
# dho = [1.0, 1.21]

# # Actual parameters
# phi = [0.5, 0.55]
# psi = [2.2, 1.8]
# Lambda = [0.19, 0.10]
# AR = [1.6, 1.6]
# ptc = [1.1, 1.1]
# n = 2
# ain = 0
dho = [1.22, 1.0]

plot = 'opt'
save = ''
start_time = time.time()
result = turbine(Po1, To1, mdot, Omega, W, t, g, phi, psi, Lambda, AR, dho, n, ptc, ain, gas)
print('Time: {} s'.format(time.time()-start_time))
print('Work = {} W'.format(result[1]))
print('Efficiency = {}'.format(result[0]))
print('Mass = {} kg'.format(result[2]))
print('Volume = {} m^3'.format(result[3]))
print('No. Blades = {}'.format(int(result[6])))
print('Axial force on rotor = {} N'.format(result[13]))
print('Average Re = {}'.format(result[14]))
# print('Cold-stat, cold-rot, warm-rot, hot-rot:', result[15])
To3=result[9]
loss_norm = To3*mdot/W

if plot == 'transient':

    H_bar = [(i[1]+i[2]+i[3])/3 for i in result[5]]
    cold_stat=[]
    cold_rot=[]
    warm=[]
    hot=[]
    for i in range(n):
        cold_stat.append(result[15][i][0]/H_bar[i])
        cold_rot.append(result[15][i][1]/H_bar[i])
        warm.append(result[15][i][2]/H_bar[i])
        hot.append(result[15][i][3]/H_bar[i])
    x = np.arange(1, n+1)  # the label locations
    width = 0.2  # the width of the bars
    fig, ax = plt.subplots()
    ax.bar(x - 3*width/2, cold_stat, width, color='tab:blue', label='Cold-static')
    ax.bar(x - width/2, cold_rot, width, color='tab:green', label='Cold-rotating')    
    ax.bar(x + width/2, warm, width, color='tab:orange', label='Warm transient')
    ax.bar(x + 3*width/2, hot, width, color='tab:red', label='Hot soaked')
    ax.set_ylabel('Change in clearance (mm)', fontsize=20)
    ax.set_xlabel('Stage number', fontsize=20)
    ax.set_xticks(x)
    ax.set_xticklabels(x)
    ax.legend(prop={'size': 18})
    ax.tick_params(axis="x", labelsize=20)
    ax.tick_params(axis="y", labelsize=20)
    tikzplotlib.save("test2.tex")
    fig.tight_layout()
    
if plot == 'mdot':

    te = []
    profile = []
    secondary = []
    tc = []
    ms = np.arange(16,61,1)
    
    for md in ms:
        
        result = turbine(Po1, To1, md, Omega, W, t, g, phi, psi, Lambda, AR, dho, n, ptc, ain, gas)
        To3=result[9]
        loss_norm = To3*mdot/W
        profile.append(loss_norm*result[7][0])
        te.append(loss_norm*result[7][1])
        secondary.append(loss_norm*result[7][2])
        tc.append(loss_norm*result[7][3])
    
    plt.figure()
    plt.plot(ms, profile, label='Profile', linewidth = 3)
    plt.plot(ms, te, label='Trailing edge', linewidth = 3)
    plt.plot(ms, secondary, label='Secondary', linewidth = 3)
    plt.plot(ms, tc, label='Shroud', linewidth = 3)
    plt.xlabel('$\\dot{m}$', fontsize=35)
    plt.ylabel('$\\frac{T_{out}\Delta s}{\Delta h_0}$', fontsize=35)
    plt.tick_params(axis="x", labelsize=30)
    plt.tick_params(axis="y", labelsize=30)
    plt.legend(prop={'size':33})
    tikzplotlib.save("test.tex")
    # plt.yscale('log')
    plt.show()

if plot == 'opt':
    
    Fx = []
    Re = []
    n_blades = []
    phi1 = []
    psi1 = []
    lam1 = []
    ar1 = []
    dh1 = []
    phi2 = []
    psi2 = []
    lam2 = []
    ar2 = []
    phim = []
    psim = []
    lamm = []
    arm = []
    eff = []
    mass = []
    clearance_labels = ['5', '10', '15', '20', '25']
    cold_stat = []
    cold_rot = []
    warm = []
    hot = []
    
    for n in range(5,26):
        
        result = turbine(Po1, To1, mdot, Omega, W, t, g, phi, psi, Lambda, AR, dho, n, ptc, ain, gas)
        
        ph, ps, Lam, ar, dh = optimise(result)
        
        turbine_data = turbine(Po1, To1, mdot, Omega, W, t, g, ph, ps, Lam, ar, dh, n, ptc, ain)
        
        if 1 == 2:
            dims = turbine_data[5]
            angs = turbine_data[10]
        
            d_hubin = [2*i[9] for i in dims]
            d_hubout = [2*(i[0]-i[2]/2) for i in dims]
            d_casout = [2*i[8] for i in dims]
            d_casin = [2*(i[0]+i[2]/2) for i in dims]
            length = [1.5*(i[4]+i[5]) for i in dims]
            a1 = [i[0] for i in angs]
            a2 = [i[1] for i in angs]
            b1 = [i[2] for i in angs]
            b2 = [i[4] for i in angs]
            
            np.savez('Dims {}'.format(n), d_hubin, d_hubout, d_casin, d_casout, length, a1, a2, b1, b2)

        
        n_blades.append(turbine_data[6])
        Re.append(turbine_data[14]/100000)
        Fx.append(turbine_data[13])
        phim.append((ph[0]+ph[-1])/2)
        psim.append((ps[0]+ps[-1])/2)
        lamm.append((Lam[0]+Lam[-1])/2)
        arm.append((ar[0]+ar[-1])/2)
        phi1.append(ph[0])
        psi1.append(ps[0])
        lam1.append(Lam[0])
        ar1.append(ar[0])
        dh1.append(dh[-1]/dh[0])
        phi2.append(ph[-1])
        psi2.append(ps[-1])
        lam2.append(Lam[-1])
        ar2.append(ar[-1])
        eff.append(turbine_data[0])
        mass.append(turbine_data[2])
        
        if np.mod(n,5) == 0:
            c1 = [i[0]*1000 for i in turbine_data[15]]
            c2 = [i[1]*1000 for i in turbine_data[15]]
            w = [i[2]*1000 for i in turbine_data[15]]
            h = [i[3]*1000 for i in turbine_data[15]]
            
            stage_max = np.argmin(w)
            
            cold_stat.append(c1[stage_max])
            cold_rot.append(c2[stage_max])
            warm.append(w[stage_max])
            hot.append(h[stage_max])
        
    
    plt.figure()
    plt.plot(range(5,26), Fx, 'x')
    plt.xlabel('Number of stages', fontsize=20)
    plt.ylabel('Axial force (N)', fontsize=20)
    plt.xticks(range(5,26,5))
    plt.tick_params(axis="x", labelsize=20)
    plt.tick_params(axis="y", labelsize=20)
    
    plt.figure()
    plt.plot(range(5,26), Re, 'x')
    plt.ylabel('Average $Re$ $(\\times 10^6)$', fontsize=20)
    plt.xlabel('Number of stages', fontsize=20)
    plt.xticks(range(5,26,5))
    plt.tick_params(axis="x", labelsize=20)
    plt.tick_params(axis="y", labelsize=20)
    
    plt.figure()
    plt.plot(range(5,26), mass, 'x')
    plt.ylabel('Mass (kg)', fontsize=20)
    plt.xlabel('Number of stages', fontsize=20)
    plt.xticks(range(5,26,5))
    plt.tick_params(axis="x", labelsize=20)
    plt.tick_params(axis="y", labelsize=20)
    
    plt.figure()
    x = np.arange(len(clearance_labels))  # the label locations
    width = 0.2  # the width of the bars
    fig, ax = plt.subplots()
    rects1 = ax.bar(x - 3*width/2, cold_stat, width, label='Cold-static')
    rects2 = ax.bar(x - width/2, cold_rot, width, label='Cold-rotating')    
    rects3 = ax.bar(x + width/2, warm, width, label='Warm transient')
    rects4 = ax.bar(x + 3*width/2, hot, width, label='Hot soaked')
    ax.set_ylabel('Maximum clearance (mm)', fontsize=20)
    ax.set_xlabel('Number of stages', fontsize=20)
    ax.set_xticks(x)
    ax.set_xticklabels(clearance_labels)
    ax.legend(prop={'size': 18})
    ax.tick_params(axis="x", labelsize=20)
    ax.tick_params(axis="y", labelsize=20)
    fig.tight_layout()
    
    plt.figure()
    plt.plot(range(5,26), n_blades, 'x')
    plt.ylabel('Number of blades', fontsize=17)
    plt.xlabel('Number of stages', fontsize=17)
    plt.tick_params(axis="x", labelsize=15)
    plt.tick_params(axis="y", labelsize=15)
    plt.xticks(range(5,26,5))
    plt.tick_params(axis="x", labelsize=20)
    plt.tick_params(axis="y", labelsize=20)
#    plt.savefig("n_blades.eps",format='eps',bbox_inches='tight')
    
    plt.figure()
    plt.plot(range(5,26), phim, label='$\\phi$', linewidth = 2)
    plt.plot(range(5,26), psim, label='$\\psi$', linewidth = 2)
    plt.plot(range(5,26), lamm, label='$\\Lambda$', linewidth = 2)
    plt.plot(range(5,26), arm, label='AR', linewidth = 2)
    plt.plot(range(5,26), dh1, color='gray', label='$\\frac{\\Delta h_{0,2}}{\\Delta h_{0,1}}$', linewidth = 2)
    plt.plot(range(5,26), eff, label='$\\eta$', linewidth = 2)
    plt.ylabel('Optimal parameter value', fontsize=20)
    plt.xlabel('Number of stages', fontsize=20)
    plt.xticks(range(5,26,5))
    plt.tick_params(axis="x", labelsize=20)
    plt.tick_params(axis="y", labelsize=20)
    plt.legend(bbox_to_anchor=(1.0,0.5), loc="center left", prop={'size': 18})
    
    plt.figure()
    plt.plot(range(5,26), phi1, color='blue', label='$\\phi_1$')
    plt.plot(range(5,26), phi2, color='blue', linestyle='--', label='$\\phi_2$')
    plt.plot(range(5,26), psi1, color='orange', label='$\\psi_1$')
    plt.plot(range(5,26), psi2, color='orange', linestyle='--', label='$\\psi_2$')
    plt.plot(range(5,26), lam1, color='green', label='$\\Lambda_1$')
    plt.plot(range(5,26), lam2, color='green', linestyle='--', label='$\\Lambda_2$')
    plt.plot(range(5,26), ar1, color='red', label='$AR_1$')
    plt.plot(range(5,26), ar2, color='red', linestyle='--', label='$AR_2$')
    plt.plot(range(5,26), dh1, color='gray', label='$\\frac{\\Delta h_{0,2}}{\\Delta h_{0,1}}$', linewidth = 2)
    plt.plot(range(5,26), eff, label='$\\eta$', linewidth = 2)
    plt.ylabel('Optimal parameter value', fontsize=17)
    plt.xlabel('Number of stages', fontsize=17)
    plt.legend(bbox_to_anchor=(1.0,0.5), loc="center left", prop={'size':15})
    plt.tick_params(axis="x", labelsize=15)
    plt.tick_params(axis="y", labelsize=15)
    # plt.savefig("Optimal.eps",format='eps',bbox_inches='tight')

    plt.xticks(range(5,26,5))
    plt.show()
        
if plot == 'mechs':
    
    te = []
    profile = []
    secondary = []
    tc = []
    
    for ph in np.arange(0.2,1.0,0.05):
        
        result = turbine(Po1, To1, mdot, Omega, W, t, g, ph, psi, Lambda, AR, dho, n, ptc, ain,gas)
        
        profile.append(loss_norm*result[7][0])
        te.append(loss_norm*result[7][1])
        secondary.append(loss_norm*result[7][2])
        tc.append(loss_norm*result[7][3])
    
    plt.figure()
    plt.plot(np.arange(0.2,1.0,0.05), profile, label='Profile', linewidth = 3)
    plt.plot(np.arange(0.2,1.0,0.05), te, label='Trailing edge', linewidth = 3)
    plt.plot(np.arange(0.2,1.0,0.05), secondary, label='Secondary', linewidth = 3)
    plt.plot(np.arange(0.2,1.0,0.05), tc, label='Shroud', linewidth = 3)
    # plt.xlabel('$\\phi$', fontsize=35)
    # plt.ylabel('$\\frac{T_{out}\Delta s}{\Delta h_0}$', fontsize=35)
    plt.xlabel('Flow coefficient', fontsize=35)
    plt.ylabel('$\\Delta \\eta = \\frac{T_{out}\Delta s}{\Delta h_0}$', fontsize=35)
    plt.tick_params(axis="x", labelsize=30)
    plt.tick_params(axis="y", labelsize=30)
    plt.legend(prop={'size':33})
    plt.show()
    
    te = []
    profile = []
    secondary = []
    tc = []
    
    for ps in np.arange(0.5,2.5,0.05):
        
        result = turbine(Po1, To1, mdot, Omega, W, t, g, phi, ps, Lambda, AR, dho, n, ptc, ain)
        
        profile.append(loss_norm*result[7][0])
        te.append(loss_norm*result[7][1])
        secondary.append(loss_norm*result[7][2])
        tc.append(loss_norm*result[7][3])
        
    plt.figure()
    plt.plot(np.arange(0.5,2.5,0.05), profile, label='Profile', linewidth = 3)
    plt.plot(np.arange(0.5,2.5,0.05), te, label='Trailing edge', linewidth = 3)
    plt.plot(np.arange(0.5,2.5,0.05), secondary, label='Secondary', linewidth = 3)
    plt.plot(np.arange(0.5,2.5,0.05), tc, label='Shroud', linewidth = 3)
    # plt.xlabel('$\\psi$', fontsize=35)
    # plt.ylabel('$\\frac{T_{out}\Delta s}{\Delta h_0}$', fontsize=35)
    plt.xlabel('Stage loading coefficient', fontsize=35)
    plt.ylabel('$\\Delta \\eta = \\frac{T_{out}\Delta s}{\Delta h_0}$', fontsize=35)
    plt.tick_params(axis="x", labelsize=30)
    plt.tick_params(axis="y", labelsize=30)
    plt.show()
    
    te = []
    profile = []
    secondary = []
    tc = []
    
    for ns in np.arange(5,26,1):
        
        result = turbine(Po1, To1, mdot, Omega, W, t, g, phi, psi, Lambda, AR, dho, ns, ptc, ain)
        
        profile.append(loss_norm*result[7][0])
        te.append(loss_norm*result[7][1])
        secondary.append(loss_norm*result[7][2])
        tc.append(loss_norm*result[7][3])
        
    plt.figure()
    plt.xticks(range(5,26,5))
    plt.plot(np.arange(5,26,1), profile, label='Profile', linewidth = 3)
    plt.plot(np.arange(5,26,1), te, label='Trailing edge', linewidth = 3)
    plt.plot(np.arange(5,26,1), secondary, label='Secondary', linewidth = 3)
    plt.plot(np.arange(5,26,1), tc, label='Shroud', linewidth = 3)
    # plt.xlabel('$n$', fontsize=35)
    # plt.ylabel('$\\frac{T_{out}\Delta s}{\Delta h_0}$', fontsize=35)
    plt.xlabel('Number of stages', fontsize=35)
    plt.ylabel('$\\Delta \\eta = \\frac{T_{out}\Delta s}{\Delta h_0}$', fontsize=35)
    plt.tick_params(axis="x", labelsize=30)
    plt.tick_params(axis="y", labelsize=30)
    plt.show()
    
    te = []
    profile = []
    secondary = []
    tc = []
    
    for ar in np.arange(0.4,2,0.05):
        
        result = turbine(Po1, To1, mdot, Omega, W, t, g, phi, psi, Lambda, ar, dho, n, ptc, ain)
        
        profile.append(loss_norm*result[7][0])
        te.append(loss_norm*result[7][1])
        secondary.append(loss_norm*result[7][2])
        tc.append(loss_norm*result[7][3])
    
    plt.figure()
    plt.plot(np.arange(0.4,2,0.05), profile, label='Profile', linewidth = 3)
    plt.plot(np.arange(0.4,2,0.05), te, label='Trailing edge', linewidth = 3)
    plt.plot(np.arange(0.4,2,0.05), secondary, label='Secondary', linewidth = 3)
    plt.plot(np.arange(0.4,2,0.05), tc, label='Shroud', linewidth = 3)
    # plt.xlabel('AR', fontsize=35)
    # plt.ylabel('$\\frac{T_{out}\Delta s}{\Delta h_0}$', fontsize=35)
    plt.xlabel('Blade aspect ratio', fontsize=35)
    plt.ylabel('$\\Delta \\eta = \\frac{T_{out}\Delta s}{\Delta h_0}$', fontsize=35)
    plt.tick_params(axis="x", labelsize=30)
    plt.tick_params(axis="y", labelsize=30)
    plt.show()

    
    
