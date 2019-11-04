"""Calculations for evaluating turbine performance"""

from Stage import turbine, free_vortex
import numpy as np
import matplotlib.pyplot as plt
import time

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

#20 stages
"""
phi = 0.3
psi = 0.8
Lambda = 0.5
AR = 1.5
ptoC = -1.0
n = 20
"""
"""
#10 stages
phi = 0.3
psi = 0.8
Lambda = 0.5
AR = 1.0
ptoC = 1.056
n = 10
"""
"""
#5 stages
phi = 0.4
psi = 1.1
Lambda = 0.5
AR = 1.0
ptoC = 1.0
n = 5
"""

#10 stages
phi = [0.3, 0.3, 0.3]
psi = [0.8, 0.8, 0.8]
Lambda = [0.5, 0.5, 0.5]
AR = 1.0
ptoC = 1.056
n = 5
dho = [1, 1 , 1]

phimin = 0.1
phimax = 0.5
psimin = 0.7
psimax = 3
Lambdamin = 0.1
Lambdamax = 0.9
ARmin = 0.5
ARmax = 3
nmin = 10
nmax = 20
size = len(np.arange(ARmin, ARmax+0.2, 0.2))

plots = ''
save = ''
start_time = time.time()
result = turbine(Po1, To1, mdot, Omega, phi, psi, Lambda, AR, W, dho, n, t, g, ptoC)
#print('Angles [a1,a2,b2,a3,b3]=', result[9])
#print('Sizes [r,Hst,Hro,Cxst,Cxro,ptcst,ptcro]=', result[5][0])
print('Work =', result[1])
print('Efficiency =', result[0])
print('Mass =', result[2])
print('Volume =', result[3])
print('No. Blades =', result[6])
#print('Free vortex angles [a1,a2,b2,b3][h,m,t]=', free_vortex(result[9][0], result[4][0][:3], phi))
print('Time:', time.time()-start_time)

if plots == 'pareto':
    start_time = time.time()
    matrix = []

    for n in np.arange(5, 11):
        print(n)
        for phi in np.arange(0.2, 1.0, 0.3):
            for psi in np.arange(0.5, 2.5, 0.5):
                for Lambda in [0.5]:
                    if abs((1+0.5*psi-Lambda)/phi)>np.tan(1.2741) or abs((0.5*psi+Lambda)/phi)>np.tan(1.2741):
                        continue
                    for AR in np.arange(1.0, 2.0, 0.5):
                        result = turbine(Po1, To1, mdot, Omega, phi, psi, Lambda, AR, W, n, t, g, ptoC)
                        eff = result[0]
                        vol = result[3]
                        matrix.append([eff, vol, n, phi, psi, Lambda, AR])

    not_pareto = []
    for i in range(len(matrix)):
        for j in range(len(matrix)):
            if matrix[i][0] < matrix[j][0] and matrix[i][1] > matrix[j][1] and matrix[i][2] >= matrix[j][2]:
                not_pareto.append(i)
                break

    pareto = []
    for i in range(len(matrix)):
        if i in not_pareto:
            continue
        else:
            pareto.append(matrix[i])

    pareto.sort(key=lambda x: (x[2], x[0]))
    volumes = []
    efficiencies = []
    stages = []
    for i in pareto:
        volumes.append(i[1])
        efficiencies.append(i[0])
        stages.append(i[2])
    a = np.asarray([volumes, efficiencies, stages])
    np.savetxt("pareto.csv",a,delimiter=",")
    labels = np.arange(10, 21)
    for l in zip(labels):
        start = next((i for i, x in enumerate(stages) if x == l), None)
        end = next((i for i, x in enumerate(stages) if x > l), None)
        plt.scatter(volumes[start:end], efficiencies[start:end], s=3, label=l)
        plt.legend(labels, loc='upper left', bbox_to_anchor=(1, 1), markerscale=2)
    plt.xlabel('Volume')
    plt.ylabel('Efficiency')
    if save == 'yes':
        plt.savefig('pareto.png', dpi=1000, bbox_inches='tight')
    plt.show()
    print(time.time()-start_time)

if plots == 'te':
    x = np.arange(0.000,0.00205,0.00005)
    size = len(x)
    TEmat = np.zeros((nmax-nmin+1, size))
    eff = np.zeros((nmax-nmin+1, size))
    for n in range(nmin, nmax+1):
        row_i = n-nmin
        col_i = 0
        for t in x:
            results = turbine(Po1, To1, mdot, Omega, phi, psi, Lambda, AR, W, n, t, g, ptoC)

            eff[row_i][col_i] = results[0]
            TEmat[row_i][col_i] = results[1][1]*To3/del_ho
            
            col_i += 1

    labels = np.arange(nmin, nmax+1)

    for i, l in zip(eff, labels):
        plt.plot(x*1000, i, label=l)
        plt.legend(labels, loc='upper left', bbox_to_anchor=(1, 1))
    plt.xlabel('Trailing edge thickness (mm)')
    plt.ylabel('Efficiency')
    plt.title('Isentropic efficiency')
    plt.show()

    for i, l in zip(TEmat, labels):
        plt.plot(x*1000, i, label=l,zorder=l)
        plt.legend(labels, loc='upper left', bbox_to_anchor=(1, 1))
    plt.xlabel('Trailing edge thickness (mm)')
    plt.ylabel('Normalised loss')
    plt.grid(linestyle='--',zorder=0)
    if save == 'yes':
        plt.savefig('te.png', dpi=1000, bbox_inches='tight')
    plt.show()
    
if plots == 'ptc':
    x = np.arange(0.5,2.5,0.05)
    size = len(x)
    PCmat = np.zeros((nmax-nmin+1, size))
    eff = np.zeros((nmax-nmin+1, size))
    for n in range(nmin, nmax+1):
        row_i = n-nmin
        col_i = 0
        for ptoC in x:
            results = turbine(Po1, To1, mdot, Omega, phi, psi, Lambda, AR, W, n, t, g, ptoC)

            eff[row_i][col_i] = results[0]
            PCmat[row_i][col_i] = results[1][1]*To3/del_ho
            
            col_i += 1

    labels = np.arange(nmin, nmax+1)

    for i, l in zip(eff, labels):
        plt.plot(x, i, label=l)
        plt.legend(labels, loc='upper left', bbox_to_anchor=(1, 1))
    plt.xlabel('Pitch to chord ratio')
    plt.ylabel('Efficiency')
    plt.title('Isentropic efficiency')
    plt.show()

    for i, l in zip(PCmat, labels):
        plt.plot(x, i, label=l,zorder=l)
        plt.legend(labels, loc='upper left', bbox_to_anchor=(1, 1))
    plt.xlabel('Pitch to chord ratio')
    plt.ylabel('Normalised loss')
    plt.grid(linestyle='--',zorder=0)
    if save == 'yes':
        plt.savefig('te.png', dpi=1000, bbox_inches='tight')
    plt.show()


if plots == 'tc':
    x = np.arange(0.0,0.0015,0.00005)
    size = len(x)
    TCmat = np.zeros((nmax-nmin+1, size))
    eff = np.zeros((nmax-nmin+1, size))
    for n in range(nmin, nmax+1):
        row_i = n-nmin
        col_i = 0
        for g in x:
            results = turbine(Po1, To1, mdot, Omega, phi, psi, Lambda, AR, W, n, t, g, ptoC)

            eff[row_i][col_i] = results[0]
            TCmat[row_i][col_i] = results[1][4]*To3/del_ho
            
            col_i += 1

    labels = np.arange(nmin, nmax+1)

    for i, l in zip(eff, labels):
        plt.plot(x*1000, i, label=l)
        plt.legend(labels, loc='upper left', bbox_to_anchor=(1, 1))
    plt.xlabel('Tip clearance (mm)')
    plt.ylabel('Efficiency')
    plt.title('Isentropic efficiency')
    plt.show()

    for i, l in zip(TCmat, labels):
        plt.plot(x*1000, i, label=l)
        plt.legend(labels, loc='upper left', bbox_to_anchor=(1, 1))
    plt.xlabel('Tip clearance (mm)')
    plt.ylabel('Normalised loss',fontsize=13)
    if save == 'yes':
        plt.savefig('tc.png', dpi=1000, bbox_inches='tight')
    plt.show()

            
if plots == 'n':
    size = len(np.arange(ARmin, ARmax+0.2, 0.2))

    te = []
    profile = []
    endwall = []
    secondary = []
    tc = []

    eff = np.zeros((nmax-nmin+1, size))
    span = np.zeros((nmax-nmin+1, nmax))
    chord_st = np.zeros((nmax-nmin+1, size))
    radius = np.zeros((nmax-nmin+1, size))
    length = np.zeros((nmax-nmin+1, size))
    volume = np.zeros((nmax-nmin+1, size))
    TEmat = np.zeros((nmax-nmin+1, size))
    profilemat = np.zeros((nmax-nmin+1, size))
    endwallmat = np.zeros((nmax-nmin+1, size))
    secondarymat = np.zeros((nmax-nmin+1, size))
    TCmat = np.zeros((nmax-nmin+1, size))

    for n in range(nmin, nmax+1):

        results = turbine(Po1, To1, mdot, Omega, phi, psi, Lambda, AR, W, n, t, g, ptoC)
        row_i = n-nmin
        col_i = 0

        profile.append(results[1][0]*To3/del_ho)
        te.append(results[1][1]*To3/del_ho)
        endwall.append(results[1][2]*To3/del_ho)
        secondary.append(results[1][3]*To3/del_ho)
        tc.append(results[1][4]*To3/del_ho)

        for i in range(n):
            span[row_i][i] = results[4][i][1]

        for ar in np.arange(ARmin, ARmax+0.2, 0.2):

            result = turbine(Po1, To1, mdot, Omega, phi, psi, Lambda, ar, W, n, t, g, ptoC)

            eff[row_i][col_i] = result[0]
            chord_st[row_i][col_i] = result[4][0][2]
            radius[row_i][col_i] = result[4][0][0]
            profilemat[row_i][col_i] = result[1][0]*To3/del_ho
            TEmat[row_i][col_i] = result[1][1]*To3/del_ho
            endwallmat[row_i][col_i] = result[1][2]*To3/del_ho
            secondarymat[row_i][col_i] = result[1][3]*To3/del_ho
            TCmat[row_i][col_i] = result[1][4]*To3/del_ho
            length[row_i][col_i] = result[2]
            volume[row_i][col_i] = result[3]

            col_i += 1

    plt.scatter(np.arange(nmin, nmax+1), te, label="TE")
    plt.scatter(np.arange(nmin, nmax+1), tc, label="TC")
    plt.scatter(np.arange(nmin, nmax+1), profile, label="Profile")
    plt.scatter(np.arange(nmin, nmax+1), endwall, label="Endwall")
    plt.scatter(np.arange(nmin, nmax+1), secondary, label="Secondary")
    plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
    plt.xlabel('Number of stages')
    plt.ylabel('Normalised loss')
    if save == 'yes':
        plt.savefig('components.png', dpi=1000, bbox_inches='tight')
    plt.show()

    labels = np.arange(nmin, nmax+1)
    x = np.arange(ARmin, ARmax+0.2, 0.2)

    for i, l in zip(eff, labels):
        plt.plot(x, i, label=l)
        plt.legend(labels, loc='upper left', bbox_to_anchor=(1, 1))
    plt.xlabel('Aspect ratio')
    plt.ylabel('Efficiency')
    plt.title('Isentropic efficiency')
    plt.show()

    for i, l in zip(volume, labels):
        plt.plot(x, i, label=l)
        plt.legend(labels, loc='upper left', bbox_to_anchor=(1, 1))
    plt.xlabel('Aspect ratio')
    plt.ylabel('Volume')
    plt.title('Volume')
    plt.show()

    for i, l in zip(radius, labels):
        plt.plot(x, i, label=l)
        plt.legend(labels, loc='upper left', bbox_to_anchor=(1, 1))
    plt.xlabel('Aspect ratio')
    plt.ylabel('Radius')
    plt.title('Radius')
    plt.show()

    for i, l in zip(length, labels):
        plt.plot(x, i, label=l)
        plt.legend(labels, loc='upper left', bbox_to_anchor=(1, 1))
    plt.xlabel('Aspect ratio')
    plt.ylabel('Length')
    plt.title('Length')
    plt.show()

    plt.plot(np.arange(1, nmax+1), span[15][:], label=20)
    plt.plot(np.arange(1, nmax), span[14][:19], label=19)
    plt.plot(np.arange(1, nmax-1), span[13][:18], label=18)
    plt.plot(np.arange(1, nmax-2), span[12][:17], label=17)
    plt.plot(np.arange(1, nmax-3), span[11][:16], label=16)
    plt.plot(np.arange(1, nmax-4), span[10][:15], label=15)
    plt.plot(np.arange(1, nmax-5), span[9][:14], label=14)
    plt.plot(np.arange(1, nmax-6), span[8][:13], label=13)
    plt.plot(np.arange(1, nmax-7), span[7][:12], label=12)
    plt.plot(np.arange(1, nmax-8), span[6][:11], label=11)
    plt.plot(np.arange(1, nmax-9), span[5][:10], label=10)
    plt.plot(np.arange(1, nmax-10), span[4][:9], label=9)
    plt.plot(np.arange(1, nmax-11), span[3][:8], label=8)
    plt.plot(np.arange(1, nmax-12), span[2][:7], label=7)
    plt.plot(np.arange(1, nmax-13), span[1][:6], label=6)
    plt.plot(np.arange(1, nmax-14), span[0][:5], label=5)
    plt.xlabel('Stage')
    plt.ylabel('Span')
    plt.title('Span variation')
    plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
    plt.show()

    for i, l in zip(chord_st, labels):
        plt.plot(x, i, label=l)
        plt.legend(labels, loc='upper left', bbox_to_anchor=(1, 1))
    plt.xlabel('Aspect ratio')
    plt.ylabel('Inlet chord')
    plt.title('Inlet chord')
    plt.show()

    for i, l in zip(profilemat, labels):
        plt.plot(x, i, label=l)
        plt.legend(labels, loc='upper left', bbox_to_anchor=(1, 1))
    plt.xlabel('Aspect ratio')
    plt.ylabel('Profile loss')
    plt.title('Profile loss')
    plt.show()

    for i, l in zip(TEmat, labels):
        plt.plot(x, i, label=l)
        plt.legend(labels, loc='upper left', bbox_to_anchor=(1, 1))
    plt.xlabel('Aspect ratio')
    plt.ylabel('Trailing edge loss')
    plt.title('Trailing edge loss')
    plt.show()

    for i, l in zip(endwallmat, labels):
        plt.plot(x, i, label=l)
        plt.legend(labels, loc='upper left', bbox_to_anchor=(1, 1))
    plt.xlabel('Aspect ratio')
    plt.ylabel('Endwall loss')
    plt.title('Endwall loss')
    plt.show()

    for i, l in zip(secondarymat, labels):
        plt.plot(x, i, label=l)
        plt.legend(labels, loc='upper left', bbox_to_anchor=(1, 1))
    plt.xlabel('Aspect ratio')
    plt.ylabel('Secondary loss')
    plt.title('Secondary loss')
    plt.show()

    for i, l in zip(TCmat, labels):
        plt.plot(x, i, label=l)
        plt.legend(labels, loc='upper left', bbox_to_anchor=(1, 1))
    plt.xlabel('Aspect ratio')
    plt.ylabel('Tip clearance loss')
    plt.title('Tip clearance loss')
    plt.show()

if plots == 'phi':
    cols = len(np.arange(ARmin, ARmax+0.2, 0.2))
    rows = len(np.arange(phimin, phimax+0.02, 0.03))

    TE = []
    profile = []
    endwall = []
    secondary = []
    TC = []

    eff = np.zeros((rows, size))
    chord_st = np.zeros((rows, size))
    radius = np.zeros((rows, size))
    length = np.zeros((rows, size))
    volume = np.zeros((rows, size))
    TEmat = np.zeros((rows, size))
    profilemat = np.zeros((rows, size))
    endwallmat = np.zeros((rows, size))
    secondarymat = np.zeros((rows, size))
    TCmat = np.zeros((rows, size))

    row_i = 0
    for phi in np.arange(phimin, phimax+0.02, 0.03):

        results = turbine(Po1, To1, mdot, Omega, phi, psi, Lambda, AR, W, n, t, g, ptoC)

        col_i = 0

        profile.append(results[1][0]*To3/del_ho)
        TE.append(results[1][1]*To3/del_ho)
        endwall.append(results[1][2]*To3/del_ho)
        secondary.append(results[1][3]*To3/del_ho)
        TC.append(results[1][4]*To3/del_ho)

        for ar in np.arange(ARmin, ARmax+0.2, 0.2):

            result = turbine(Po1, To1, mdot, Omega, phi, psi, Lambda, AR, W, n, t, g, ptoC)

            eff[row_i][col_i] = result[0]
            chord_st[row_i][col_i] = result[4][0][2]
            radius[row_i][col_i] = result[4][0][0]
            profilemat[row_i][col_i] = result[1][0]*To3/del_ho
            TEmat[row_i][col_i] = result[1][1]*To3/del_ho
            endwallmat[row_i][col_i] = result[1][2]*To3/del_ho
            secondarymat[row_i][col_i] = result[1][3]*To3/del_ho
            TCmat[row_i][col_i] = result[1][4]*To3/del_ho
            length[row_i][col_i] = result[2]
            volume[row_i][col_i] = result[3]

            col_i += 1

        row_i += 1

    plt.scatter(np.arange(phimin, phimax+0.02, 0.03), TE, label="TE")
    plt.scatter(np.arange(phimin, phimax+0.02, 0.03), TC, label="TC")
    plt.scatter(np.arange(phimin, phimax+0.02, 0.03), profile, label="Profile")
    plt.scatter(np.arange(phimin, phimax+0.02, 0.03), endwall, label="Endwall")
    plt.scatter(np.arange(phimin, phimax+0.02, 0.03), secondary, label="Secondary")
    plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
    plt.xlabel('Flow coefficient')
    plt.ylabel('Efficiency penalty')
    plt.title('Loss breakdown')
    plt.show()

    labels = np.arange(phimin, phimax+0.02, 0.03)
    x = np.arange(ARmin, ARmax+0.2, 0.2)

    for i, l in zip(eff, labels):
        plt.plot(x, i, label=l)
        plt.legend(labels, loc='upper left', bbox_to_anchor=(1, 1))
    plt.xlabel('Aspect ratio')
    plt.ylabel('Efficiency')
    plt.title('Isentropic efficiency')
    plt.show()

    for i, l in zip(volume, labels):
        plt.plot(x, i, label=l)
        plt.legend(labels, loc='upper left', bbox_to_anchor=(1, 1))
    plt.xlabel('Aspect ratio')
    plt.ylabel('Volume')
    plt.title('Volume')
    plt.show()

    for i, l in zip(radius, labels):
        plt.plot(x, i, label=l)
        plt.legend(labels, loc='upper left', bbox_to_anchor=(1, 1))
    plt.xlabel('Aspect ratio')
    plt.ylabel('Radius')
    plt.title('Radius')
    plt.show()

    for i, l in zip(length, labels):
        plt.plot(x, i, label=l)
        plt.legend(labels, loc='upper left', bbox_to_anchor=(1, 1))
    plt.xlabel('Aspect ratio')
    plt.ylabel('Length')
    plt.title('Length')
    plt.show()

    for i, l in zip(chord_st, labels):
        plt.plot(x, i, label=l)
        plt.legend(labels, loc='upper left', bbox_to_anchor=(1, 1))
    plt.xlabel('Aspect ratio')
    plt.ylabel('Inlet chord')
    plt.title('Inlet chord')
    plt.show()

    for i, l in zip(profilemat, labels):
        plt.plot(x, i, label=l)
        plt.legend(labels, loc='upper left', bbox_to_anchor=(1, 1))
    plt.xlabel('Aspect ratio')
    plt.ylabel('Profile loss')
    plt.title('Profile loss')
    plt.show()

    for i, l in zip(TEmat, labels):
        plt.plot(x, i, label=l)
        plt.legend(labels, loc='upper left', bbox_to_anchor=(1, 1))
    plt.xlabel('Aspect ratio')
    plt.ylabel('Trailing edge loss')
    plt.title('Trailing edge loss')
    plt.show()

    for i, l in zip(endwallmat, labels):
        plt.plot(x, i, label=l)
        plt.legend(labels, loc='upper left', bbox_to_anchor=(1, 1))
    plt.xlabel('Aspect ratio')
    plt.ylabel('Endwall loss')
    plt.title('Endwall loss')
    plt.show()

    for i, l in zip(secondarymat, labels):
        plt.plot(x, i, label=l)
        plt.legend(labels, loc='upper left', bbox_to_anchor=(1, 1))
    plt.xlabel('Aspect ratio')
    plt.ylabel('Secondary loss')
    plt.title('Secondary loss')
    plt.show()

    for i, l in zip(TCmat, labels):
        plt.plot(x, i, label=l)
        plt.legend(labels, loc='upper left', bbox_to_anchor=(1, 1))
    plt.xlabel('Aspect ratio')
    plt.ylabel('Tip clearance loss')
    plt.title('Tip clearance loss')
    plt.show()

if plots == 'psi':
    cols = len(np.arange(ARmin, ARmax+0.2, 0.2))
    rows = len(np.arange(psimin, psimax, 0.2))

    TE = []
    profile = []
    endwall = []
    secondary = []
    TC = []

    eff = np.zeros((rows, size))
    chord_st = np.zeros((rows, size))
    radius = np.zeros((rows, size))
    length = np.zeros((rows, size))
    volume = np.zeros((rows, size))
    TEmat = np.zeros((rows, size))
    profilemat = np.zeros((rows, size))
    endwallmat = np.zeros((rows, size))
    secondarymat = np.zeros((rows, size))
    TCmat = np.zeros((rows, size))

    row_i = 0
    for psi in np.arange(psimin, psimax, 0.2):

        results = turbine(Po1, To1, mdot, Omega, phi, psi, Lambda, AR, W, n, t, g, ptoC)

        col_i = 0

        profile.append(results[1][0]*To3/del_ho)
        TE.append(results[1][1]*To3/del_ho)
        endwall.append(results[1][2]*To3/del_ho)
        secondary.append(results[1][3]*To3/del_ho)
        TC.append(results[1][4]*To3/del_ho)

        for ar in np.arange(ARmin, ARmax+0.2, 0.2):

            result = turbine(Po1, To1, mdot, Omega, phi, psi, Lambda, AR, W, n, t, g, ptoC)

            eff[row_i][col_i] = result[0]
            chord_st[row_i][col_i] = result[4][0][2]
            radius[row_i][col_i] = result[4][0][0]
            profilemat[row_i][col_i] = result[1][0]*To3/del_ho
            TEmat[row_i][col_i] = result[1][1]*To3/del_ho
            endwallmat[row_i][col_i] = result[1][2]*To3/del_ho
            secondarymat[row_i][col_i] = result[1][3]*To3/del_ho
            TCmat[row_i][col_i] = result[1][4]*To3/del_ho
            length[row_i][col_i] = result[2]
            volume[row_i][col_i] = result[3]

            col_i += 1

        row_i += 1

    plt.scatter(np.arange(psimin, psimax, 0.2), TE, label="TE")
    plt.scatter(np.arange(psimin, psimax, 0.2), TC, label="TC")
    plt.scatter(np.arange(psimin, psimax, 0.2), profile, label="Profile")
    plt.scatter(np.arange(psimin, psimax, 0.2), endwall, label="Endwall")
    plt.scatter(np.arange(psimin, psimax, 0.2), secondary, label="Secondary")
    plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
    plt.xlabel('Loading coefficient')
    plt.ylabel('Efficiency penalty')
    plt.title('Loss breakdown')
    plt.show()

    labels = np.arange(psimin, psimax, 0.2)
    x = np.arange(ARmin, ARmax+0.2, 0.2)

    for i, l in zip(eff, labels):
        plt.plot(x, i, label=l)
        plt.legend(labels, loc='upper left', bbox_to_anchor=(1, 1))
    plt.xlabel('Aspect ratio')
    plt.ylabel('Efficiency')
    plt.title('Isentropic efficiency')
    plt.show()

    for i, l in zip(volume, labels):
        plt.plot(x, i, label=l)
        plt.legend(labels, loc='upper left', bbox_to_anchor=(1, 1))
    plt.xlabel('Aspect ratio')
    plt.ylabel('Volume')
    plt.title('Volume')
    plt.show()

    for i, l in zip(radius, labels):
        plt.plot(x, i, label=l)
        plt.legend(labels, loc='upper left', bbox_to_anchor=(1, 1))
    plt.xlabel('Aspect ratio')
    plt.ylabel('Radius')
    plt.title('Radius')
    plt.show()

    for i, l in zip(length, labels):
        plt.plot(x, i, label=l)
        plt.legend(labels, loc='upper left', bbox_to_anchor=(1, 1))
    plt.xlabel('Aspect ratio')
    plt.ylabel('Length')
    plt.title('Length')
    plt.show()

    for i, l in zip(chord_st, labels):
        plt.plot(x, i, label=l)
        plt.legend(labels, loc='upper left', bbox_to_anchor=(1, 1))
    plt.xlabel('Aspect ratio')
    plt.ylabel('Inlet chord')
    plt.title('Inlet chord')
    plt.show()

    for i, l in zip(profilemat, labels):
        plt.plot(x, i, label=l)
        plt.legend(labels, loc='upper left', bbox_to_anchor=(1, 1))
    plt.xlabel('Aspect ratio')
    plt.ylabel('Profile loss')
    plt.title('Profile loss')
    plt.show()

    for i, l in zip(TEmat, labels):
        plt.plot(x, i, label=l)
        plt.legend(labels, loc='upper left', bbox_to_anchor=(1, 1))
    plt.xlabel('Aspect ratio')
    plt.ylabel('Trailing edge loss')
    plt.title('Trailing edge loss')
    plt.show()

    for i, l in zip(endwallmat, labels):
        plt.plot(x, i, label=l)
        plt.legend(labels, loc='upper left', bbox_to_anchor=(1, 1))
    plt.xlabel('Aspect ratio')
    plt.ylabel('Endwall loss')
    plt.title('Endwall loss')
    plt.show()

    for i, l in zip(secondarymat, labels):
        plt.plot(x, i, label=l)
        plt.legend(labels, loc='upper left', bbox_to_anchor=(1, 1))
    plt.xlabel('Aspect ratio')
    plt.ylabel('Secondary loss')
    plt.title('Secondary loss')
    plt.show()

    for i, l in zip(TCmat, labels):
        plt.plot(x, i, label=l)
        plt.legend(labels, loc='upper left', bbox_to_anchor=(1, 1))
    plt.xlabel('Aspect ratio')
    plt.ylabel('Tip clearance loss')
    plt.title('Tip clearance loss')
    plt.show()

if plots == 'Lambda':
    cols = len(np.arange(ARmin, ARmax+0.2, 0.2))
    rows = len(np.arange(Lambdamin, Lambdamax+0.1, 0.1))

    TE = []
    profile = []
    endwall = []
    secondary = []
    TC = []

    eff = np.zeros((rows, size))
    chord_st = np.zeros((rows, size))
    radius = np.zeros((rows, size))
    length = np.zeros((rows, size))
    volume = np.zeros((rows, size))
    TEmat = np.zeros((rows, size))
    profilemat = np.zeros((rows, size))
    endwallmat = np.zeros((rows, size))
    secondarymat = np.zeros((rows, size))
    TCmat = np.zeros((rows, size))

    row_i = 0
    for Lambda in np.arange(Lambdamin, Lambdamax+0.1, 0.1):

        results = turbine(Po1, To1, mdot, Omega, phi, psi, Lambda, AR, W, n, t, g, ptoC)

        col_i = 0

        profile.append(results[1][0]*To3/del_ho)
        TE.append(results[1][1]*To3/del_ho)
        endwall.append(results[1][2]*To3/del_ho)
        secondary.append(results[1][3]*To3/del_ho)
        TC.append(results[1][4]*To3/del_ho)

        for ar in np.arange(ARmin, ARmax+0.2, 0.2):

            result = turbine(Po1, To1, mdot, Omega, phi, psi, Lambda, AR, W, n, t, g, ptoC)

            eff[row_i][col_i] = result[0]
            chord_st[row_i][col_i] = result[4][0][2]
            radius[row_i][col_i] = result[4][0][0]
            profilemat[row_i][col_i] = result[1][0]*To3/del_ho
            TEmat[row_i][col_i] = result[1][1]*To3/del_ho
            endwallmat[row_i][col_i] = result[1][2]*To3/del_ho
            secondarymat[row_i][col_i] = result[1][3]*To3/del_ho
            TCmat[row_i][col_i] = result[1][4]*To3/del_ho
            length[row_i][col_i] = result[2]
            volume[row_i][col_i] = result[3]

            col_i += 1

        row_i += 1

    plt.scatter(np.arange(Lambdamin, Lambdamax+0.1, 0.1), TE, label="TE")
    plt.scatter(np.arange(Lambdamin, Lambdamax+0.1, 0.1), TC, label="TC")
    plt.scatter(np.arange(Lambdamin, Lambdamax+0.1, 0.1), profile, label="Profile")
    plt.scatter(np.arange(Lambdamin, Lambdamax+0.1, 0.1), endwall, label="Endwall")
    plt.scatter(np.arange(Lambdamin, Lambdamax+0.1, 0.1), secondary, label="Secondary")
    plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
    plt.xlabel('Reaction')
    plt.ylabel('Efficiency penalty')
    plt.title('Loss breakdown')
    plt.show()

    labels = np.arange(Lambdamin, Lambdamax+0.1, 0.1)
    x = np.arange(ARmin, ARmax+0.2, 0.2)

    for i, l in zip(eff, labels):
        plt.plot(x, i, label=l)
        plt.legend(labels, loc='upper left', bbox_to_anchor=(1, 1))
    plt.xlabel('Aspect ratio')
    plt.ylabel('Efficiency')
    plt.title('Isentropic efficiency')
    plt.show()

    for i, l in zip(volume, labels):
        plt.plot(x, i, label=l)
        plt.legend(labels, loc='upper left', bbox_to_anchor=(1, 1))
    plt.xlabel('Aspect ratio')
    plt.ylabel('Volume')
    plt.title('Volume')
    plt.show()

    for i, l in zip(radius, labels):
        plt.plot(x, i, label=l)
        plt.legend(labels, loc='upper left', bbox_to_anchor=(1, 1))
    plt.xlabel('Aspect ratio')
    plt.ylabel('Radius')
    plt.title('Radius')
    plt.show()

    for i, l in zip(length, labels):
        plt.plot(x, i, label=l)
        plt.legend(labels, loc='upper left', bbox_to_anchor=(1, 1))
    plt.xlabel('Aspect ratio')
    plt.ylabel('Length')
    plt.title('Length')
    plt.show()

    for i, l in zip(chord_st, labels):
        plt.plot(x, i, label=l)
        plt.legend(labels, loc='upper left', bbox_to_anchor=(1, 1))
    plt.xlabel('Aspect ratio')
    plt.ylabel('Inlet chord')
    plt.title('Inlet chord')
    plt.show()

    for i, l in zip(profilemat, labels):
        plt.plot(x, i, label=l)
        plt.legend(labels, loc='upper left', bbox_to_anchor=(1, 1))
    plt.xlabel('Aspect ratio')
    plt.ylabel('Profile loss')
    plt.title('Profile loss')
    plt.show()

    for i, l in zip(TEmat, labels):
        plt.plot(x, i, label=l)
        plt.legend(labels, loc='upper left', bbox_to_anchor=(1, 1))
    plt.xlabel('Aspect ratio')
    plt.ylabel('Trailing edge loss')
    plt.title('Trailing edge loss')
    plt.show()

    for i, l in zip(endwallmat, labels):
        plt.plot(x, i, label=l)
        plt.legend(labels, loc='upper left', bbox_to_anchor=(1, 1))
    plt.xlabel('Aspect ratio')
    plt.ylabel('Endwall loss')
    plt.title('Endwall loss')
    plt.show()

    for i, l in zip(secondarymat, labels):
        plt.plot(x, i, label=l)
        plt.legend(labels, loc='upper left', bbox_to_anchor=(1, 1))
    plt.xlabel('Aspect ratio')
    plt.ylabel('Secondary loss')
    plt.title('Secondary loss')
    plt.show()

    for i, l in zip(TCmat, labels):
        plt.plot(x, i, label=l)
        plt.legend(labels, loc='upper left', bbox_to_anchor=(1, 1))
    plt.xlabel('Aspect ratio')
    plt.ylabel('Tip clearance loss')
    plt.title('Tip clearance loss')
    plt.show()
