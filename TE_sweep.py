from Stage import turbine, angles, velocities, Y_po, p_w, losses
import numpy as np

Po1 = 145*10**5
To1 = 950
W = 17*10**6
mdot = 16
Omega = 6782*2*np.pi/60
t = 0.0003
g = 0.0005

#20 stages
"""
phi = 0.3
psi = 0.8
Lambda = 0.5
AR = 1.5
ptoC = 1.1
n = 20

"""
#10 stages
phi = 0.3
psi = 0.8
Lambda = 0.5
AR = 1.0
ptoC = 1.056
n = 10
r = 0.513
H_st = 0.00626
H_ro = 0.0065

"""
#5 stages
phi = 0.4
psi = 1.1
Lambda = 0.5
AR = 1.0
ptoC = 1.0
n = 5
"""

result = turbine(Po1, To1, mdot, Omega, phi, psi, Lambda, AR, W, n, t, g, ptoC)
r = result[4][0][0]
H_st = result[4][0][1]
H_ro = result[4][0][2]
psi = result[10][0]
g=0

#Gas properties
cp = 5187
gamma = 1.6625
R = 2067
#Set the specific enthalpy drop for each stage
del_ho = W/(mdot*n)
#Mean radius and blade speed, constant throughout
U = r*Omega
Vx = U*phi #Axial velocity fixed by turbine parameters
#Calculate angles
a1, a2, b2, a3, b3 = [i for i in angles(phi, psi, Lambda)]
#Create array of angles needed in the loss function
angs = [a1, a2, b2, b3]
#Stage sizing
Cx_st = H_st/AR
Cx_ro = H_ro/AR
w_st = p_w(Cx_st, a1, a2, ptoC)[0]
w_ro = p_w(Cx_ro, b2, b3, ptoC)[0]
#Calculate velocities
V1, V2, W2, V3, W3 = [i for i in velocities([a1, a2, b2, a3, b3], Vx)]
#Create array of velocities needed in the loss function
vels = [V1, V2, W2, W3]
#Temperatures, speed of sound and viscosity through the stage
T1 = To1 - 0.5*(V1**2)/cp
c1 = np.sqrt(gamma*R*T1)
To2 = To1 #No stagnation temperature drop over stator
T2 = To2 - 0.5*(V2**2)/cp
c2 = np.sqrt(gamma*R*T2)
mu2 = (3.674*10**(-7))*(T2**0.7)
To3 = To2 - del_ho/cp #Stagnation temperature drop depending on stage work
T3 = To3 - 0.5*(V3**2)/cp
c3 = np.sqrt(gamma*R*T3)
mu3 = (3.674*10**(-7))*(T3**0.7)

#Flow is incompressible however compressible relation requires no
#knowledge of the density to find the static pressure
P1 = Po1*(1+(gamma-1)*((V1/c1)**2)/2)**(-gamma/(gamma-1))
rho1 = P1/(R*T1)
for t in np.arange(0,0.00105,0.00005):
     #Initialise loss coefficients
    Y_st_guess = 0.01
    Y_ro_guess = 0.01
    Y_st = 0
    Y_ro = 0
    #Iterate until the loss coefficients are stable
    while abs(Y_st_guess-Y_st)/Y_st_guess > 0.001 and abs(Y_ro_guess-Y_ro)/Y_ro_guess > 0.001:
        #Update the loss coefficients
        Y_st = Y_st_guess
        Y_ro = Y_ro_guess
        #Stator to rotor pressure change
        Po2 = Y_po(Po1, Y_st, V2, c2, gamma) #Stagnation pressure after loss
        P2 = Po2*(1+(gamma-1)*((V2/c2)**2)/2)**(-gamma/(gamma-1))
        #Move into relative frame for the rotor
        Po2rel = P2*(1+(gamma-1)*((W2/c2)**2)/2)**(gamma/(gamma-1))
        rho2 = P2/(R*T2)
        #Exit
        Po3rel = Y_po(Po2rel, Y_ro, W3, c3, gamma)
        #Move out of the relative frame to find exit states
        P3 = Po3rel*(1+(gamma-1)*((W3/c3)**2)/2)**(-gamma/(gamma-1))
        Po3 = P3*(1+(gamma-1)*((V3/c3)**2)/2)**(gamma/(gamma-1))
        rho3 = P3/(R*T3)
        #Create array for the loss function
        dimensions = [t, g, r, H_st, H_ro, Cx_st, Cx_ro, w_st, w_ro]    
        #Create array for the loss function
        states = [T1, T2, T3, rho1, rho2, rho3, mu2, mu3, mdot]
        #Find the losses across the stator and rotor
        loss_calc = losses(angs, vels, states, dimensions)
        loss_st = loss_calc[1] #Stator entropy rise
        loss_ro = loss_calc[2] #Rotor entropy rise
        #Use entropy rise to find stagnation pressure loss coefficients
        #assuming entropy loss coefficient is equivalent at low Mach numbers
        Y_st_guess = loss_st*T2/(0.5*V2**2)
        Y_ro_guess = loss_ro*T3/(0.5*W3**2)
    
    loss = loss_calc[0]
    eff = del_ho/(del_ho+To3*loss)

    print(eff)
    

