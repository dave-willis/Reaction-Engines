import numpy as np
import matplotlib.pyplot as plt
import Stage as f

#Parameters
Poin = 145*10**5 #Pa
Toin = 950 #K
Wout = 17*10**6 #W
mdot = 16 #kg/s
Omega = 6782*2*np.pi/60 #rad/s
Lambda = 0.5

cp = 5190 #J/K
gamma = 1.667
R = 2080 #J/K
rhoin = 7.6 #kg/m^3
mu = 31*10**-6 #Ns/m^2

t = 0.001
g = 0.0005
CD = 0.002

phi = 0.45
psi = 2.4
AR = 1.4

a1 = np.arctan((1-0.5*psi-Lambda)/phi)
b3 = np.arctan(np.tan(a1)-1/phi)
a2 = -1*b3
b2 = np.arctan(np.tan(a2)-1/phi)
a3 = a1

delta_ho = Wout/mdot

TE = []
profile = []
endwall = []
secondary = []
TC = []
eff = []

for n in range(5,21):
    
    delta_ho_stage = delta_ho/n
    delta_To_stage = delta_ho_stage/cp
    To_stages = [Toin - i*delta_To_stage for i in range(1,n+1)]
        
    r = (Wout/(mdot*psi*n*Omega**2))**0.5
    Vx = Omega*r*phi
    Hin = mdot/(rhoin*2*np.pi*r*Vx)
    Cxin = Hin/AR
    
    V1 = Vx/np.cos(a1)
    V2 = Vx/np.cos(a2)
    W2 = Vx/np.cos(b2)
    V3 = Vx/np.cos(a3)
    W3 = Vx/np.cos(b3)

    Tin = Toin - 0.5*V1**2/cp

    T_stators_s = [i - 0.5*V2**2/cp + delta_To_stage for i in To_stages]
    T_rotors_s = [i - 0.5*V3**2/cp for i in To_stages]
    rho_stators_s = [rhoin*(Tin/i)**(-1/(gamma-1)) for i in T_stators_s]
    rho_rotors_s = [rhoin*(Tin/i)**(-1/(gamma-1)) for i in T_rotors_s]
    
    H_rotors = [rhoin*Hin/i for i in rho_stators_s] 
    H_stators = [rhoin*Hin/i for i in rho_rotors_s]
    H_stators.insert(0,Hin)
    del H_stators[n]
    Cx_rotors = [i/AR for i in H_rotors]
    Cx_stators = [i/AR for i in H_stators]
    w_stators = [f.p_w(i,a1,a2) for i in Cx_stators]
    w_rotors = [f.p_w(i,b2,b3) for i in Cx_rotors]
    
    TEn = []
    profilen = []
    endwalln = []
    secondaryn = []
    TCn = []
    
    for i in range(n):
        
        Re_stator = rho_stators_s[i]*V2*Cx_stators[i]/mu
        Re_rotor = rho_rotors_s[i]*W3*Cx_rotors[i]/mu
        
        stator_profile = f.profile(Re_stator,a1,a2,V2,T_stators_s[i])[0]
        rotor_profile = f.profile(Re_rotor,b2,b3,W3,T_rotors_s[i])[0]
        stator_zeta = f.profile(Re_stator,a1,a2,V2,T_stators_s[i])[1]
        rotor_zeta = f.profile(Re_rotor,b2,b3,W3,T_rotors_s[i])[1]
        
        profilen.append(stator_profile+rotor_profile)      
        
        rotor_TE = f.TE(t,w_rotors[i],rotor_zeta,W3,T_rotors_s[i])
        stator_TE = f.TE(t,w_stators[i],stator_zeta,V2,T_stators_s[i])
        
        TEn.append(rotor_TE+stator_TE)
       
        if i == 0:
           inlet = 2*rhoin*CD*V1**3*2*np.pi*r*0.25*Cxin/Tin
        
        else:
            inlet = 2*rho_rotors_s[i-1]*CD*V1**3*2*np.pi*r*0.25*Cx_stators[i]/T_rotors_s[i-1]
        
        mid = 2*rho_stators_s[i]*CD*(V2**3+W2**3)*2*np.pi*r*0.25*0.5*(Cx_stators[i]+Cx_rotors[i])/T_stators_s[i]
        
        outlet = 2*rho_rotors_s[i]*CD*W3**3*2*np.pi*r*0.25*Cx_stators[i]/T_rotors_s[i]
       
        endwalln.append(inlet+mid+outlet)
        
        stator_secondary = f.secondary(a1,a2,Cx_stators[i],H_stators[i],V2,T_stators_s[i])
        rotor_secondary = f.secondary(b2,b3,Cx_rotors[i],H_rotors[i],W3,T_rotors_s[i])

        secondaryn.append(rotor_secondary+stator_secondary)
        
        tip_rotor = f.tip(g,H_rotors[i],b2,b3,W3,T_rotors_s[i])
        tip_stator = f.tip(g,H_stators[i],a1,a2,V2,T_stators_s[i])
        
        TCn.append(tip_rotor+tip_stator)
        
    TE.append(float(np.sum(TEn))*To_stages[n-1]/delta_ho)
    profile.append(float(np.sum(profilen))*To_stages[n-1]/delta_ho)
    endwall.append(float(np.sum(endwalln))*To_stages[n-1]/delta_ho)
    secondary.append(float(np.sum(secondaryn))*To_stages[n-1]/delta_ho)
    TC.append(float(np.sum(TCn))*To_stages[n-1]/delta_ho)
    
    loss = np.sum(TEn)+np.sum(profilen)+np.sum(endwalln)+np.sum(secondaryn)+np.sum(TCn)
    eff.append(delta_ho/(delta_ho+loss*To_stages[n-1]))

def turbine(Po1, To1, mdot, Omega, phi, psi, Lambda, AR, W, dho, n, t, g, ptoC=-1):
    """Return the turbine performance and sizing"""

    #Set hard limit of 73º on exit angles. Don't account for increase in psi
    #due to leakage here, can remove designs if they're unacceptable later.
    if abs(np.degrees(angles(phi, psi, Lambda)[1])) > 73 or abs(np.degrees(angles(phi, psi, Lambda)[4])) > 73:
        #Set efficiency to 0 and volume to 100 so pareto ignores them
        eff = 0
        volume = 100
        #Zeros keep indexing the same, indicate it's interrupted at the start
        return eff, 0, 0, volume, 0, 0, 0, 0, 0, 0, 0, 0

    #Gas properties
    cp = 5187
    gamma = 1.6625
    R = 2067
    gas_props = [cp, gamma, R]
    #Set the specific enthalpy drop for each stage
    del_ho = W/(mdot*n)
    #Mean radius and blade speed, constant throughout
    r = np.sqrt(del_ho/(psi*Omega**2))
    U = r*Omega
    Vx = U*phi #Axial velocity fixed by turbine parameters
    #Create size array
    sizes = [t, g, r]
    #Create array of parameters needed for the stage
    params = [mdot, AR, U, Vx, psi, phi, Lambda, ptoC]
    #Initialise quantities
    Poi = Po1
    Toi = To1
    loss = 0
    loss_comp = [0, 0, 0, 0, 0] #Entropy rise from profile, trailing edge, endwalls, secondary flow and tip leakage
    length = 0
    volume = 0
    work = 0 #To check that total work from each stage does add up to required work
    dims = [] #Array containing turbine radius, stator span and stator chord
    Yps = [] #Loss coefficients through the turbine
    m_ls = [] #Leakage mass flow through the turbine
    psis = [] #Stage loading of each stage
    angle = [] #Stage angles
    #Increment through every stage
    i = 0
    while i < n:
        #Update the input conditions for the next stage
        Poin = Poi
        Toin = Toi
        stage_calc = stage(Poin, Toin, del_ho, params, sizes, gas_props)
        #Find the exit conditions from this stage to pass to the next
        Toi = stage_calc[0]
        Poi = stage_calc[1]
        #Find the required quantities and store them
        psi_stage = stage_calc[11]
        #Check the exit angles again now that they've been properly calculated
        if abs(np.degrees(angles(phi, psi_stage, Lambda)[1])) > 73 or abs(np.degrees(angles(phi, psi_stage, Lambda)[4])) > 73:
            #Add enough loss to set efficiency effectively to 0
            loss += del_ho
            #Similar for volume
            volume += 100
            #set i to stop the while loop
            i = n+1
        psis.append(psi_stage)
        angle.append([np.degrees(i) for i in angles(phi, psi_stage, Lambda)])
        loss += stage_calc[2]
        loss_comp[0] += stage_calc[8][0] #Profile
        loss_comp[1] += stage_calc[8][1] #Trailing edge
        loss_comp[2] += stage_calc[8][2] #Endwall
        loss_comp[3] += stage_calc[8][3] #Secondary flow
        loss_comp[4] += stage_calc[8][4] #Tip leakage
        Yps.append(stage_calc[4])
        m_ls.append(stage_calc[9])
        length += stage_calc[5]
        volume += stage_calc[6]
        work += stage_calc[10]
        r = stage_calc[7][2]
        H_st = stage_calc[7][3]
        H_ro = stage_calc[7][4]
        Cx_st = stage_calc[7][5]
        Cx_ro = stage_calc[7][6]
        ptc_st = stage_calc[7][9]
        ptc_ro = stage_calc[7][10]
        dims.append([r, H_st, H_ro, Cx_st, Cx_ro, ptc_st, ptc_ro])
        #Move to the next stage
        i += 1
    #Find the efficiency using the overall loss
    eff = del_ho*n/(del_ho*n+To1*loss)
    #More outputs can be added for whatever is needed
    return eff, loss_comp, length, volume, dims, Yps, Poi, Toi, m_ls, work*mdot, psis, angle
