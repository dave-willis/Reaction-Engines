"""
Module contains functions to evaluate helium turbine performance.
Helium gas properties set in turbine function
"""
#Import required modules
import numpy as np
from scipy import interpolate as sciint

##########################
###MAIN TURBINE ROUTINE###
##########################

def turbine(Po1, To1, mdot, Omega, W, t, g, phi, psi, Lambda, AR, dho, n, ptc=-1, a1i=0):
    """Return the turbine performance and sizing"""

    #If the inputs for phi, psi, Lambda, dho, AR, or ptc ints or floats, assume they're
    #constant through the turbine
    if isinstance(phi, (float, int)):
        phi = [phi, phi]
    if isinstance(psi, (float, int)):
        psi = [psi, psi]
    if isinstance(Lambda, (float, int)):
        Lambda = [Lambda, Lambda]
    if isinstance(dho, (float, int)):
        dho = [dho, dho]
    if isinstance(AR, (float, int)):
        AR = [AR, AR]
    if isinstance(ptc, (float, int)):
        ptc = [ptc, ptc]
    #Ensure the inputs aren't longer than the number of stages, if they are cut them off
    if len(phi) > n:
        phi = phi[:n]
    if len(psi) > n:
        psi = psi[:n]
    if len(Lambda) > n:
        Lambda = Lambda[:n]
    if len(dho) > n:
        dho = dho[:n]
    if len(AR) > n:
        AR = AR[:n]
    if len(ptc) > n:
        ptc = ptc[:n]
    #Store the inputs for later use
    inits = [Po1, To1, mdot, Omega, W, t, g, phi, psi, Lambda, AR, dho, n, ptc, a1i]
    #If all loadings are constant treat as repeating stages
    rep = False
    if all(i == phi[0] for i in phi) and all(i == psi[0] for i in psi) and all(i == Lambda[0] for i in Lambda):
        find_angs = repeating
        rep = True
    else:
        find_angs = angles
    #Create splines of the stage loadings if n>1
    if n > 1:
        phi = spline(n, phi)
        psi = spline(n, psi)
        Lambda = spline(n, Lambda)
        dho = spline(n, dho)
        AR = spline(n, AR)
        ptc = spline(n, ptc)
        #Scale the enthalpy drops to ensure the total work is as desired
        sum_dho = sum(dho)
        dho = [W/(mdot*sum_dho)*i for i in dho]
    #If n=1 then scale the first vale of dho accordingly
    if n == 1:
        dho = [W/mdot, 0]
    #Set hard limit of 73º on exit angles. Don't account for increase in psi
    #due to leakage here, can remove designs if they're unacceptable later.
    ang_1 = 0 #Initialise the inlet angle
    angle_warning = [False, 0] #Initialise the angle warning
    for l in range(len(phi)):
        ang_check = find_angs(phi[l], psi[l], Lambda[l], ang_1)
        ang_2 = np.degrees(ang_check[1])
        ang_3 = np.degrees(ang_check[4])
        ang_1 = ang_check[3]
        if abs(ang_2) > 73 or abs(ang_3) > 73:
            angle_warning = [True, max(abs(ang_2), abs(ang_3))]
        elif angle_warning[0] == False:
            angle_warning = [False, max(abs(ang_2), abs(ang_3))]
    #Gas properties
    cp = 5187
    gamma = 1.6625
    R = 2067
    gas_props = [cp, gamma, R]
    #Overall specific enthalpy drop
    del_ho = W/mdot
    #Initialise quantities
    Poi = Po1
    Toi = To1
    a1i = np.radians(a1i)
    loss = 0
    loss_comp = [0, 0, 0, 0, 0] #Entropy rise from profile, TE, EW, secondary flow and tip
    length = 0
    volume = 0
    mass = 0
    work = 0 #To check that total work from each stage does add up to required work
    dims = [] #Array containing turbine radius, stator span and stator chord
    angle = [] #Stage angles
    n_blades = 0 #Total blades in the turbine
    #Increment through every stage
    i = 0
    while i < n:
        #Mean radius and blade speed, constant throughout
        r = np.sqrt(dho[i]/(psi[i]*Omega**2))
        U = r*Omega
        Vx = U*phi[i] #Axial velocity fixed by turbine parameters
        #Create size array
        sizes = [t, g, r]
        #Create array of parameters needed for the stage
        params = [mdot, U, Vx, psi[i], phi[i], Lambda[i], a1i, AR[i], ptc[i], rep]
        #Update the input conditions for the next stage
        Poin = Poi
        Toin = Toi
        stage_calc = stage(Poin, Toin, dho[i], params, sizes, gas_props)
        #Store angles before updating a1
        angle.append([np.degrees(x) for x in find_angs(phi[i], psi[i], Lambda[i], a1i)])
        #Find the exit conditions from this stage to pass to the next
        Toi = stage_calc[0]
        Poi = stage_calc[1]
        a1i = stage_calc[10]
        #Find the required quantities and store them
        loss += stage_calc[8]
        loss_comp[0] += stage_calc[9][0] #Profile
        loss_comp[1] += stage_calc[9][1] #Trailing edge
        loss_comp[2] += stage_calc[9][2] #Endwall
        loss_comp[3] += stage_calc[9][3] #Secondary flow
        loss_comp[4] += stage_calc[9][4] #Tip leakage
        length += stage_calc[6]
        volume += stage_calc[4]
        mass += stage_calc[3]
        work += stage_calc[5]
        r = stage_calc[7][2]
        n_blades += stage_calc[11]
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
    eff = del_ho/(del_ho+To1*loss)
    #More outputs can be added for whatever is needed
    return eff, work*mdot, mass, volume, length, dims, n_blades, loss_comp, Poi, Toi, angle, inits, angle_warning

###################
###STAGE ROUTINE###
###################

def stage(Po1, To1, del_ho, params, sizes, gas_props):
    """Return the stage losses and exit conditions"""

    #Gas properties
    cp, gamma, R = gas_props
    #Stage parameters
    mdot, U, Vx, psi, phi, Lambda, a1, AR, ptc, rep = params
    #Set the angle function:
    if rep:
        find_angs = repeating
    else:
        find_angs = angles
    #Major sizes
    t, g, r = sizes
    #Initial guess for the stage work output accounting for tip leakage
    work = del_ho#*0.95
    #Iterate until the stage produces the requird work
#    while work/del_ho < 0.999 or work/del_ho > 1.001:
#        psi = psi*del_ho/work
    #Calculate angles
    a1, a2, b2, a3, b3 = find_angs(phi, psi, Lambda, a1)
    #Create array of angles needed in the loss function
    angs = [a1, a2, b2, b3]
    #Calculate velocities
    V1, V2, W2, V3, W3 = velocities([a1, a2, b2, a3, b3], Vx)
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
        states = [T1, T2, T3, rho1, rho2, rho3, mu2, mu3, mdot]
        #Stage sizing
        H1 = mdot/(rho1*2*np.pi*r*Vx)
        H2 = mdot/(rho2*2*np.pi*r*Vx)
        H3 = mdot/(rho3*2*np.pi*r*Vx)
        #Use average blade span to find axial chord and throat width
        H_st = (H1+H2)/2
        H_ro = (H2+H3)/2
        Cx_st = H_st/AR
        Cx_ro = H_ro/AR
        w_st = p_w(Cx_st, a1, a2, ptc)[0]
        w_ro = p_w(Cx_ro, b2, b3, ptc)[0]
        #Create array for the loss function
        dimensions = [t, g, r, H_st, H_ro, Cx_st, Cx_ro, w_st, w_ro]
        #Find the losses across the stator and rotor
        loss_calc = losses(angs, vels, states, dimensions)
        loss_st = loss_calc[1] #Stator entropy rise
        loss_ro = loss_calc[2] #Rotor entropy rise
        #Use entropy rise to find stagnation pressure loss coefficients
        #assuming entropy loss coefficient is equivalent at low Mach numbers
        Y_st_guess = loss_st*T2/(0.5*V2**2)
        Y_ro_guess = loss_ro*T3/(0.5*W3**2)
    #Find the actual work output given the rotor leakage
    m_l = np.amax(loss_calc[4])
    work = psi*(1-m_l)*U**2
    #Stage performance
    loss_array = loss_calc[3]
    loss = loss_calc[0]
    eff = del_ho/(del_ho+To3*loss)
    #Stage sizes
    length = 1.5*Cx_st + 1.5*Cx_ro
    vol_st = np.pi*Cx_st*(0.25*((r+H1/2)**2)+((r+(H1+H2)/4)**2)+0.25*(r+H2/2)**2)
    vol_ro = np.pi*Cx_ro*(0.25*((r+H2/2)**2)+((r+(H2+H3)/4)**2)+0.25*(r+H3/2)**2)
    volume = vol_st + vol_ro
    ptc_st = p_w(Cx_st, a1, a2, ptc)[1]
    ptc_ro = p_w(Cx_ro, b2, b3, ptc)[1]
    dimensions.append(ptc_st)
    dimensions.append(ptc_ro)
    #Stage mass
    m_vessel_st = (vessel_mass(P1, r+H1/2, 1.5*Cx_st)+vessel_mass(P2, r+H2/2, 1.5*Cx_st))/2
    m_vessel_ro = (vessel_mass(P2, r+H2/2, 1.5*Cx_ro)+vessel_mass(P3, r+H3/2, 1.5*Cx_ro))/2
    m_drum = drum_mass(r-H_st/2, 1.5*Cx_st)+drum_mass(r-H_ro/2, 1.5*Cx_ro)
    m_blade = blade_mass(r, ptc_st, H_st, Cx_st)[0]+blade_mass(r, ptc_ro, H_ro, Cx_ro)[0]
    mass = m_vessel_st+m_vessel_ro+m_drum+m_blade
    #Total number of blades in the stage
    n_blades = blade_mass(r, ptc_st, H_st, Cx_st)[1]+blade_mass(r, ptc_ro, H_ro, Cx_ro)[1]

    return To3, Po3, eff, mass, volume, work, length, dimensions, loss, loss_array, a3, n_blades

######################
###VORTEX FUNCTIONS###
######################

def free_vortex(angs, size, phi):
    """Return the hub, midspan and tip metal angles for a free vortex stage"""

    #Pull angles out of the array, converting them to radians
    a1, a2, b2, a3, b3 = np.radians(angs)
    #Pull out radius and spans
    r, hst, hro = size
    #Free vortex has constant Vx so keeping rVt constant is simple
    a1h = np.arctan(r*np.tan(a1)/(r-hst/2))
    a1t = np.arctan(r*np.tan(a1)/(r+hst/2))
    a2h = np.arctan(r*np.tan(a2)/(r-hst/2))
    a2t = np.arctan(r*np.tan(a2)/(r+hst/2))
    #New relative flow angles depend on the flow coefficient at the span
    b2h = np.arctan(np.tan(a2h)-(r-hro/2)/(phi*r))
    b2t = np.arctan(np.tan(a2t)-(r+hro/2)/(phi*r))
    b3h = np.arctan(np.tan(a1h)-(r-hro/2)/(phi*r))
    b3t = np.arctan(np.tan(a1t)-(r+hro/2)/(phi*r))
    a3 = a3

    return np.degrees([a1h, a1, a1t]), np.degrees([a2h, a2, a2t]), np.degrees([b2h, b2, b2t]), np.degrees([b3h, b3, b3t])

####################
###LOSS FUNCTIONS###
####################

def Y_po(po1, Y, V, c, gamma):
    """Return the exit stagnation pressure after loss"""

    M = V/c #Flow Mach number
    fac = (1+(gamma-1)*(M**2)/2)**(-gamma/(gamma-1)) #Makes next line neater
    #Find new stagnation pressure
    po2 = po1/(Y+1-Y*fac)

    return po2

def losses(angs, vels, states, dimensions):
    """Return the entropy rises for given inputs in Stage"""

    #Create values from input arrays
    a1, a2, b2, b3 = angs
    V1, V2, W2, W3 = vels
    T1, T2, T3, rho1, rho2, rho3, mu2, mu3, mdot = states
    t, g, r, H_st, H_ro, Cx_st, Cx_ro, w_st, w_ro = dimensions
    #Calculate tip clearance effect first so that mass flow can be
    #accounted for in profile and trailing edge loss
    TC_calcs_st = TC(g, H_st, a1, a2, V2, T2)
    TC_st = TC_calcs_st[0] #Entropy rise
    m_st = TC_calcs_st[1] #Stator leakage mass flow fraction
    TC_calcs_ro = TC(g, H_ro, b2, b3, W3, T3)
    TC_ro = TC_calcs_ro[0]
    m_ro = TC_calcs_ro[1]
    TC_loss = TC_st+TC_ro #Overall tip leakage entropy rise
    m_ls = [m_st, m_ro] #Leakage mass flow fraction for use in analysis
    #Reynolds numbers calculated using exit conditions
    Re_st = (rho2)*(V2)*Cx_st/mu2/4
    Re_ro = (rho3)*(W3)*Cx_ro/mu3/4
    #Profile loss calculations
    profile_calcs_st = profile(Re_st, a1, a2, V2, T2)
    profile_calcs_ro = profile(Re_ro, b2, b3, W3, T3)
    profile_st = profile_calcs_st[0]*(1-m_st)
    profile_ro = profile_calcs_ro[0]*(1-m_ro)
    profile_loss = profile_st+profile_ro
    #Find the profile loss coefficients for the momentum boundary layer thickness
    zeta_st = profile_calcs_st[1]
    zeta_ro = profile_calcs_ro[1]
    #Trailing edge loss calculations
    TE_st = TE(t, w_st, zeta_st, V2, T2)*(1-m_st)
    TE_ro = TE(t, w_ro, zeta_ro, W3, T3)*(1-m_ro)
    TE_loss = TE_st+TE_ro
    #Endwall loss calculations account for varying chord and conditions
    EW_st = 0#EW(rho1, V1, r, Cx_st, T1, mdot)+EW(rho2, V2, r, Cx_st, T2, mdot)
    EW_ro = 0#EW(rho2, W2, r, Cx_ro, T2, mdot)+EW(rho3, W3, r, Cx_ro, T3, mdot)
    EW_loss = EW_st+EW_ro
    #Secondary flow loss calculations
    secondary_st = secondary(a1, a2, Cx_st, H_st, V2, T2)
    secondary_ro = secondary(b2, b3, Cx_ro, H_ro, W3, T3)
    secondary_loss = secondary_st+secondary_ro
    #Add up entropy rises for the stator and rotor and find the overall loss
    loss_st = profile_st+TE_st+EW_st+secondary_st+TC_st
    loss_ro = profile_ro+TE_ro+EW_ro+secondary_ro+TC_ro
    loss = loss_st+loss_ro
    loss_comp = [profile_loss, TE_loss, EW_loss, secondary_loss, TC_loss] #Array of loss components

    return loss, loss_st, loss_ro, loss_comp, m_ls

def profile(Re, a1, a2, V, T):
    """Return the profile entropy rise for the stage"""

    v_frac = 1/np.sqrt(3) #Fix the idealised velocity fraction
    #Calculate dissipation coefficient using correlation
    Cd = 0.002*(Re/500000)**(-0.2)
    #Calculate entropy loss coefficient for use in trailing edge loss function
    zeta = Cd*(2/v_frac+6*v_frac)*abs(np.tan(a2)-np.tan(a1))
    #Calculate entropy rise
    entropy = (zeta*0.5*V**2)/T

    return entropy, zeta

def TE(t, w, zeta_profile, V, T):
    """Return the trailing edge entropy rise for the stage"""

    SF = 1.4 #Fix the shape function for a turbulent boundary layer
    Cpb = -0.15 #Fix the base pressure coefficient
    #Calculate the momentum boundary layer thickness using profile loss coefficient
    BL_mom = zeta_profile*w/2
    #Calculate the displacement boundary layer thickness
    BL_dis = SF*BL_mom
    #Calculate the entropy loss coefficient
    zeta = ((t+BL_dis)/w)**2-Cpb*t/w
    #Calculate the entropy rise
    entropy = (zeta*0.5*V**2)/T

    return entropy

def EW(rho, V, r, Cx, T, mdot):
    """Return the endwall entropy rise for the stage"""

    CD = 0.002 #Fix endwall dissipation coefficient
    #Calculate total entropy rise, factor of 2 accounts for hub and casing
    entropy = 2*rho*CD*(V**3)*2*np.pi*r*0.25*Cx/T
    #Convert to specific entropy
    entropy_sp = entropy/mdot

    return entropy_sp

def secondary(a1, a2, Cx, H, V, T):
    """Return the secondary flow entropy rise for the stage"""

    #Calculate the vector mean air angle
    am = np.arctan(0.5*(np.tan(a1)+np.tan(a2)))
    #Calculate the stagnation pressure loss coefficient
    Y = 0.375*0.1336*(Cx/H)*(np.cos(a2)**3)*((np.tan(a1)-np.tan(a2))**2)/(np.sqrt(np.cos(a1))*np.cos(am))
    #Convert to entropy rise assuming incompressible flow
    entropy = (Y*0.5*V**2)/T

    return entropy

def TC(g, H, a1, a2, V2, T):
    """Return the tip leakage entropy rise for the stage"""

    Cc = 0.6 #Fix the shroud contraction coefficient
    #Calculate the leakage mass flow fraction
    m = g*Cc*np.sqrt(abs((1/np.cos(a2))**2-np.tan(a1)**2))/H
    #Calculate the entropy rise
    entropy = m*(V2**2)*(1-np.tan(a1)*(np.sin(a2)**2)/np.tan(a2))/T

    return entropy, m

########################
###GEOMETRY FUNCTIONS###
########################

def angles(phi, psi, Lambda, a1):
    """Return the flow angles for the specified loading and non-repeating stages"""

    #Define coefficients in the quadratic for tan(a3) to keep it clear
    a = Lambda*phi**2
    b = 2*psi*phi
    c = 2*psi*Lambda-Lambda*(phi*np.tan(np.radians(a1)))**2-2*psi+psi**2
    #Calculate angles, using positive root for a3
    a3 = np.arctan((-b+np.sqrt(b**2-4*a*c))/(2*a))
    a2 = np.arctan(np.tan(a3)+psi/phi)
    b2 = np.arctan(np.tan(a2)-1/phi)
    b3 = np.arctan(np.tan(a3)-1/phi)

    return a1, a2, b2, a3, b3

def repeating(phi, psi, Lambda, a1):
    """Return the flow angles for specified loading and repeating stage"""

    #Calculate angles
    a1 = np.arctan((1-0.5*psi-Lambda)/phi)
    a2 = np.arctan(np.tan(a1)+psi/phi)
    b2 = np.arctan(np.tan(a2)-1/phi)
    b3 = np.arctan(np.tan(a1)-1/phi)
    a3 = a1

    return a1, a2, b2, a3, b3

def velocities(angs, Vx):
    """ Return the velocities for the specified angles"""

    #Extract angles from array
    a1, a2, b2, a3, b3 = [i for i in angs]
    #Calculate velocities
    V1 = Vx/np.cos(a1)
    V2 = Vx/np.cos(a2)
    W2 = Vx/np.cos(b2)
    V3 = Vx/np.cos(a3)
    W3 = Vx/np.cos(b3)

    return V1, V2, W2, V3, W3


def p_w(Cx, a1, a2, ptoC):
    """Return the pitch and throat width for the stage"""
    #If pitch to chord is specified, use it
    if ptoC > 0:
        ptc = ptoC
        p = Cx*ptc
    #Otherwise use Zweifel criterion
    else:
        Z = 0.8 #Fix the Zweifel coefficient
        #Calculate pitch to chord ratio
        ptc = 0.5*Z/((np.cos(a2)**2)*abs(np.tan(a2)-np.tan(a1)))
        #Calculate pitch
        p = Cx*ptc
    #Calculate throat width
    w = p*np.cos(a2)

    return w, ptc

####################
###MASS FUNCTIONS###
####################

def vessel_mass(Pin, r, l):
    """Return the pressure vessel mass required using Tresca"""

    sigma_yield = 500*10**6 #Yield stress of steel in Pa
    rho_m = 7840 #Density of steel in kg/m^3
    safety_fac = 1.5 #Use a safety factor
    #Calculate wall thickness required with hoop stress
    t = safety_fac*Pin*r/sigma_yield
    #Calculate the mass of the vessel assuming thin walled
    m = rho_m*2*np.pi*r*l*t

    return m

def drum_mass(r, l):
    """Return the mass of the drum for the blade row"""

    rho_m = 7840 #Density of steel in kg/m^3
    thk = 0.1*r #Drum wall thickness
    #Calculate the mass of the drum
    m = rho_m*l*2*np.pi*r*thk

    return m

def blade_mass(rm, ptc, H, Cx):
    """Return the mass of the blades in the blade row"""

    rho_m = 8970 #Density of Inconel 188 in kg/m^3
    t_av = Cx*0.1 #Define an average blade thickness relative to the chord
    #Calculate the total number of blades
    n = round(2*np.pi*rm/(ptc*Cx))
    #Approximate blade as flate plate to calculate mass
    m = rho_m*H*Cx*t_av*n

    return m, n

######################
###SPLINE FUNCTION###
######################

def spline(n, y):
    """Return spline of loadings using specified control points"""

    #Array of stages
    stages = [0]*len(y)
    stages[0] = 1
    stages[-1] = n
    if len(stages) > 2:
        for i in range(1, len(stages)-1):
            stages[i] = 1+i*(n-1)/(len(stages)-1)
    order = len(stages)-1
    #Spline the points
    tck = sciint.splrep(stages, y, k=order)
    xnew = np.arange(1, n+1, 1)
    ynew = sciint.splev(xnew, tck)

    return ynew
