"""
Module contains functions to evaluate helium turbine performance.
Helium gas properties set in turbine function, would require new values of
cp and gamma and a new viscosity law
"""
#Import required modules
import numpy as np
from scipy import interpolate as sciint
import scipy.integrate as integrate
from scipy.optimize import minimize
from Profile import blade_dims
import csv

#Declare a few global options
print_warnings = False #If True, will print warnings when material limits are exceeded
fast_dimensions = True #If True, will interpolate from tables for blade dimensions

##########################
###MAIN TURBINE ROUTINE###
##########################

def turbine(Po1, To1, mdot, Omega, W, t, g, phi, psi, Lambda, AR, dho, n, ptc=-1, ain=0, gas='He'):
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
    #If all loadings are constant treat as repeating stages
    rep = False
    if all(i == phi[0] for i in phi) and all(i == psi[0] for i in psi) and all(i == Lambda[0] for i in Lambda):
        find_angs = repeating
        rep = True
        ain = np.degrees(find_angs(phi[0], psi[0], Lambda[0], ain)[0])
    else:
        find_angs = angles
    #Store the inputs for later use
    inits = [Po1, To1, mdot, Omega, W, t, g, phi, psi, Lambda, AR, dho, n, ptc, ain, gas]
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
        if np.round(abs(ang_2), 2) > 73.0 or np.round(abs(ang_3), 2) > 73.0:
            angle_warning = [True, max(abs(ang_2), abs(ang_3))]
        elif not angle_warning[0]:
            angle_warning = [False, max(abs(ang_2), abs(ang_3))]
    #Rotor and stator materials
    rotor_mat = 'inconel'
    stator_mat = 'inconel'
    materials = [stator_mat, rotor_mat]
    #Overall specific enthalpy drop
    del_ho = W/mdot
    #Initialise quantities
    Poi = Po1
    Toi = To1
    a1i = np.radians(ain)
    loss = 0
    loss_comp = [0, 0, 0, 0] #Entropy rise from profile, TE, secondary flow and tip
    length = 0
    volume = 0
    mass = 0
    work = 0 #To check that total work from each stage does add up to required work
    dims = [] #Array containing turbine radius, stator span and stator chord
    angle = [] #Stage angles
    n_blades = 0 #Total blades in the turbine
    Fx = 0 #Axial force on rotor
    Re = 0 #Average Reynolds number
    expansion_lims = [] #Stator and rotor expansion, relative to static-cold, and thermal time
    #Select the function for calculating blade dimensions depending on speed wanted
    global fast_dimensions
    dim_fs = []
    if fast_dimensions:
        #Load the grid of normalised suction surface lengths and create the interpolation function
        xin = []
        xout = []
        sslen = []
        areas = []
        with open("ss_grid.csv") as csvfile:
            reader = csv.reader(csvfile, quoting=csv.QUOTE_NONNUMERIC) # change contents to floats
            line = 0
            for row in reader: # each row is a list
                if line == 0:
                    xout = row[1:]
                    line += 1
                else:
                    xin.append(row[0])
                    sslen.append(row[1:])
                    line += 1
        with open("Profile_areas.csv") as csvfile:
            reader = csv.reader(csvfile, quoting=csv.QUOTE_NONNUMERIC) # change contents to floats
            line = 0
            for row in reader: # each row is a list
                if line == 0:
                    line += 1
                else:
                    areas.append(row[1:])
                    line += 1
        #Make the lists numpy arrays for the interpolation
        xout = np.asarray(xout)
        xin = np.asarray(xin)
        sslen = np.asarray(sslen)
        areas = np.asarray(areas)
        #Make the interpolation functions
        ss_length = sciint.RectBivariateSpline(xin, xout, sslen, kx=1, ky=1)
        profile_area = sciint.RectBivariateSpline(xin, xout, areas, kx=1, ky=1)
        dim_fs = [ss_length, profile_area]
    #Increment through every stage
    i = 0
    while i < n:
        #Mean radius and blade speed, constant throughout
        r = np.sqrt(dho[i]/(psi[i]*Omega**2))
        U = r*Omega
        Vx = U*phi[i] #Axial velocity fixed by turbine parameters
        #Create size array, including the function for ss length
        sizes = [t, g, r, dim_fs]
        #Gas properties
        gamma, cp, R = thermo_props(Toi, gas)
        gas_props = [cp, gamma, R, gas]
        #Create array of parameters needed for the stage
        params = [mdot, Omega, U, Vx, psi[i], phi[i], Lambda[i], a1i, AR[i], ptc[i], rep, i+1, n]
        #Update the input conditions for the next stage
        Poin = Poi
        Toin = Toi
        stage_calc = stage(Poin, Toin, dho[i], params, sizes, gas_props, materials)
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
        loss_comp[2] += stage_calc[9][2] #Secondary flow
        loss_comp[3] += stage_calc[9][3] #Tip leakage
        length += stage_calc[6]
        volume += stage_calc[4]
        mass += stage_calc[3]
        work += stage_calc[5]
        r = stage_calc[7][2]
        n_blades += stage_calc[11]
        Fx += stage_calc[12]
        H_st = stage_calc[7][3]
        H_ro = stage_calc[7][4]
        Cx_st = stage_calc[7][5]
        Cx_ro = stage_calc[7][6]
        ptc_st = stage_calc[7][9]
        ptc_ro = stage_calc[7][10]
        Ro_st = stage_calc[7][11]
        Ri_ro = stage_calc[7][12]
        dims.append([r, H_st, H_ro, Cx_st, Cx_ro, ptc_st, ptc_ro, Ro_st, Ri_ro])
        Re += (stage_calc[13][0]+stage_calc[13][1])/(2*n)
        expansion_lims.append(stage_calc[14])
        #Move to the next stage
        i += 1
    #Find the efficiency using the overall loss
    eff = del_ho/(del_ho+To1*loss)
    #More outputs can be added for whatever is needed
    return eff, work*mdot, mass, volume, length, dims, n_blades, loss_comp, Poi, Toi, angle, inits, angle_warning, Fx, Re, expansion_lims

###################
###STAGE ROUTINE###
###################

def stage(Po1, To1, del_ho, params, sizes, gas_props, materials):
    """Return the stage losses and exit conditions"""

    global fast_dimensions, print_warnings
    #Gas properties
    cp, gamma, R, gas = gas_props
    #Stage parameters
    mdot, Omega, U, Vx, psi, phi, Lambda, a1, AR, ptc, rep, stage_n, n_stages = params
    #Set the angle function:
    if rep:
        find_angs = repeating
    else:
        find_angs = angles
    #Major sizes
    t, g, r, dim_fs = sizes
    #Initial guess for the stage work output accounting for tip leakage
    work = del_ho#*0.95
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
    mu2 = viscosity(T2, gas)
    To3 = To2 - del_ho/cp #Stagnation temperature drop depending on stage work
    T3 = To3 - 0.5*(V3**2)/cp
    c3 = np.sqrt(gamma*R*T3)
    mu3 = viscosity(T3, gas)
    #Use compressible relation to find static pressure
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
        #Calculate the Reynolds numbers
        #Run the profile generator for SS length and section area
        if fast_dimensions:
            ss_st = Cx_st*dim_fs[0](np.degrees(a1), np.degrees(a2))[0][0]
            ss_ro = Cx_ro*dim_fs[0](np.degrees(b2), np.degrees(b3))[0][0]
        else:
            stator_dims = blade_dims(np.degrees(a1), np.degrees(a2), t, Cx_st)
            rotor_dims = blade_dims(np.degrees(b2), np.degrees(b3), t, Cx_ro)
            stator_A = stator_dims[0]
            ss_st = stator_dims[1]
            rotor_A = rotor_dims[0]
            ss_ro = rotor_dims[1]
        #Calculate Re based on SS length
        Res = [rho2*V2*ss_st/mu2, rho3*W3*ss_ro/mu3]
        #Find the losses across the stator and rotor
        loss_calc = losses(angs, vels, states, dimensions, Res)
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
    #Stage mass and blades
    if fast_dimensions:
        stator_A = Cx_st**2*dim_fs[1](np.degrees(a1), np.degrees(a2))[0][0]
        rotor_A = Cx_ro**2*dim_fs[1](np.degrees(b2), np.degrees(b3))[0][0]
    blade_areas = [stator_A, rotor_A]
    mass_calc = stage_mass(To1, Po1, Po3, dimensions, Omega, blade_areas, materials, stage_n, n_stages)
    mass = mass_calc[0]
    n_blades = mass_calc[1]
    Ro_st = mass_calc[2]
    Ri_ro = mass_calc[3]
    dimensions.append(Ro_st)
    dimensions.append(Ri_ro)
    #Axial force on rotor
    Fx = blade_force(P2, P3, r, H2, H3)
    #Thermal start up
    vels.append(V3)
    pressures = [P1, P2, P3]
    thermal_calc = start_up(vels, states, pressures, dimensions, materials, gas)
    dg_c = mass_calc[4][0]-mass_calc[4][1]
    dg_h = mass_calc[5][0]-mass_calc[5][1]
    #The heating time constants can be used to get an estimate of a worst
    #case scenrio where one component is fully heated and the other only
    #heated to some fraction governed by the relative heating times
    if thermal_calc[0] < thermal_calc[1]:
        #Scenario where stator heats up faster
        dr_st = mass_calc[5][0]
        dr_ro = mass_calc[4][1]+thermal_calc[0]/thermal_calc[1]*(mass_calc[5][1]-mass_calc[4][1])
        dg_warm = dr_st-dr_ro
    if thermal_calc[0] >= thermal_calc[1]:
        #Scenario where rotor heats up faster
        dr_ro = mass_calc[5][1]
        dr_st = mass_calc[4][0]+thermal_calc[1]/thermal_calc[0]*(mass_calc[5][0]-mass_calc[4][0])
        dg_warm = dr_st-dr_ro
    #Depending on the relative values of the hot and cold strains, the tip
    #clearance when cold-static, cold-rotating warm-rotating hot-rotating will change
    if dg_c >= 0 and dg_h >= 0:
        cold_stat_g = 0.0
        cold_rot_g = dg_c
        warm_rot_g = dg_warm
        hot_rot_g = dg_h
    elif dg_c >= 0 and dg_h < 0:
        cold_stat_g = abs(dg_h)
        cold_rot_g = dg_c+abs(dg_h)
        warm_rot_g = dg_warm+abs(dg_h)
        hot_rot_g = 0.0
    elif dg_h >= 0 and dg_c < 0:
        cold_stat_g = abs(dg_c)
        cold_rot_g = 0.0
        warm_rot_g = dg_warm+abs(dg_c)
        hot_rot_g = dg_h+abs(dg_c)
    elif dg_c < 0 and dg_h < dg_c:
        cold_stat_g = abs(dg_h)
        cold_rot_g = dg_c-dg_h
        warm_rot_g = dg_warm+abs(dg_h)
        hot_rot_g = 0.0
    elif dg_h < 0 and dg_c <= dg_h:
        cold_stat_g = abs(dg_c)
        cold_rot_g = 0.0
        warm_rot_g = dg_warm+abs(dg_c)
        hot_rot_g = dg_h-dg_c
    if warm_rot_g < 0:
        if print_warnings:
            print('WARNING: STAGE {}: POSSIBLE ROTOR-STATOR CONTACT IN TRANSIENT HEATING'.format(stage_n))
    expansion_lims = cold_stat_g, cold_rot_g, warm_rot_g, hot_rot_g

    return To3, Po3, eff, mass, volume, work, length, dimensions, loss, loss_array, a3, n_blades, Fx, Res, expansion_lims

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

def losses(angs, vels, states, dimensions, Res):
    """Return the entropy rises for given inputs in Stage"""

    #Create values from input arrays
    a1, a2, b2, b3 = angs
    V1, V2, W2, W3 = vels
    T1, T2, T3, rho1, rho2, rho3, mu2, mu3, mdot = states
    t, g, r, H_st, H_ro, Cx_st, Cx_ro, w_st, w_ro = dimensions
    Re_st, Re_ro = Res
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
    #Secondary flow loss calculations
    secondary_st = secondary(a1, a2, Cx_st, H_st, V2, T2)
    secondary_ro = secondary(b2, b3, Cx_ro, H_ro, W3, T3)
    secondary_loss = secondary_st+secondary_ro
    #Add up entropy rises for the stator and rotor and find the overall loss
    loss_st = profile_st+TE_st+secondary_st+1.0*TC_st
    loss_ro = profile_ro+TE_ro+secondary_ro+TC_ro
    loss = loss_st+loss_ro
    loss_comp = [profile_loss, TE_loss, secondary_loss, TC_loss] #Array of loss components

    return loss, loss_st, loss_ro, loss_comp, m_ls

def profile(Re, a1, a2, V, T):
    """Return the profile entropy rise for the stage"""

    #Denton 1993
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

    #Denton 1993
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

def secondary(a1, a2, Cx, H, V, T):
    """Return the secondary flow entropy rise for the stage"""

    #Dunham & Came 1970
    #Calculate the vector mean air angle
    am = np.arctan(0.5*(np.tan(a1)+np.tan(a2)))
    #Calculate the stagnation pressure loss coefficient
    Y = 0.375*0.1336*(Cx/H)*(np.cos(a2)**3)*((np.tan(a1)-np.tan(a2))**2)/(np.sqrt(np.cos(a1))*np.cos(am))
    #Convert to entropy rise assuming incompressible flow
    entropy = (Y*0.5*V**2)/T

    return entropy

def TC(g, H, a1, a2, V2, T):
    """Return the tip leakage entropy rise for the stage"""

    #Denton 1993
    Cc = 0.6 #Fix the shroud contraction coefficient
    #Calculate the leakage mass flow fraction
    m = g*Cc*np.sqrt(abs((1/np.cos(a2))**2-np.tan(a1)**2))/H
    #Calculate the entropy rise
    entropy = m*(V2**2)*(1-np.tan(a1)*(np.sin(a2)**2)/np.tan(a2))/T

    return entropy, m

####################
###GAS PROPERTIES###
####################

def thermo_props(T, gas):
    """Calculate gamma, cp and R for the gas"""

    #Helium
    if gas == 'He':
        gamma = 1.6625
        cp = 5187
        R = cp-cp/gamma
    #Combustion gases with A1 jet fuel (AFR = 50)
    if gas == 'A1':
        gamma = 1.41-8.49*10**(-5)*T
        cp = 0.22*T+951
        R = cp-cp/gamma

    return gamma, cp, R

def viscosity(T, gas):
    """Calculate viscosity of the gas"""

    #Helium
    if gas == 'He':
        mu = (3.674*10**(-7))*(T**0.7)
    #Combustion gases with A1 jet fuel (AFR = 50)
    if gas == 'A1':
        mu = (1.31*10**(-5)+0.0059*T-1.71*10**(-6)*T**2)*10**(-5)

    return mu

def Prandtl(T, P, gas):
    """Calculate Prandtl number of the gas"""

    #Helium
    if gas == 'He':
        Pr = 0.6728/(1+P*2.7*10**(-9))*(T/273)**(-(0.01-P*1.42*10**(-9)))
    #Combustion gases with A1 jet fuel (AFR = 50)
    if gas == 'A1':
        Pr = 0.716-7.66*10**(-6)*T

    return Pr

def conductivity(T, gas):
    """Calculate thermal conductivity of the gas"""

    #Helium
    if gas == 'He':
        k = 0.14789*(T/273)**0.6958
    #Combustion gases with A1 jet fuel (AFR = 50)
    if gas == 'A1':
        k = (0.0059*T+0.94)*10**(-2)

    return k

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
    a1, a2, b2, a3, b3 = angs
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

########################################
###MASS, THERMAL AND STRESS FUNCTIONS###
########################################

def steel(T):
    """Return the properties of steel at temperature T"""

    #Poission's ratio (Eurocode 3)
    nu_m = 0.3
    #Density (Eurocode 3)
    rho_m = 7850
    #Thermal strain (Eurocode 3)
    strain = 1.2*10**(-5)*(T-273)+4*10**(-9)*(T-273)**2-2.416*10**(-4)
    #Specific heat capacity (Eurocode 3)
    if T < 873:
        c_m = 425+7.73*10**(-1)*(T-273)-1.69*10**(-3)*(T-273)**2+2.22*10**(-6)*(T-273)**3
    else:
        c_m = 666+13002/(1011-T)
    #Thermal conductivity (Eurocode 3)
    k_m = 54-3.33*(10**(-2)*(T-273))
    #Yield stress (Wang, Liu & Kodur 2013)
    if T < 723:
        sigma_m = 500*10**6
    else:
        sigma_m = (4.32*np.exp(-(T-273)/880)-1.6)*500*10**6
    #Young's modulus
    E_m = (1.02-0.035*np.exp((T-273)/280))*210*10**9

    return nu_m, rho_m, strain, c_m, k_m, sigma_m, E_m

def inconel(T):
    """Return the properties of Inconel 718 at temperature T"""

    #Poission's ratio (Díaz-Álvarez et al. 2017)
    nu_m = 0.3
    #Density (Díaz-Álvarez et al. 2017)
    rho_m = 8300
    #Thermal strain (Lewandowski & Overfelt 1999)
    strain = integrate.quad(lambda x: (1.28+5*10**(-4)*(x-366))*10**(-5), 293, T)
    #Specific heat capacity (Díaz-Álvarez et al. 2017)
    c_m = 400+0.15*np.exp((T-300)/90)
    #Thermal conductivity (Díaz-Álvarez et al. 2017)
    k_m = (1/60)*T+5
    #Yield stress (Woo & Lee 2019)
    sigma_m = (1250-np.exp((T-273)/128))*10**6
    #Young's modulus (Thomas et el. 2006)
    E_m = 2*(1+nu_m)*(1-0.5*(T-300)/1673)*8.31*10**10

    return nu_m, rho_m, strain[0], c_m, k_m, sigma_m, E_m

def blade_mass(blade_sizes, rho_m, row_type, shroud=True):
    """Return the mass of the blades in the blade row"""

    rm, ptc, H, Cx, A = blade_sizes
    #Calculate the total number of blades
    N = round(2*np.pi*rm/(ptc*Cx))
    #Calculate blade mass
    m = rho_m*H*A*N
    #Include shroud mass if specified
    m_shroud = 0
    if shroud:
        #Shroud thickness of 1.5mm
        t_shroud = 0.001
        #Calculation different on rotor or stator
        if row_type == 'rotor':
            #Single blade shroud mass
            m_shroud = rho_m*Cx*t_shroud*2*np.pi*(rm+H/2)/N
        else:
            #Single blade shroud mass
            m_shroud = rho_m*Cx*t_shroud*2*np.pi*(rm-H/2)/N
        #Add to total
        m += m_shroud

    return m, N, m_shroud

def blade_stress(rm, H, omega, Cx, N, rho_m, A, m_shroud, sigma_y, stage_n):
    """Return the blade root stress"""

    #Find the hub and rip radii
    rt = rm+H/2
    rh = rm-H/2
    #Find the stress
    sigma_root = 0.5*rho_m*rt**2*omega**2*(1-(rh/rt)**2)+m_shroud*rt*omega**2/A
    #Warn if root stress exceeds yield
    if sigma_root > sigma_y:
        global print_warnings
        if print_warnings:
            print('WARNING: STAGE {}: MATERIAL LIMITS EXCEEDED AT BLADE ROOT'.format(stage_n))
    #Find the apparent pressure on the rotor ring from the blades
    P_blades = N*sigma_root*A/(2*np.pi*rh*1.5*Cx)

    return P_blades

def rotor_mass(rm, H_ro, H_st, Cx_ro, Cx_st, sigma_y, Po1, P_blades, omega, rho_m, nu, stage_n, n_stages):
    """Return the mass of the rotor ring"""

    global print_warnings
    #Find the outer radius of the hub
    Ro = rm-H_st/2
    #First find the rotating condition with no external pressure
    a = rho_m*omega**2*(1-nu)
    b = 2*(rho_m*omega**2*Ro**2*(1+nu)-2*sigma_y)
    c = 4*sigma_y*Ro**2-8*P_blades*Ro**2-rho_m*omega**2*(3+nu)*Ro**4
    #If the material is not strong enough then the equation breaks down. Set Ri to 0
    #and post a warning
    if b**2-4*a*c < 0:
        Ri_rot = 0
        if print_warnings:
            print('WARNING: STAGE {}: MATERIAL LIMITS EXCEEDED IN HUB: BURST'.format(stage_n))
    elif (-b-np.sqrt(b**2-4*a*c))/(2*a) < 0:
        Ri_rot = 0
        if print_warnings:
            print('WARNING: STAGE {}: MATERIAL LIMITS EXCEEDED IN HUB: BURST'.format(stage_n))
    elif stage_n == 1 or stage_n == n_stages:
        Ri_rot = 0
    else:
        Ri_rot = np.sqrt((-b-np.sqrt(b**2-4*a*c))/(2*a))
    #Find the static condition with only external pressure. Again check if material
    #limits are reached
    if Ro**2*(sigma_y-2*Po1) < 0:
        Ri_P = 0
        if print_warnings:
            print('WARNING: STAGE {}: MATERIAL LIMITS EXCEEDED HUB: PRESSURE'.format(stage_n))
    else:
        Ri_P = np.sqrt(Ro**2*(sigma_y-2*Po1)/sigma_y)
    #Take whichever inner radius is smaller
    if Ri_rot < Ri_P:
        Ri = Ri_rot
    else:
        Ri = Ri_P
    #Calculate a thickness
    thk = Ro-Ri
    #Calculate mass for spearate possibilities
    if rm-H_ro/2-thk < 0:
        #Mean hub radius of stage
        R_m = (2*rm-(H_st+H_ro)/2)/2
        #Calculate the mass of the rotor ring, which has to cover the rotor too
        m = rho_m*1.5*(Cx_ro+Cx_st)*np.pi*(R_m**2-(Ri/2)**2)
    else:
        #Mean hub radius of stage
        R_m = (2*rm-(H_st+H_ro)/2)/2-thk/2
        #Calculate the mass of the rotor ring, which has to cover the rotor too
        m = rho_m*1.5*(Cx_ro+Cx_st)*2*np.pi*R_m*thk

    return m, Ri

def stator_mass(rm, H_st, H_ro, Cx_st, Cx_ro, sigma_y, Po1, rho_m, stage_n):
    """Return the mass of the stator ring"""

    global print_warnings
    #Find the inner radius of the case
    Ri = rm+H_ro/2
    #Find the outer radius required to resist the pressure
    #Check for material failure first
    if 2*Po1 > sigma_y:
        if print_warnings:
            print('WARNING: STAGE {}: MATERIAL LIMITS EXCEEDED IN CASE: PRESSURE'.format(stage_n))
        Ro = 1.1*Ri
    else:
        Ro = np.sqrt(sigma_y*Ri**2/(sigma_y-2*Po1))
    #Calculate thickness
    thk = Ro-Ri
    #Mean case radius of stage
    R_m = (2*rm+(H_st+H_ro)/2)/2+thk/2
    #Calculate the mass of the case, which has to cover the rotor too
    m = rho_m*1.5*(Cx_st+Cx_ro)*2*np.pi*R_m*thk

    return m, Ro

def stage_mass(To1, Po1, Po3, dimensions, Omega, blade_areas, materials, stage_n, n_stages):
    """Return the total mass of the stage"""

    global print_warnings
    #Extract values from list
    t, g, r, H_st, H_ro, Cx_st, Cx_ro, w_st, w_ro, ptc_st, ptc_ro = dimensions
    A_st, A_ro = blade_areas
    #Set the materials for the rotor and stator
    stator_mat, rotor_mat = materials
    #Set safety factors
    SF_st = 1.0
    SF_ro = 1.0
    #Find the properties for the stator material:
    if stator_mat == 'inconel':
        nu_m, rho_m, strain, c_m, k_m, sigma_m, E_m = inconel(To1)
    if stator_mat == 'steel':
        nu_m, rho_m, strain, c_m, k_m, sigma_m, E_m = steel(To1)
    if sigma_m <= 0 or E_m <= 0:
        if print_warnings:
            print('WARNING: STAGE {}: MATERIAL THERMAL LIMITS REACHED'.format(stage_n))
    #Find the mass and number of stator blades
    stator_blade_sizes = [r, ptc_st, H_st, Cx_st, A_st]
    stator_blade_calc = blade_mass(stator_blade_sizes, rho_m, 'stator')
    stator_blade_mass = stator_blade_calc[0]
    stator_blade_N = stator_blade_calc[1]
    #Find the mass of the stator case and sum
    stator_case_calc = stator_mass(r, H_st, H_ro, Cx_st, Cx_ro, sigma_m/SF_st, Po1, rho_m, stage_n)
    stator_case_mass = stator_case_calc[0]
    stator_Ro = stator_case_calc[1]
    stator_total_mass = stator_case_mass+stator_blade_mass
    #Find the properties for the rotor material:
    if rotor_mat != stator_mat:
        if rotor_mat == 'inconel':
            nu_m, rho_m, strain, c_m, k_m, sigma_m, E_m = inconel(To1)
        if rotor_mat == 'steel':
            nu_m, rho_m, strain, c_m, k_m, sigma_m, E_m = steel(To1)
    if sigma_m <= 0 or E_m <= 0:
        if print_warnings:
            print('WARNING: STAGE {}: MATERIAL THERMAL LIMITS REACHED'.format(stage_n))
    #Find the mass and number of rotor blades
    rotor_blade_sizes = [r, ptc_ro, H_ro, Cx_ro, A_ro]
    rotor_blade_calc = blade_mass(rotor_blade_sizes, rho_m, 'rotor')
    rotor_blade_mass = rotor_blade_calc[0]
    rotor_blade_N = rotor_blade_calc[1]
    m_shroud = rotor_blade_calc[2]
    #Find the stress from the rotor blade roots
    P_blades = blade_stress(r, H_ro, Omega, Cx_ro, rotor_blade_N, rho_m, A_ro, m_shroud, sigma_m/SF_ro, stage_n)
    #Find the mass of the rotor ring
    rotor_ring_calc = rotor_mass(r, H_ro, H_st, Cx_ro, Cx_st, sigma_m/SF_ro, Po1, P_blades, Omega, rho_m, nu_m, stage_n, n_stages)
    rotor_ring_mass = rotor_ring_calc[0]
    rotor_Ri = rotor_ring_calc[1]
    rotor_total_mass = rotor_ring_mass+rotor_blade_mass
    #Sum for total mass and blades
    stage_mass = rotor_total_mass+stator_total_mass
    stage_blades = rotor_blade_N+stator_blade_N
    #Changes in radii from cold-static to cold-rotating
    stator_Rs = [r+H_st/2, stator_Ro]
    rotor_Rs = [rotor_Ri, r-H_ro/2]
    dr_cold = elongation(materials, stator_Rs, rotor_Rs, m_shroud, A_ro, Omega, Po1, Po3, P_blades, 293)
    #Changes in radii from cold-static to hot-rotating
    dr_hot = elongation(materials, stator_Rs, rotor_Rs, m_shroud, A_ro, Omega, Po1, Po3, P_blades, To1)

    return stage_mass, stage_blades, stator_Ro, rotor_Ri, dr_cold, dr_hot


def blade_force(P1, P2, r, H1, H2):
    """Return the axial force on the blade row"""

    #SFME with constant axial velocity: just pressure difference
    Fx = 2*np.pi*r*(H2*P2-H1*P1+(H2-H1)*(P1+P2)/2)

    return Fx

def elongation(materials, stator_Rs, rotor_Rs, m_shroud, A_ro, omega, Po1, Po3, P_blades, T):
    """Return the elongation of components"""

    #Unpack variables
    stator_mat, rotor_mat = materials
    Ri_st, Ro_st = stator_Rs
    Ri_ro, Ro_ro = rotor_Rs
    #Get material properties
    if stator_mat == 'inconel':
        nu_st, rho_st, strain_st, c_st, k_st, sigma_st, E_st = inconel(T)
    if stator_mat == 'steel':
        nu_st, rho_st, strain_st, c_st, k_st, sigma_st, E_st = steel(T)
    if rotor_mat == 'inconel':
        nu_ro, rho_ro, strain_ro, c_ro, k_ro, sigma_ro, E_ro = inconel(T)
    if rotor_mat == 'steel':
        nu_ro, rho_ro, strain_ro, c_ro, k_ro, sigma_ro, E_ro = steel(T)
    #Stator
    #Increase in inner radius due to pressure
    dr_P_st = Po1*Ri_st/E_st*((Ro_st**2+Ri_st**2)/(Ro_st**2-Ri_st**2)+nu_st)
    #Thermal expansion
    dr_T_st = strain_st*Ri_st
    #Total
    dr_st = dr_P_st+dr_T_st
    #Rotor
    #Increase in outer radius due to rotation
    dr_r_ro = rho_ro*omega**2*Ro_ro/(4*E_ro)*((1-nu_ro)*Ro_ro**2+(3+nu_ro)*Ri_ro**2)
    #Increase in outer radius due to combined blade and gas pressure
    dr_P_ro = (P_blades+Po3-Po1)*Ro_ro/E_st*((Ro_ro**2+3*Ri_ro**2)/(Ro_ro**2-Ri_ro**2)-nu_ro)
    #Increase in blade height from blade mass
    dH_blade = rho_ro*Ri_st**2*omega**2/(2*E_ro)*(2*Ri_st/3-Ro_ro+Ro_ro**3/(3*Ri_st**2))
    #Increase in blade height from shroud mass
    dH_shroud = m_shroud*Ri_st*omega**2/(E_ro*A_ro)*(Ri_st-Ro_ro)
    #Thermal expansion (outer radius being the tip of the blade)
    dr_T_ro = strain_ro*Ri_st
    #Total
    dr_ro = dr_r_ro+dr_P_ro+dr_T_ro+dH_blade+dH_shroud

    return dr_st, dr_ro

def start_up(vels, states, pressures, dimensions, materials, gas):
    """Return thermal behaviour of the stage at start up"""

    #Extract stage parameters
    V1, V2, W2, W3, V3 = vels
    T1, T2, T3, rho1, rho2, rho3, mu2, mu3, mdot = states
    P1, P2, P3 = pressures
    t, g, r, H_st, H_ro, Cx_st, Cx_ro, w_st, w_ro, ptc_st, ptc_ro, Ro_st, Ri_ro = dimensions
    stator_mat, rotor_mat = materials
    #Biot number calculation
    #First need a characteristic length s=V/A
    s_st = (Ro_st**2-(r+H_st/2)**2)/(2*(r+H_st/2))
    s_ro = ((r-H_ro/2)**2-Ri_ro**2)/(2*(r-H_ro/2))
    #Find the heat transfer coefficient for a turbulent boundary layer
    #Need average properties across the stage
    mu1 = viscosity(T1, gas)
    Re_av_st = 1.5*Cx_st*(rho1*V1/mu1+rho2*V2/mu2)/2+1.5*Cx_ro*(rho2*V2/mu2+rho3*V3/mu3)/2
    Re_av_ro = 1.5*Cx_st*(rho1*W3/mu1+rho2*W2/mu2)/2+1.5*Cx_ro*(rho2*W2/mu2+rho3*W3/mu3)/2
    Pr_av = (Prandtl(T1, P1, gas)+2*Prandtl(T2, P2, gas)+Prandtl(T3, P3, gas))/4
    k_gas_av = (conductivity(T1, gas)+2*conductivity(T2, gas)+conductivity(T3, gas))/4
    #Calculate Nusselt numbers
    Nu_st = 0.0296*Re_av_st**(4/5)*Pr_av**(1/3)
    Nu_ro = 0.0296*Re_av_ro**(4/5)*Pr_av**(1/3)
    #Calculate heat transfer coefficients
    h_st = Nu_st*k_gas_av/(1.5*(Cx_st+Cx_ro))
    h_ro = Nu_ro*k_gas_av/(1.5*(Cx_st+Cx_ro))
    #Material properties
    if stator_mat == 'inconel':
        nu_st, rho_st, strain_st, c_st, k_st, sigma_st, E_st = inconel(293)
    if stator_mat == 'steel':
        nu_st, rho_st, strain_st, c_st, k_st, sigma_st, E_st = steel(293)
    if rotor_mat == 'inconel':
        nu_ro, rho_ro, strain_ro, c_ro, k_ro, sigma_ro, E_ro = inconel(293)
    if rotor_mat == 'steel':
        nu_ro, rho_ro, strain_ro, c_ro, k_ro, sigma_ro, E_ro = steel(293)
    #Calculate the Biot numbers
    Bi_st = h_st*s_st/k_st
    Bi_ro = h_ro*s_ro/k_ro
    #Calculate lumped mass time constants. To account for different
    #conductivities, the time constant is scaled by the ratio of the rotor
    #and stator materials conductivities
    tau_st = (k_ro/k_st)*c_st*rho_st*s_st/h_st
    tau_ro = (k_st/k_ro)*c_ro*rho_ro*s_ro/h_ro

    return tau_st, tau_ro, Bi_st, Bi_ro

#############################
###INTERPOLATION FUNCTIONs###
#############################

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

#######################
###OPTIMISE FUNCTION###
#######################

def optimise(start_turbine):
    """Find an optimum design point from a given start"""

    #Extract starting inputs
    Po1, To1, mdot, Omega, W, t, g, phi, psi, Lambda, AR, dho, n, ptc, ain, gas = start_turbine[11]
    #Define starting points
    phi0 = phi[0]
    psi0 = psi[0]
    Lam0 = Lambda[0]
    AR0 = AR[0]
    dho0 = dho[0]
    phi1 = phi[-1]
    psi1 = psi[-1]
    Lam1 = Lambda[-1]
    AR1 = AR[-1]
    dho1 = dho[-1]
    #Define limits
    phi_lim = (0.1, 1.5)
    psi_lim = (0.4, 3.0)
    Lam_lim = (0, 1)
    AR_lim = (0.1, 5)
    dh_lim = (1, 5)

    def turbine_calcs(args):
        """Takes a list of arguments and feeds them to the turbine function"""
        #Extract arguments
        ph1, ph2, ps1, ps2, L1, L2, AR1, AR2, dh1, dh2 = [i for i in args]
        #Form them into lists
        phi = [ph1, ph2]
        psi = [ps1, ps2]
        Lambda = [L1, L2]
        dho = [dh1, dh2]
        AR = [AR1, AR2]
        #Pass them to the turbine to find efficiency
        eff = turbine(Po1, To1, mdot, Omega, W, t, g, phi, psi, Lambda, AR, dho, n, ptc, ain, gas)[0]
        #Return the negative of efficiency to allow for minimising
        return -eff

    def constraint_a2(args):
        """Constraint on a2"""
        #Extract and form lists
        ph1, ph2, ps1, ps2, L1, L2, AR1, AR2, dh1, dh2 = [i for i in args]
        phi = [ph1, ph2]
        psi = [ps1, ps2]
        Lambda = [L1, L2]
        #Spline over n
        phi = spline(n, phi)
        psi = spline(n, psi)
        Lambda = spline(n, Lambda)
        #Calculate angles
        a2_max = 0
        a1 = ain
        for i in range(len(phi)):
            ang_check = angles(phi[i], psi[i], Lambda[i], a1)
            a2 = np.degrees(ang_check[1])
            a1 = ang_check[3]
            if abs(a2) > a2_max:
                a2_max = abs(a2)
        #Return a value that should be greater than zero
        return 73-a2_max

    def constraint_b3(args):
        """Constraint on b2"""
        #Extract and form lists
        ph1, ph2, ps1, ps2, L1, L2, AR1, AR2, dh1, dh2 = [i for i in args]
        phi = [ph1, ph2]
        psi = [ps1, ps2]
        Lambda = [L1, L2]
        #Spline over n
        phi = spline(n, phi)
        psi = spline(n, psi)
        Lambda = spline(n, Lambda)
        #Calculate angles
        b3_max = 0
        a1 = ain
        for i in range(len(phi)):
            ang_check = angles(phi[i], psi[i], Lambda[i], a1)
            b3 = np.degrees(ang_check[4])
            a1 = ang_check[3]
            if  abs(b3) > b3_max:
                b3_max = abs(b3)
        #Return a value that should be greater than zero
        return 73-b3_max
    #Form the starting point list
    x0 = [phi0, phi1, psi0, psi1, Lam0, Lam1, AR0, AR1, dho0, dho1]
    #Form the tuple of bounds
    bnds = (phi_lim, phi_lim, psi_lim, psi_lim, Lam_lim, Lam_lim, AR_lim, AR_lim, dh_lim, dh_lim)
    #Form the tuple of constraints
    cons = cons = ({'type': 'ineq', 'fun': constraint_a2}, {'type': 'ineq', 'fun': constraint_b3})
    #Find the minimum
    res = minimize(turbine_calcs, x0, method='SLSQP', bounds=bnds, constraints=cons)
    #Extract the optimal variable and return them
    phi1, phi2, psi1, psi2, Lam1, Lam2, AR1, AR2, dho1, dho2 = res['x']
    return [phi1, phi2], [psi1, psi2], [Lam1, Lam2], [AR1, AR2], [dho1, dho2]
