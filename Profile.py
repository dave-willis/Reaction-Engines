"""Module contains functions generating blade profiles from CJC and plotting
functions from WD"""

#Import required modules
import math
import numpy as np
import scipy.optimize as SciOpt
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, TextBox
from Stage import turbine

#################################
###PROFILE GENERATOR FUNCTIONS###
#################################

def binomial(i, n):
    """Binomial coefficient"""

    return math.factorial(n) / float(math.factorial(i) * math.factorial(n - i))

def bernstein(t, i, n):
    """Bernstein polynom"""

    return binomial(i, n) * (t ** i) * ((1 - t) ** (n - i))

def bezier(t, points):
    """Calculate coordinate of a point in the bezier curve"""

    n = len(points) - 1
    x = y = 0
    for i, pos in enumerate(points):
        bern = bernstein(t, i, n)
        x += pos[0] * bern
        y += pos[1] * bern

    return x, y

def bezier_curve_range(n, points):
    """Range of points in a curve bezier"""

    for i in range(n):
        t = i / float(n - 1)
        yield bezier(t, points)

def GEN_TURNING(Xi1, Xi2, controls):
    """Generate x,y of a bezier from control points"""

    spacings = np.linspace(0, 1, len(controls)+2)[1:-1]
    CPs = []
    CPs.append([0.0, 0.0])
    for i in range(len(controls)):
        CPs.append([spacings[i]+min(spacings[i], 1.-spacings[i])*controls[i], spacings[i]-min(spacings[i], 1.-spacings[i])*controls[i]])
    CPs.append([1.0, 1.0])

    steps = 100
    x = []
    y = []
    for point in bezier_curve_range(steps, CPs):
        x.append(point[0])
        y.append(point[1])

    ### integrate profile
    xcam = [0.]
    ycam = [0.]

    DXi = Xi2-Xi1
    for i in range(1, len(x)):
        Xi = Xi1+DXi*(y[i]+y[i-1])*0.5
        xcam.append(xcam[-1]+(x[i]-x[i-1])*np.cos(Xi))
        ycam.append(ycam[-1]+(x[i]-x[i-1])*np.sin(Xi))

    DX = xcam[-1]
    ycam = ycam/DX
    xcam = xcam/DX

    return(np.array(xcam), np.array(ycam))

def GEN_BZ_2D(controls_x, controls_y):
    """Generate x,y of a bezier from control points"""

    CPs = []
    for i in range(len(controls_x)):
        CPs.append([controls_x[i], controls_y[i]])

    steps = 100
    controlPoints = CPs
    x = []
    y = []
    for point in bezier_curve_range(steps, controlPoints):
        x.append(point[0])
        y.append(point[1])

    return(np.array(x), np.array(y))

def dist_vino(n, xst, xen, s1, s2):
    """Vinokur distribution"""

    ### INPUT DESCRIPTION ###
    ### n - number of points
    ### xst - xstart
    ### xen - xend
    ### s1 - start spacing
    ### s2 - end spacing
    #####

    s1 = s1/(xen - xst)
    s2 = s2/(xen - xst)
    eta = np.arange(n)
    a = math.sqrt(s2)/math.sqrt(s1)
    b = 1.0/((n-1)*math.sqrt(s1*s2))
    nm1 = float(n-1)
    trans = lambda delta: math.sinh(delta)/delta - b
    left = 0.00001
    delta = SciOpt.brentq(trans, left, 100)
    u = 0.5*(1.0 + np.tanh(delta*(eta/nm1 - 0.5))/(np.tanh(delta*0.5)))
    s = u/(a + (1-a)*u)

    return xst + (xen-xst)*s


def F_CF(psi):
    """Function to generate class function of psi"""

    C = (psi**0.5)*(1.0-psi)
    return C

def F_TF_bez(controls_cam_x, controls_cam_y, TE, Oval_frac):
    """Function that returns thickness contributions due to parameters"""

    psi = dist_vino(200, 0, 1, 1./2000, 1./500)
    X, S = GEN_BZ_2D(controls_cam_x, controls_cam_y)
    S = np.interp(psi, X, S)

    S = S+Oval_frac*S[0]*(1.-psi)**18
    C = F_CF(psi)		### calculate class function

    TF = S*C+psi*TE/2.	### combine to form thickness funtion

    return TF


def calc_norm(x_cam, y_cam):
    """Calculate the normal to the camber line"""

    s = np.zeros((len(x_cam)))		### initialise camber length array
    grad = np.zeros((len(x_cam), 2))	### initialise camber gradient array
    norm = np.zeros((len(x_cam), 2))	### initialise camber normals array

    for i in range(1, len(x_cam)):
        s[i] = s[i-1]+((x_cam[i]-x_cam[i-1])**2+(y_cam[i]-y_cam[i-1])**2)**0.5	### calculate cumsum length of camber

    grad[:-1, 0] = (x_cam[1:]-x_cam[:-1])/(s[1:]-s[:-1])	### calculate dx/ds
    grad[:-1, 1] = (y_cam[1:]-y_cam[:-1])/(s[1:]-s[:-1])	### calculate dy/ds
    grad[-1, :] = grad[-2, :]							### deal with final point
    norm[:, 0] = -grad[:, 1]							### calculate normal 1 as -1/grad2
    norm[:, 1] = grad[:, 0]

    return norm

def F_Make(Xi1, Xi2, controls_cam, controls_thk_x, controls_thk_y, Tte, Oval_frac):
    """Function that compiles a blade from reduced variables"""

    ### turning distribution
    Xi1 = np.radians(Xi1)
    Xi2 = np.radians(Xi2)
    (x_cam, y_cam) = GEN_TURNING(Xi1, Xi2, controls_cam)      ### generate camber line from bezier curve

    s2 = np.zeros((len(x_cam)))			### initialise camber length array
    for i in range(1, len(x_cam)):
        s2[i] = s2[i-1]+((x_cam[i]-x_cam[i-1])**2+(y_cam[i]-y_cam[i-1])**2)**0.5
    psi = s2/s2.max()

    psi_new = dist_vino(200, 0, 1, 1./2000, 1./500)

    x_cam = np.interp(psi_new, psi, x_cam)
    y_cam = np.interp(psi_new, psi, y_cam)
    psi = psi_new
    norm = calc_norm(x_cam, y_cam)

    thk = F_TF_bez(controls_thk_x, controls_thk_y, Tte, Oval_frac)

    Z_U = y_cam+thk*norm[:, 1]	### apply thickness onto camber line for upper
    Z_L = y_cam-thk*norm[:, 1]	### apply thickness onto camber line for lower
    X_U = x_cam+thk*norm[:, 0]	### repeat upper for x rather than y
    X_L = x_cam-thk*norm[:, 0]

    return(X_L, Z_L, X_U, Z_U)


def Profile(X1, X2, TKTE, Cx, points=500):
    """Function based on profgen.f (JDD) - GENERATES SIMPLE BLADE PROFILES"""

    ### THIS FUNCTION SHOULD BE REPLACED AT SOME POINT TO USE A BETTER BLADE SECTION GENERATOR
    ### EXCUSE THE WEIRD FORMULATION, THIS IS CONVERSION FROM FORTRAN
    ### USER HARD CODED VALUES BELOW - COULD BE TAKEN OUTSIDE FUNCTION CJC

    controls_cam = [-0.0, -0.5, -0.8]
    Rle = 0.05
    Tte = TKTE
    Beta_te = 4.0
    Oval_frac = 0.3
    #controls_cam = [0,0.2,0.0]
    controls_thk_x = [0.0, 0.5, 1.0]
    controls_thk_y = [(2.*Rle)**0.5*(1.-Oval_frac), 0.25, np.tan(np.radians(Beta_te))+Tte/2.]

    (XSIN, YSIN, XPIN, YPIN) = F_Make(X1, X2, controls_cam, controls_thk_x, controls_thk_y, Tte, Oval_frac)

    XIN = (XPIN+XSIN)/2.
    YIN = (YPIN+YSIN)/2.
    XScale = max(max(XPIN), max(XSIN))
    XPIN = XPIN/XScale
    XSIN = XSIN/XScale
    XIN = XIN/XScale

    yoffset = np.mean(YIN)
    YPIN = (YPIN-yoffset)/XScale    ### center on peak thickness - DEFAULTS STACKING TO PEAK THICKNESS
    YSIN = (YSIN-yoffset)/XScale    ### center on peak thickness - DEFAULTS STACKING TO PEAK THICKNESS
    YIN = (YIN - yoffset)/XScale    ### center on peak thickness - DEFAULTS STACKING TO PEAK THICKNESS

    Xnew = np.zeros((len(XPIN)+len(XSIN)-1))
    Xnew[:200] = XPIN[::-1]
    Xnew[199:] = XSIN[:]

    Ynew = np.zeros((len(YPIN)+len(YSIN)-1))
    Ynew[:200] = YPIN[::-1]
    Ynew[199:] = YSIN[:]

    Snew = np.zeros((len(Ynew)))
    for i in range(1, len(Xnew)):
        Snew[i] = Snew[i-1]+((Xnew[i]-Xnew[i-1])**2.+(Ynew[i]-Ynew[i-1])**2.)**0.5

    for i in range(1, len(Xnew)):
        if Xnew[i] == Xnew.min():
            Sle = Snew[i]
            Xle = Xnew[i]
            break

    XPIN = np.interp(np.linspace(0, Sle, points), Snew, Xnew)[::-1]
    YPIN = np.interp(np.linspace(0, Sle, points), Snew, Ynew)[::-1]
    XSIN = np.interp(np.linspace(Sle, Snew[-1], points), Snew, Xnew)
    YSIN = np.interp(np.linspace(Sle, Snew[-1], points), Snew, Ynew)

    XScale = max(max(XPIN), max(XSIN))-Xle
    XPIN = (XPIN-Xle)/XScale
    XSIN = (XSIN-Xle)/XScale
    XIN = (XIN-Xle)/XScale

    YPIN = (YPIN-yoffset)/XScale	### center on peak thickness - DEFAULTS STACKING TO PEAK THICKNESS
    YSIN = (YSIN-yoffset)/XScale	### center on peak thickness - DEFAULTS STACKING TO PEAK THICKNESS
    YIN = (YIN - yoffset)/XScale ### center on peak thickness - DEFAULTS STACKING TO PEAK THICKNESS

    TE_circ_r = 0.5*((abs(XSIN[-1]-XPIN[-1]))**2+(abs(YSIN[-1]-YPIN[-1]))**2)**0.5
    TE_circ_cx = 0.5*(XSIN[-1]+XPIN[-1])
    TE_circ_cy = 0.5*(YSIN[-1]+YPIN[-1])
    TEcx = np.linspace(XSIN[-1], TE_circ_cx+TE_circ_r-TE_circ_r/(points*100), points)
    TEcx2 = np.linspace(XPIN[-1], TE_circ_cx+TE_circ_r-TE_circ_r/(points*100), points)
    TEcy = np.asarray([TE_circ_cy-(TE_circ_r**2-(i-TE_circ_cx)**2)**0.5 for i in TEcx])
    TEcy2 = np.asarray([TE_circ_cy+(TE_circ_r**2-(i-TE_circ_cx)**2)**0.5 for i in TEcx2])

    maxes = np.array([np.amax(XIN), np.amax(XPIN), np.amax(XSIN), np.amax(TEcx), np.amax(TEcx2)])
    mins = np.array([np.amin(XIN), np.amin(XPIN), np.amin(XSIN)])
    cx = Cx*(np.amax(maxes)-np.amin(mins))

    XP = [cx*i for i in XPIN]
    YP = [cx*i for i in YPIN]
    XS = [cx*i for i in XSIN]
    YS = [cx*i for i in YSIN]
    TEcx = [cx*i for i in TEcx]
    TEcy = [cx*i for i in TEcy]
    TEcx2 = [cx*i for i in TEcx2]
    TEcy2 = [cx*i for i in TEcy2]

    return(XP, YP, XS, YS, TEcx, TEcy, TEcx2, TEcy2)

#######################
###GEOMETRY PLOTTING###
#######################

def b2b_data(turbine_data):
    """Return blade-to-blade profiles for the whole turbine"""

    #Extract values from the turbine function output
    angs = [[i[0], i[1], i[2], i[4]] for i in turbine_data[10]]
    chords = [[i[3], i[4]] for i in turbine_data[5]]
    ptcs = [[i[5], i[6]] for i in turbine_data[5]]
    t = turbine_data[11][11]
    #Number of points per blade plot
    points = 500
    #Scale the outputs for m/cm/mm etc
    scale = 1000
    #Initialise values
    l = 0
    xp = np.empty([4*len(angs), points])
    yp = np.empty([4*len(angs), points])
    xs = np.empty([4*len(angs), points])
    ys = np.empty([4*len(angs), points])
    tecx = np.empty([4*len(angs), points])
    tecy = np.empty([4*len(angs), points])
    tecx2 = np.empty([4*len(angs), points])
    tecy2 = np.empty([4*len(angs), points])
    #Loop over every stage
    for i in range(len(angs)):
        #Extract stage parameters
        a1, a2, b2, b3 = [j for j in angs[i]]
        Cxst, Cxro = [j for j in chords[i]]
        ptcst, ptcro = [j for j in ptcs[i]]
        #Determine start point depending on b2b spacing
        if i == 0:
            x0st = 0
            x0ro = Cxst*1.25+Cxro*0.25
        else:
            x0st = l+Cxst*0.25
            x0ro = l+Cxst*1.5+Cxro*0.25
        #Calculate the length along the turbine
        if i == 0:
            l += Cxst*1.25+Cxro*1.5
        else:
            l += Cxst*1.5+Cxro*1.5
        #Pass parameters to profile function for stator
        XP, YP, XS, YS, TEcx, TEcy, TEcx2, TEcy2 = Profile(a1, a2, t/Cxst, Cxst, points)
        #Store results
        xp[4*i] = [j+x0st for j in XP]
        yp[4*i] = YP
        xs[4*i] = [j+x0st for j in XS]
        ys[4*i] = YS
        tecx[4*i] = [j+x0st for j in TEcx]
        tecy[4*i] = TEcy
        tecx2[4*i] = [j+x0st for j in TEcx2]
        tecy2[4*i] = TEcy2
        xp[4*i+1] = [j+x0st for j in XP]
        yp[4*i+1] = [j+Cxst*ptcst for j in YP]
        xs[4*i+1] = [j+x0st for j in XS]
        ys[4*i+1] = [j+Cxst*ptcst for j in YS]
        tecx[4*i+1] = [j+x0st for j in TEcx]
        tecy[4*i+1] = [j+Cxst*ptcst for j in TEcy]
        tecx2[4*i+1] = [j+x0st for j in TEcx2]
        tecy2[4*i+1] = [j+Cxst*ptcst for j in TEcy2]
        #Pass parameters to profile function for rotor
        XP, YP, XS, YS, TEcx, TEcy, TEcx2, TEcy2 = Profile(b2, b3, t/Cxro, Cxro, points)
        #Store results
        xp[4*i+2] = [j+x0ro for j in XP]
        yp[4*i+2] = YP
        xs[4*i+2] = [j+x0ro for j in XS]
        ys[4*i+2] = YS
        tecx[4*i+2] = [j+x0ro for j in TEcx]
        tecy[4*i+2] = TEcy
        tecx2[4*i+2] = [j+x0ro for j in TEcx2]
        tecy2[4*i+2] = TEcy2
        xp[4*i+3] = [j+x0ro for j in XP]
        yp[4*i+3] = [j+Cxro*ptcro for j in YP]
        xs[4*i+3] = [j+x0ro for j in XS]
        ys[4*i+3] = [j+Cxro*ptcro for j in YS]
        tecx[4*i+3] = [j+x0ro for j in TEcx]
        tecy[4*i+3] = [j+Cxro*ptcro for j in TEcy]
        tecx2[4*i+3] = [j+x0ro for j in TEcx2]
        tecy2[4*i+3] = [j+Cxro*ptcro for j in TEcy2]

    return scale*xp, scale*yp, scale*xs, scale*ys, scale*tecx, scale*tecy, scale*tecx2, scale*tecy2

def b2b_plot(turbine_data):
    """Plot blade-to-blade profiles for the whole turbine"""

    #Use data from turbine function to get data for profiles
    data = b2b_data(turbine_data)
    #Plot each bit of data, which are in the form of 2D arrays
    plt.figure()
    for i in range(4):
        plt.plot(data[2*i].T, data[2*i+1].T, 'black', lw=2)
    plt.axis('equal')
    plt.xlabel('Distance along turbine (mm)')
    plt.ylabel('Tangential distance (mm)')
    plt.show()

def b2b_variable(turbine_data):
    """Plot blade-to-blade profiles with variable loadings"""

    #Extract parameters from turbine data
    Po1, To1, mdot, Omega, phi, psi, Lambda, AR, W, dho, n, t, g, ptoC, a1i = turbine_data[11]
    phi01 = phi[0]
    phi02 = phi[-1]
    psi01 = psi[0]
    psi02 = psi[-1]
    Lambda01 = Lambda[0]
    Lambda02 = Lambda[-1]
    dho01 = dho[0]
    dho02 = dho[-1]
    AR01 = AR[0]
    AR02 = AR[-1]
    ptc01 = ptoC[0]
    ptc02 = ptoC[-1]
    #Use function to get data for b2b plot
    data = b2b_data(turbine_data)
    #Initialise plots
    fig, ax = plt.subplots()
    plt.subplots_adjust(bottom=0.25)
    pres = plt.plot(data[0].T, data[1].T, 'black', lw=2)
    suc = plt.plot(data[2].T, data[3].T, 'black', lw=2)
    te1 = plt.plot(data[4].T, data[5].T, 'black', lw=2)
    te2 = plt.plot(data[6].T, data[7].T, 'black', lw=2)
    plt.axis('equal')
    plt.xlabel('Distance along turbine (mm)')
    plt.ylabel('Tangential distance (mm)')
    #Create axes for position of sliders [left, bottom, height, width]
    axphi1 = plt.axes([0.12, 0.16, 0.3, 0.02])
    axphi2 = plt.axes([0.6, 0.16, 0.3, 0.02])
    axpsi1 = plt.axes([0.12, 0.13, 0.3, 0.02])
    axpsi2 = plt.axes([0.6, 0.13, 0.3, 0.02])
    axLambda1 = plt.axes([0.12, 0.1, 0.3, 0.02])
    axLambda2 = plt.axes([0.6, 0.1, 0.3, 0.02])
    axdho1 = plt.axes([0.12, 0.07, 0.3, 0.02])
    axdho2 = plt.axes([0.6, 0.07, 0.3, 0.02])
    axAR1 = plt.axes([0.12, 0.04, 0.3, 0.02])
    axAR2 = plt.axes([0.6, 0.04, 0.3, 0.02])
    axptc1 = plt.axes([0.12, 0.01, 0.3, 0.02])
    axptc2 = plt.axes([0.6, 0.01, 0.3, 0.02])
    #Create sliders for each variable
    sphi1 = Slider(axphi1, '$\phi_1$', 0.1, 1.0, valinit=phi01)
    sphi2 = Slider(axphi2, '$\phi_2$', 0.1, 1.0, valinit=phi02)
    spsi1 = Slider(axpsi1, '$\psi_1$', 0.5, 2.5, valinit=psi01)
    spsi2 = Slider(axpsi2, '$\psi_2$', 0.5, 2.5, valinit=psi02)
    sLambda1 = Slider(axLambda1, '$\Lambda_1$', 0.01, 0.99, valinit=Lambda01)
    sLambda2 = Slider(axLambda2, '$\Lambda_2$', 0.01, 0.99, valinit=Lambda02)
    sdho1 = Slider(axdho1, '$\Delta h_{o1}$', 1, 3, valinit=dho01)
    sdho2 = Slider(axdho2, '$\Delta h_{o2}$', 1, 3, valinit=dho02)
    sAR1 = Slider(axAR1, '$AR_1$', 0.9, 2.0, valinit=AR01)
    sAR2 = Slider(axAR2, '$AR_2$', 0.9, 2.0, valinit=AR02)
    sptc1 = Slider(axptc1, '$p/C_1$', 0.5, 1.5, valinit=ptc01)
    sptc2 = Slider(axptc2, '$p/C_2$', 0.5, 1.5, valinit=ptc02)
    #Text box showing turbine efficiency
    effax = plt.axes([0.85, 0.9, 0.09, 0.04])
    eff = TextBox(effax, '', 'Efficiency: {}%'.format(np.round(100*turbine_data[0], 2)), color='1.0')
    #Text box showing maximum angle
    angle_maxax = plt.axes([0.72, 0.9, 0.11, 0.04])
    angle_max = TextBox(angle_maxax, '', 'Maximum angle: {}º'.format(np.round(turbine_data[12][1], 2)), color='g', hovercolor='g')
    #Update the graphs and buttons when the sliders are changed
    def update(val):
        """Update the plots"""
        #Update the inputs with the slider values
        phi = [sphi1.val, sphi2.val]
        psi = [spsi1.val, spsi2.val]
        Lambda = [sLambda1.val, sLambda2.val]
        dho = [sdho1.val, sdho2.val]
        AR = [sAR1.val, sAR2.val]
        ptoC = [sptc1.val, sptc2.val]
        #Recalculate turbine performance
        new_turbine = turbine(Po1, To1, mdot, Omega, phi, psi, Lambda, AR, W, dho, n, t, g, ptoC, a1i)
        #Get new profile data
        data = b2b_data(new_turbine)
        #Update plots
        for i in range(n):
            for j in range(4):
                pres[4*i+j].set_xdata(data[0].T[:, 4*i+j])
                pres[4*i+j].set_ydata(data[1].T[:, 4*i+j])
                suc[4*i+j].set_xdata(data[2].T[:, 4*i+j])
                suc[4*i+j].set_ydata(data[3].T[:, 4*i+j])
                te1[4*i+j].set_xdata(data[4].T[:, 4*i+j])
                te1[4*i+j].set_ydata(data[5].T[:, 4*i+j])
                te2[4*i+j].set_xdata(data[6].T[:, 4*i+j])
                te2[4*i+j].set_ydata(data[7].T[:, 4*i+j])
        #Update the efficiency and angle boxes
        eff.set_val('Efficiency: {}%'.format(np.round(100*new_turbine[0], 2)))
        angle_max.set_val('Maximum angle: {}º'.format(np.round(new_turbine[12][1], 2)))
        #Change the colour of the angle box if needed
        if new_turbine[12][0]:
            angle_max.color = 'r'
            angle_max.hovercolor = 'r'
        else:
            angle_max.color = 'g'
            angle_max.hovercolor = 'g'
        #Update the figure
        ax.set_xbound
        fig.canvas.draw_idle()
    #When any of the sliders are changed, update the graph
    sphi1.on_changed(update)
    sphi2.on_changed(update)
    spsi1.on_changed(update)
    spsi2.on_changed(update)
    sLambda1.on_changed(update)
    sLambda2.on_changed(update)
    sdho1.on_changed(update)
    sdho2.on_changed(update)
    sAR1.on_changed(update)
    sAR2.on_changed(update)
    sptc1.on_changed(update)
    sptc2.on_changed(update)
    plt.axis('equal')
    #Reset button to return sliders to initial values
    resetax = plt.axes([0.47, 0.05, 0.08, 0.04])
    reset_button = Button(resetax, 'Reset', color='1.0', hovercolor='0.5')
    def reset(event):
        """Reset the sliders"""
        sphi1.reset()
        sphi2.reset()
        spsi1.reset()
        spsi2.reset()
        sLambda1.reset()
        sLambda2.reset()
        sdho1.reset()
        sdho2.reset()
        sAR1.reset()
        sAR2.reset()
        sptc1.reset()
        sptc2.reset()
    reset_button.on_clicked(reset)
    reset_button.on_clicked(update)
    #Button that sets the exit loading equal to inlet to create repeating stages
    repeatingax = plt.axes([0.47, 0.1, 0.08, 0.04])
    repeating_button = Button(repeatingax, 'Repeating', color='1.0', hovercolor='0.5')
    def repeat(event):
        """Set output loading to input"""
        sphi2.set_val(sphi1.val)
        spsi2.set_val(spsi1.val)
        sLambda2.set_val(sLambda1.val)
        sdho2.set_val(sdho1.val)
        sAR2.set_val(sAR1.val)
        sptc2.set_val(sptc1.val)
    repeating_button.on_clicked(repeat)
    repeating_button.on_clicked(update)
    #These create dummy variables to ensure buttons can be referenced outside of function
    repeatingax._button = repeating_button
    resetax._button = reset_button

def annulus(turbine_data, close=True, scale=True):
    """Plot annulus size for stators through the turbine, return the plotting points"""

    #Extract values form turbine function output
    chords = [[i[3], i[4]] for i in turbine_data[5]]
    Hst = [i[1] for i in turbine_data[5]]
    rm = [i[0] for i in turbine_data[5]]
    #Initialise lists
    l = 0
    x = []
    r_hub = []
    r_cas = []
    #Loop over every stage
    for i in range(len(rm)):
        #Extract stage parameters
        Cxst, Cxro = [j for j in chords[i]]
        #Calculate the length along the turbine
        if i == 0:
            l += Cxst*1.25+Cxro*1.5
        l += Cxst*1.5+Cxro*1.5
        x.append(l)
        #Calculate the hub and case radii
        r_hub.append(rm[i]-Hst[i]/2)
        r_cas.append(rm[i]+Hst[i]/2)
    #Plot the results
    if close:
        plt.figure()
        plt.plot(x, r_hub, 'black')
        plt.plot(x, r_cas, 'black')
        plt.xlabel('Length along turbine (m)')
        plt.ylabel('Radius (m)')
        plt.show()
    if scale:
        plt.figure()
        plt.plot(x, r_hub, 'black', linewidth=0.5)
        plt.plot(x, r_cas, 'black', linewidth=0.5)
        plt.xlabel('Length along turbine (m)')
        plt.ylabel('Radius (m)')
        plt.ylim(0, 1.1*np.amax(r_cas))
        plt.show()

    return x, r_hub, r_cas
