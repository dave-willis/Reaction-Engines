#Import required modules
import math
import numpy as np
import scipy.optimize as SciOpt
import matplotlib.pyplot as plt

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

    return(X_L, Z_L, X_U, Z_U, thk)


def Profile(X1, X2, TKTE, points=500, TE_points=200):
    """Function based on profgen.f (JDD) - GENERATES SIMPLE BLADE PROFILES"""

    ### THIS FUNCTION SHOULD BE REPLACED AT SOME POINT TO USE A BETTER BLADE SECTION GENERATOR
    ### EXCUSE THE WEIRD FORMULATION, THIS IS CONVERSION FROM FORTRAN
    ### USER HARD CODED VALUES BELOW - COULD BE TAKEN OUTSIDE FUNCTION CJC
    ### NOTE: XP/XS ETC DON'T ALWAYS REFER TO THE PRESSURE/SUCTION SURFACEC=S
    ### 'S' VARIABLES ARE THE LOWER SURFACE AND 'P' THE UPPER, SO IF EXIT ANGLE
    ### IS GREATER THAN INLET THIS WILL BE THE SUCTION AND PRESSURE SURFACES

    controls_cam = [-0.0, -0.5, -0.8]
    Rle = 0.05
    Tte = TKTE
    Beta_te = 4.0
    Oval_frac = 0.3
    controls_thk_x = [0.0, 0.5, 1.0]
    controls_thk_y = [(2*Rle)**0.5*(1-Oval_frac), 0.25, np.tan(np.radians(Beta_te))+Tte/2]

    (XLIN, YLIN, XUIN, YUIN, thk) = F_Make(X1, X2, controls_cam, controls_thk_x, controls_thk_y, Tte, Oval_frac)

    XIN = (XUIN+XLIN)/2
    YIN = (YUIN+YLIN)/2
    XScale = max(max(XUIN), max(XLIN))
    XUIN = XUIN/XScale
    XLIN = XLIN/XScale
    XIN = XIN/XScale

    yoffset = np.mean(YIN)
    YUIN = (YUIN-yoffset)/XScale    ### center on peak thickness - DEFAULTS STACKING TO PEAK THICKNESS
    YLIN = (YLIN-yoffset)/XScale    ### center on peak thickness - DEFAULTS STACKING TO PEAK THICKNESS
    YIN = (YIN - yoffset)/XScale    ### center on peak thickness - DEFAULTS STACKING TO PEAK THICKNESS

    Xnew = np.zeros((len(XUIN)+len(XLIN)-1))
    Xnew[:200] = XUIN[::-1]
    Xnew[199:] = XLIN[:]

    Ynew = np.zeros((len(YUIN)+len(YLIN)-1))
    Ynew[:200] = YUIN[::-1]
    Ynew[199:] = YLIN[:]

    Snew = np.zeros((len(Ynew)))
    for i in range(1, len(Xnew)):
        Snew[i] = Snew[i-1]+((Xnew[i]-Xnew[i-1])**2.+(Ynew[i]-Ynew[i-1])**2.)**0.5

    for i in range(1, len(Xnew)):
        if Xnew[i] == Xnew.min():
            Sle = Snew[i]
            Xle = Xnew[i]
            break

    XUIN = np.interp(np.linspace(0, Sle, points), Snew, Xnew)[::-1]
    YUIN = np.interp(np.linspace(0, Sle, points), Snew, Ynew)[::-1]
    XLIN = np.interp(np.linspace(Sle, Snew[-1], points), Snew, Xnew)
    YLIN = np.interp(np.linspace(Sle, Snew[-1], points), Snew, Ynew)

    XScale = max(max(XUIN), max(XLIN))-Xle
    XUIN = (XUIN-Xle)/XScale
    XLIN = (XLIN-Xle)/XScale
    XIN = (XIN-Xle)/XScale

    YUIN = (YUIN-yoffset)/XScale	### center on peak thickness - DEFAULTS STACKING TO PEAK THICKNESS
    YLIN = (YLIN-yoffset)/XScale	### center on peak thickness - DEFAULTS STACKING TO PEAK THICKNESS
    YIN = (YIN - yoffset)/XScale ### center on peak thickness - DEFAULTS STACKING TO PEAK THICKNESS

    UTE_m = (XUIN[-2]-XUIN[-1])/(YUIN[-1]-YUIN[-2])
    LTE_m = (XLIN[-2]-XLIN[-1])/(YLIN[-1]-YLIN[-2])
    TE_circ_cx = (YUIN[-1]-YLIN[-1]+LTE_m*XLIN[-1]-UTE_m*XUIN[-1])/(LTE_m-UTE_m)
    TE_circ_cy = YUIN[-1]+UTE_m*(TE_circ_cx-XUIN[-1])
    TE_circ_r = ((XLIN[-1]-TE_circ_cx)**2+(YLIN[-1]-TE_circ_cy)**2)**0.5
    U_theta = np.arccos((XUIN[-1]-TE_circ_cx)/TE_circ_r)
    L_theta = -np.arccos((XLIN[-1]-TE_circ_cx)/TE_circ_r)
    
    TE_points = 200
    TE_theta = np.linspace(0.95*U_theta, L_theta, TE_points)
    TEx = TE_circ_cx - TE_circ_r*np.cos(np.pi-TE_theta)
    TEy = TE_circ_cy + TE_circ_r*np.sin(np.pi-TE_theta)

    maxes = np.array([np.amax(XIN), np.amax(XUIN), np.amax(XLIN), np.amax(TEx)])
    mins = np.array([np.amin(XIN), np.amin(XUIN), np.amin(XLIN)])
    cx = np.amax(maxes)-np.amin(mins)

    XU = [cx*i for i in XUIN]
    YU = [cx*i for i in YUIN]
    XL = [cx*i for i in XLIN]
    YL = [cx*i for i in YLIN]
    TEx = [cx*i for i in TEx]
    TEy = [cx*i for i in TEy]
    
    te_index = int(points+TE_points/2)
    X = XU + TEx + XL[::-1]
    Y = YU + TEy + YL[::-1]

    XU = X[:te_index]
    YU = Y[:te_index]
    XL = X[te_index:][::-1]
    YL = Y[te_index:][::-1]
    TE = [TE_circ_cx, TE_circ_cy, TE_circ_r, TEx, TEy]

    return(XU, YU, XL, YL, XIN, YIN, X, Y, TE, thk)

X1 = 40
X2 = -60
TKTE = 0.04
points = 500
TE_points=100

XU, YU, XL, YL, XIN, YIN, X, Y, TE, thk = Profile(X1, X2, TKTE, points, TE_points)
TE_circ_cx, TE_circ_cy, TE_circ_r, TEcx, TEcy = TE

cx = np.amax(X)-np.amin(X)
ptc = 0.7

if 1 == 1: #Blade geometry labels 
    
    i=25 #Incidence angle
    d=25 #Deviation angle
    arrows=0.015*cx
    n = 100
    lw = 0.3 #Plot line width
    dash = (0, (5, 30))
    
    plt.plot(X,Y,'black',linewidth=2)
    plt.plot(X,[i+ptc*cx for i in Y],'black',linewidth=2)
    plt.plot(XIN,[i+ptc*cx for i in YIN],'black',linestyle='--', linewidth=1)
    plt.axis('off')
    plt.margins(0.05,0.05)
    plt.axis('equal')
    
    ##### THROAT LINES #####
    if 1 == 1:
        dist_min = 100
        throat = 0 #Index of throat
        for j in range(len(YU)):
            if np.sqrt((XL[points]-XU[j])**2+(YL[points]+ptc*cx-YU[j])**2)<dist_min:
                dist_min = np.sqrt((XL[points]-XU[j])**2+(YL[points]+ptc*cx-YU[j])**2)
                throat = j
        
        ##### THROAT #####
        if 1 == 1:
            throat_dx = XL[points]-XU[throat]
            throat_dy = YL[points]+ptc*cx-YU[throat]
            guide1x = np.linspace(XU[throat],XL[points],n)
            guide1y = [YL[points]+ptc*cx+(throat_dy/throat_dx)*(q-XL[points]) for q in guide1x]
            plt.plot(guide1x,guide1y,'black',linewidth=0.1)
        
        ##### UNCOVERED TURNING ANGLE #####
        if 2 == 1:
            throat_m = (YU[throat]-YU[throat+1])/(XU[throat]-XU[throat+1])
            uncov_max_x = XU[-1]+(YU[throat]+throat_m*(XU[-1]-XU[throat])-YU[-1])*np.cos(abs(np.arctan(throat_m)))*np.sin(abs(np.arctan(throat_m)))
            uncov_x = np.linspace(XU[throat],uncov_max_x, n)
            plt.plot(uncov_x,[YU[throat]+throat_m*(i-XU[throat]) for i in uncov_x],'black',linestyle=dash,linewidth=lw)
    
    ##### INLET LINES #####
    if 1 == 1:
        cam1m = (YIN[1]-YIN[0])/(XIN[1]-XIN[0]) #Inlet camber gradient
        
        ##### INLET AXIAL LINE #####
        if 1 == 1:
            ax1 = np.linspace(XIN[0]-0.2*cx,XIN[0],n)
            plt.plot(ax1,[YIN[0]+ptc*cx for i in ax1],'black',linestyle=dash,linewidth=lw)
        
        ##### INLET METAL ANGLE #####
        if 1 == 1:
            cam1x = np.linspace(XIN[0]-0.2*cx*np.cos(np.arctan(cam1m)),XIN[0],n)
            plt.plot(cam1x,[YIN[0]+ptc*cx+cam1m*(i-XIN[0]) for i in cam1x],'black',linestyle=dash,linewidth=lw)
            
        ##### INLET VELOCITY #####
        if 2 == 1:
            v1x = XIN[0]-0.2*cx*np.cos(np.arctan(cam1m)+np.radians(i))
            v1y = YIN[0]+ptc*cx-0.2*cx*np.sin(np.arctan(cam1m)+np.radians(i))
            v1dx = abs(v1x-XIN[0])-arrows*np.cos(np.arctan(cam1m)+np.radians(i))
            v1dy = abs(v1y-YIN[0]-ptc*cx)-arrows*np.sin(np.arctan(cam1m)+np.radians(i))
            guide2x = np.linspace(XIN[0]-0.2*cx*np.cos(np.arctan(cam1m)+np.radians(i)),XIN[0],n)
            guide2y = [YIN[0]+ptc*cx+(v1dy/v1dx)*(h-XIN[0]) for h in guide2x]
            plt.plot(guide2x,guide2y,'black',linewidth=0.1)
        
        ##### LEADING EDGE THICKNESS #####
        if 2 == 1:
            frac = 0.1
            l_frac = 0.1
            lexl = np.linspace(XL[int(points*frac)]-l_frac*cx,XL[int(points*frac)],n)
            perp_x = np.linspace(XIN[0],lexl[-1],n)
            perp_y = [YL[int(points*frac)]+(-1/cam1m)*(l-XL[int(points*frac)]) for l in perp_x]
            min_dis = 100
            po = 0
            for j in range(len(YU[:np.argmax(YU)])):
                for k in range(len(perp_x)):
                    if np.sqrt((perp_x[k]-XU[j])**2+(perp_y[k]-YU[j])**2)<min_dis:
                        po = j
                        min_dis = np.sqrt((perp_x[k]-XU[j])**2+(perp_y[k]-YU[j])**2)
            
            lexu = np.linspace(XU[po]-l_frac*cx,XU[po],n)
            plt.plot(lexl,[YL[int(points*frac)]+cam1m*(i-XL[int(points*frac)]) for i in lexl],'black',linestyle=dash,linewidth=lw)
            plt.plot(lexu,[YU[po]+cam1m*(i-XU[po]) for i in lexu],'black',linestyle=dash,linewidth=lw)

    ##### EXIT LINES #####
    if 1 == 1:
        cam2m = (YIN[-1]-YIN[-2])/(XIN[-1]-XIN[-2])
        endx = XL[-1]
        endy = YL[-1]
        
        ##### STAGGER #####
        if 2 == 1:
            stagm = (YIN[0]-endy)/(XIN[0]-endx)
            stagx = np.linspace(XIN[0]-0.2*cx*np.cos(np.arctan(stagm)),XIN[0],n)
            plt.plot(stagx,[YIN[0]+ptc*cx+stagm*(i-XIN[0]) for i in stagx],'black',linestyle=dash,linewidth=lw)

        ##### EXIT METAL ANGLE #####
        if 1 == 1:
            cam2x = np.linspace(endx,endx+0.2*cx*np.cos(np.arctan(cam2m)),n)
            plt.plot(cam2x,[YL[-1]+ptc*cx+cam2m*(i-XL[-1]) for i in cam2x],'black',linestyle=dash,linewidth=lw)
    
        ##### EXIT VELOCITY #####
        if 2 == 1:
            v2x = endx-0.2*cx*np.cos(np.arctan(cam2m)+np.radians(d))
            v2y = endy+ptc*cx-0.2*cx*np.sin(np.arctan(cam2m)+np.radians(d))
            v2dx = abs(v2x-endx)-arrows*np.cos(np.arctan(cam2m)+np.radians(d))
            v2dy = abs(v2y-endy-ptc*cx)+arrows*np.sin(np.arctan(cam2m)+np.radians(d))
            guide3x = np.linspace(endx,endx+0.2*cx*np.cos(np.arctan(cam2m)+np.radians(d)),n)
            guide3y = [endy+ptc*cx-(v2dy/v2dx)*(h-endx) for h in guide3x]
            plt.plot(guide3x,guide3y,'black',linewidth=0.1)
    
        ##### EXIT AXIAL LINE #####
        if 1 == 1:
            ax2 = np.linspace(endx,endx+(0.2+0.02)*cx,n)
            plt.plot(ax2,[endy+ptc*cx for i in ax2],'black',linestyle=dash,linewidth=lw)
            plt.plot(ax2,[endy for i in ax2],'black',linestyle=dash,linewidth=lw)
            
        ##### TRAILING EDGE THICKNESS #####
        if 1 == 1:
            texl = np.linspace(XL[points],XL[points]+0.05*cx,n)
            texu = np.linspace(XU[points],XU[points]+0.05*cx,n)
            plt.plot(texl,[YL[points]+cam2m*(i-XL[points]) for i in texl],'black',linestyle=dash,linewidth=lw)
            plt.plot(texu,[YU[points]+cam2m*(i-XU[points]) for i in texu],'black',linestyle=dash,linewidth=lw)

    ##### MAXIMUM THICKNESS CIRCLE #####
    if 2 == 1:
        thkx = XIN[np.argmax(thk)]
        thky = YIN[np.argmax(thk)]
        thk_circ_r = 0.92*thk[np.argmax(thk)]
        theta = np.linspace(0, np.pi*2, 500)
        thkcx = [thk_circ_r*np.cos(i)+thkx for i in theta]
        thkcy = [thk_circ_r*np.sin(i)+thky for i in theta]
        plt.plot(thkcx,thkcy,'black', linewidth=0.5)
        
    ##### MAXIMUM CAMBER POINT #####
    if 2 == 1:
        cam_max = np.argmax(YIN)
        cam_maxy = np.linspace(YIN[0],YIN[cam_max],n)
        plt.plot([XIN[cam_max]for i in cam_maxy],cam_maxy,'black',linestyle='--',dashes=(5,4),linewidth=lw)

    ##### AXIAL CHORD #####
    if 1 == 1:
        cx1 = np.linspace(YU[np.argmin(XU)]+ptc*cx,np.amax(YU)+(ptc+0.06)*cx,n)
        cx2 = np.linspace(TEcy[np.argmax(TEcx)]+ptc*cx,np.amax(YU)+(ptc+0.06)*cx,n)
        plt.plot([np.amin(X) for i in cx1],cx1,'0.5',linestyle=dash,linewidth=lw)
        plt.plot([np.amax(X) for i in cx2],cx2,'black',linestyle=dash,linewidth=lw)
    
    plt.show()
    plt.savefig('Blades.eps',format='eps',bbox_inches='tight',transparent=True)

if 2==1:			### Plot of a stage
    
    #Define lenghts
    vx = 0.5*cx
    edge_space = 0.2
    sx = 1.8*cx
    sy = 0.8
    ptc = 0.9
    
    ##### AXIAL VELOCITY LINES #####
    
    vx1 = np.linspace(np.amin(X)-vx-edge_space*cx,np.amin(X)-edge_space*cx,n)
    vx2 = np.linspace(np.amax(X)+(sx-cx-vx)/2, np.amin([i+sx for i in X])-(sx-cx-vx)/2, n)
    vx3 = np.linspace(np.amax(X)+edge_space*cx+sx,np.amax(X)+vx+edge_space*cx+sx,n)
    y_pos = YIN[0]+0.4*ptc*cx #y position where lines are aligned
    
    ##### OTHER VELOCITY LINES #####
    
    v1y = [y_pos + np.tan(-np.radians(X1))*(i-vx1[0]) for i in vx1]
    v2y = [y_pos + np.tan(-np.radians(X2))*(i-vx2[0]) for i in vx2]
    v3y = [y_pos + np.tan(-np.radians(X1))*(i-vx3[0]) for i in vx3]
    w2y = [y_pos + np.tan(np.radians(X1))*(i-vx2[0]) for i in vx2]
    w3y = [y_pos + np.tan(np.radians(X2))*(i-vx3[0]) for i in vx3]
    
    ##### BLADE VELOCITY LINE #####
    
    U = abs(v2y[-1]-w2y[-1])
    uy = np.linspace(np.amax(YU)+sy+ptc*cx+0.05*cx,np.amax(YU)+sy+ptc*cx+0.05*cx+U,n)
    ux = [XU[np.argmax(YU)]+cx*sx for i in uy]
    
    ##### +/- LINES #####
    
    anghx = np.linspace(np.amin(X)+0.2*cx,np.amin(X)+0.55*cx,n)
    anghy = [YIN[0]+sy+ptc*cx for i in anghx]
    angpx = np.linspace(np.amin(X)+0.2*cx,np.amin(X)+0.45*cx,n)
    angpy = [anghy[0]+np.tan(np.radians(40))*(i-anghx[0]) for i in angpx]
    angnx = np.linspace(np.amin(X)+0.2*cx,np.amin(X)+0.45*cx,n)
    angny = [anghy[0]+np.tan(np.radians(-40))*(i-anghx[0]) for i in angnx]
    
    plt.axis('equal')
    plt.margins(0.05,0.05)
    plt.plot(vx1,[y_pos for i in vx1],'black',linewidth=0.1)
    plt.plot(vx2,[y_pos for i in vx2],'black',linewidth=0.1)
    plt.plot(vx3,[y_pos for i in vx3],'black',linewidth=0.1)
    plt.plot(vx1,v1y,'black',linewidth=0.1)
    plt.plot(vx2,v2y,'black',linewidth=0.1)
    plt.plot(vx3,v3y,'black',linewidth=0.1)
    plt.plot(vx2,w2y,'black',linewidth=0.1)
    plt.plot(vx3,w3y,'black',linewidth=0.1)
    plt.plot(ux,uy,'black',linewidth=0.1)
    plt.plot(anghx,anghy,'black',linewidth=0.1)
    plt.plot(angpx,angpy,'black',linewidth=0.1)
    plt.plot(angnx,angny,'black',linewidth=0.1)

    plt.plot(X,[-i for i in Y],'black',linewidth=2)
    plt.plot(X,[i+ptc*cx for i in [-i for i in Y]],'black',linewidth=2)
    plt.plot([i+sx for i in X],[i+sy for i in Y],'black',linewidth=2)
    plt.plot([i+sx for i in X],[i+ptc*cx+sy for i in Y],'black',linewidth=2)
    plt.axis('off')
    plt.show()
    plt.savefig('Stage.eps',format='eps',bbox_inches='tight',transparent=True)
