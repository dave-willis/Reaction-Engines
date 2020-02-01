"""Module contains functions generating blade profiles from CJC and length, area
and plotting functions from WFD"""

#Import required modules
import math
import numpy as np
import scipy.optimize as SciOpt
from scipy.integrate import simps
import csv

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
    ### NOTE: XP/XS ETC DON'T ALWAYS REFER TO THE PRESSURE/SUCTION SURFACEC=S
    ### 'S' VARIABLES ARE THE LOWER SURFACE AND 'P' THE UPPER, SO IF EXIT ANGLE
    ### IS GREATER THAN INLET THIS WILL BE THE SUCTION AND PRESSURE SURFACES

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

################################
###SS LENGTH AND SECTION AREA###
################################

def blade_dims(X1, X2, TE, Cx):
    """Return the suction surface length and section area of the normalised blade"""

    #Check if exit angle is less than the inlet angle. If so, invert them so
    #the Profile function returns XS and YS as the SS not the PS
    if X2 < X1:
        X2 = -X2
        X1 = -X1
    #Set number of points on the surfaces
    n_points = 40
    #Apply profile function to get the surface points
    surface = Profile(X1, X2, TE/Cx, Cx, points=n_points)
    #Extract the two surfaces and offset them upwards so all y values are positive
    min_y = np.amin(np.concatenate((surface[1], surface[3])))
    XP = surface[0]
    YP = [i+1.01*abs(min_y) for i in surface[1]]
    XS = surface[2]
    YS = [i+1.01*abs(min_y) for i in surface[3]]
    #Calculate the blade area by calculating the area under the pressure and
    #suction surfaces with Simpson's rule then finding the difference
    AP = simps(YP, x=XP)
    AS = simps(YS, x=XS)
    A_blade = AP-AS
    #Calculate the length of the suction surface
    ss_len = 0
    for i in range(n_points-1):
        ss_len += np.hypot(XS[i+1]-XS[i], YS[i+1]-YS[i])

    return A_blade, ss_len

def make_grid():
    """Make .csv files with grids of normalised SS lengths and section areas"""

    #Define lists of variables
    xin = []
    xout = []
    sslen = []
    profile_area = []
    #Loop over reasonable ranges of inlet and exit angle, ensure exit is 
    #greater than inlet for convenience
    for a1 in range(-60,62,2):
        for a2 in range(-80,82,2):
            x1 = a1
            x2 = a2
            if a2<a1:
                x1 = -a1
                x2 = -a2
            #Calculate SS length and blade area, with a representative TE/Cx
            calc = blade_dims(x1, x2, 0.0003/0.01, 1.0)
            xin.append(a1)
            xout.append(a2)
            sslen.append(calc[1])
            profile_area.append(calc[0])
    #Add a blank space to the start of xin for the top left corner of the grid
    xout.insert(0,'')
    #Create a grid for the SS lengths
    with open('ss_grid.csv', mode='w') as ss_grid:
        table_writer = csv.writer(ss_grid, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        #The top row of exit angles
        table_writer.writerow(xout[:82])
        #Create other rows
        for i in range(int(len(xin)/81)):
            row = sslen[81*i:81*(i+1)]
            row.insert(0, xin[81*i])
            table_writer.writerow(row)
    #Create a grid for the areas
    with open('Profile_areas.csv', mode='w') as areas:
        table_writer = csv.writer(areas, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        #The top row of exit angles
        table_writer.writerow(xout[:82])
        #Create other rows
        for i in range(int(len(xin)/81)):
            row = profile_area[81*i:81*(i+1)]
            row.insert(0, xin[81*i])
            table_writer.writerow(row)

# make_grid()
from scipy import interpolate as sciint
#Load the grid of normalised suction surface lengths and create the interpolation function
xin = []
xout = []
sslen = []
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
xout = np.asarray(xout)
xin = np.asarray(xin)
sslen = np.asarray(sslen)
ss_length = sciint.RectBivariateSpline(xin, xout, sslen, kx=1, ky=1)

print(0.006*ss_length(33,75))
print(blade_dims(33,75,0.0003,0.006)[1])

# with open("blade_table.csv") as csvfile:
#     reader = csv.reader(csvfile, quoting=csv.QUOTE_NONNUMERIC) # change contents to floats
#     for row in reader: # each row is a list
#         table.append(row)
#         xin.append(row[0])
#         xout.append(row[1])
#         sslen.append(row[2])
#         profile_area.append(row[3])
        
# table = []       
# with open("ss_grid.csv") as csvfile:
#     reader = csv.reader(csvfile, quoting=csv.QUOTE_NONNUMERIC) # change contents to floats
#     for row in reader: # each row is a list
#         table.append(row)
        
# for i in range(1, len(table)):
#     for j in range(1, len(table[0])):
#         if table[i][0] > table[0][j]:
#             table[i][j] = table[len(table)-i][len(table[0])-j]

# with open('ss_grid.csv', mode='w') as ss_grid:
#     table_writer = csv.writer(ss_grid, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
#     for row in table:
#         table_writer.writerow(row)

# with open("ss_grid.csv") as csvfile:
#     reader = csv.reader(csvfile, quoting=csv.QUOTE_NONNUMERIC) # change contents to floats
#     line = 0
#     for row in reader: # each row is a list
#         table.append(row)
#         if line == 0:
#             xout = row[1:]
#             line += 1
#         else:
#             xin.append(row[0])
#             sslen.append(row[1:])
#             line += 1
