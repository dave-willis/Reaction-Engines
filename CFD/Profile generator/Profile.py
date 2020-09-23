import numpy as np
import matplotlib.pyplot as plt
import math
import scipy.optimize as SciOpt

def binomial(i, n):
    """Binomial coefficient"""

    return math.factorial(n) / float(
        math.factorial(i) * math.factorial(n - i))

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
    
def GEN_TURNING(Xi1,Xi2,controls,plotting):
    """Generate x,y of a bezier from control points"""

    spacings = np.linspace(0,1,len(controls)+2)[1:-1]
    CPs = []
    CPs.append([0.0,0.0])   
    for i in range(len(controls)):
        CPs.append([spacings[i]+min(spacings[i],1.-spacings[i])*controls[i],spacings[i]-min(spacings[i],1.-spacings[i])*controls[i]])
    CPs.append([1.0,1.0])

    steps = 100
    controlPoints=CPs
    oldPoint = controlPoints[0]
    x=[]
    y=[]
    for point in bezier_curve_range(steps, controlPoints):
        ##printpoint
        x.append(point[0])
        y.append(point[1])
   
    ### integrate profile
    xcam = [0.]
    ycam = [0.]
    
    DXi = Xi2-Xi1
    for i in range(1,len(x)):
        Xi = Xi1+DXi*(y[i]+y[i-1])*0.5
        xcam.append(xcam[-1]+(x[i]-x[i-1])*np.cos(Xi))
        ycam.append(ycam[-1]+(x[i]-x[i-1])*np.sin(Xi))

    DX=xcam[-1]
    ycam=ycam/DX
    xcam=xcam/DX

    return(np.array(xcam),np.array(ycam))

def GEN_BZ_2D(controls_x,controls_y,plotting):
    """Generate x,y of a bezier from control points"""

    CPs = []
    for i in range(len(controls_x)):
        CPs.append([controls_x[i],controls_y[i]])

    steps = 100
    controlPoints=CPs
    oldPoint = controlPoints[0]
    x=[]
    y=[]
    for point in bezier_curve_range(steps, controlPoints):
        ##printpoint
        x.append(point[0])
        y.append(point[1])          

    return(np.array(x),np.array(y),CPs)

def dist_vino(n, xst, xen, s1, s2):
    """vinokur distribution"""

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
    """function to generate class function of psi"""

    C = (psi**0.5)*(1.0-psi)
    return(C)
    
def F_TF_bez(controls_cam_x,controls_cam_y,TE,Oval_frac,plotting,spacing):
    """function that returns thickness contributions due to parameters"""

    psi = dist_vino(200, 0, 1, 1./2000, 1./500)
    X,S,cp = GEN_BZ_2D(controls_cam_x,controls_cam_y,plotting)
    S = np.interp(psi,X,S)
    
    S=S+Oval_frac*S[0]*(1.-psi)**18
    C=F_CF(psi)		### calculate class function

    TF =S*C	+psi*TE/2.	### combine to form thickness funtion

    return(TF)
        
def calc_norm(x_cam,y_cam,plotting):
    """Calculate the normal to the camber line"""
    
    s=np.zeros((len(x_cam)))			### initialise camber length array
    grad = np.zeros((len(x_cam),2))		### initialise camber gradient array
    norm = np.zeros((len(x_cam),2))		### initialise camber normals array
    for i in range(1,len(x_cam)):
        s[i]=s[i-1]+((x_cam[i]-x_cam[i-1])**2+(y_cam[i]-y_cam[i-1])**2)**0.5		### calculate cumsum length of camber
    grad[:-1,0]=(x_cam[1:]-x_cam[:-1])/(s[1:]-s[:-1])					### calculate dx/ds
    grad[:-1,1]=(y_cam[1:]-y_cam[:-1])/(s[1:]-s[:-1])				### calculate dy/ds
    grad[-1,:]=grad[-2,:]								### deal with final point
    norm[:,0]=-grad[:,1]								### calculate normal 1 as -1/grad2
    norm[:,1]=grad[:,0]

    return(norm)

def F_Make(Xi1,Xi2,controls_cam,controls_thk_x,controls_thk_y,Tte,Oval_frac,plotting):
    """function that compiles a blade from reduced variables"""

    ### turning distribution
    Xi1 = np.radians(Xi1)
    Xi2 = np.radians(Xi2)
    (x_cam,y_cam)=GEN_TURNING(Xi1,Xi2,controls_cam,plotting)      ### generate camber line from bezier curve

    Gamma = np.degrees(np.arctan((y_cam[-1])/x_cam[-1]))
    #print'Xi1:',Xi1,'Xi2:',Xi2,'Gamma:',Gamma

    s2=np.zeros((len(x_cam)))			### initialise camber length array
    for i in range(1,len(x_cam)):
            s2[i]=s2[i-1]+((x_cam[i]-x_cam[i-1])**2+(y_cam[i]-y_cam[i-1])**2)**0.5	
    psi = s2/s2.max()
    
    psi_new = dist_vino(200, 0, 1, 1./2000, 1./500)   
    
    x_cam = np.interp(psi_new,psi,x_cam)    
    y_cam = np.interp(psi_new,psi,y_cam)   
    psi =psi_new
    norm = calc_norm(x_cam,y_cam,plotting=False)
    
    thk = F_TF_bez(controls_thk_x,controls_thk_y,Tte,Oval_frac,plotting=0,spacing=None)
    
    Z_U = y_cam+thk*norm[:,1]	### apply thickness onto camber line for upper
    Z_L = y_cam-thk*norm[:,1]	### apply thickness onto camber line for lower
    X_U = x_cam+thk*norm[:,0]	### repeat upper for x rather than y
    X_L = x_cam-thk*norm[:,0]

    X = np.append(X_L[::-1],X_U[1:])	### combine axial distributions upper and lower
    Z = np.append(Z_L[::-1],Z_U[1:])	### combine tangential distributions upper and lower

    return(X_L,Z_L,X_U,Z_U,thk)

def Profile2(X1,X2,TKTE):
    """function based on profgen.f (JDD) - GENERATES SIMPLE BLADE PROFILES"""

    ### THIS FUNCTION SHOULD BE REPLACED AT SOME POINT TO USE A BETTER BLADE SECTION GENERATOR
    ### EXCUSE THE WEIRD FORMULATION, THIS IS CONVERSION FROM FORTRAN
    ### USER HARD CODED VALUES BELOW - COULD BE TAKEN OUTSIDE FUNCTION CJC

    controls_cam = [-0.0,-0.5,-0.8]
    Rle = 0.05
    Tte = TKTE
    Beta_te = 4.0
    Oval_frac = 0.3
    #controls_cam = [0,0.2,0.0]
    controls_thk_x = [0.0,0.5,1.0]
    controls_thk_y = [(2.*Rle)**0.5*(1.-Oval_frac),0.25,np.tan(np.radians(Beta_te))+Tte/2.]

    (XSIN,YSIN,XPIN,YPIN,thk)= F_Make(X1,X2,controls_cam,controls_thk_x,controls_thk_y,Tte,Oval_frac,plotting=False)
    
    XIN = (XPIN+XSIN)/2.
    YIN = (YPIN+YSIN)/2.
    XScale = max(max(XPIN),max(XSIN))
    XPIN = XPIN/XScale
    XSIN = XSIN/XScale
    XIN = XIN/XScale
			
    yoffset = np.mean(YIN)
    YPIN=(YPIN-yoffset)/XScale	### center on peak thickness - DEFAULTS STACKING TO PEAK THICKNESS
    YSIN=(YSIN-yoffset)/XScale	### center on peak thickness - DEFAULTS STACKING TO PEAK THICKNESS
    YIN = (YIN - yoffset)/XScale ### center on peak thickness - DEFAULTS STACKING TO PEAK THICKNESS
    
    Xnew =np.zeros((len(XPIN)+len(XSIN)-1))
    Xnew[:200]=XPIN[::-1]
    Xnew[199:]=XSIN[:]

    Ynew =np.zeros((len(YPIN)+len(YSIN)-1))
    Ynew[:200]=YPIN[::-1]
    Ynew[199:]=YSIN[:]

    Snew = np.zeros((len(Ynew)))
    for i in range(1,len(Xnew)):
        Snew[i]=Snew[i-1]+((Xnew[i]-Xnew[i-1])**2.+(Ynew[i]-Ynew[i-1])**2.)**0.5

    for i in range(1,len(Xnew)):
        if Xnew[i]==Xnew.min():
            Sle=Snew[i]
            Xle=Xnew[i]
            break
    
    points = 500
    XPIN = np.interp(np.linspace(0,Sle,points),Snew,Xnew)[::-1]
    YPIN = np.interp(np.linspace(0,Sle,points),Snew,Ynew)[::-1]
    XSIN = np.interp(np.linspace(Sle,Snew[-1],points),Snew,Xnew)
    YSIN = np.interp(np.linspace(Sle,Snew[-1],points),Snew,Ynew)

    XScale = max(max(XPIN),max(XSIN))-Xle
    XPIN = (XPIN-Xle)/XScale
    XSIN = (XSIN-Xle)/XScale
    XIN = (XIN-Xle)/XScale
			
    YPIN=(YPIN-yoffset)/XScale	### center on peak thickness - DEFAULTS STACKING TO PEAK THICKNESS
    YSIN=(YSIN-yoffset)/XScale	### center on peak thickness - DEFAULTS STACKING TO PEAK THICKNESS
    YIN =(YIN - yoffset)/XScale ### center on peak thickness - DEFAULTS STACKING TO PEAK THICKNESS

    TE_circ_r = 0.5*((abs(XSIN[-1]-XPIN[-1]))**2+(abs(YSIN[-1]-YPIN[-1]))**2)**0.5
    TE_circ_cx = 0.5*(XSIN[-1]+XPIN[-1])
    TE_circ_cy = 0.5*(YSIN[-1]+YPIN[-1])
    TEcx = np.arange(XSIN[-1],TE_circ_cx+TE_circ_r,0.00001)
    TEcx2 = np.arange(XPIN[-1],TE_circ_cx+TE_circ_r,0.00001)
    TEcy = np.asarray([TE_circ_cy-(TE_circ_r**2-(i-TE_circ_cx)**2)**0.5 for i in TEcx])
    TEcy2 = np.asarray([TE_circ_cy+(TE_circ_r**2-(i-TE_circ_cx)**2)**0.5 for i in TEcx2])

    maxes = np.array([np.amax(XIN),np.amax(XPIN),np.amax(XSIN),np.amax(TEcx),np.amax(TEcx2)])
    mins = np.array([np.amin(XIN),np.amin(XPIN),np.amin(XSIN)])
    cx = np.amax(maxes)-np.amin(mins)
    ptc = 1.1
    plt.show()
    
    if 1==1:            ### Plot with lines
        
        plt.plot(XPIN,YPIN,'black',linewidth=2)
        plt.plot(XSIN,YSIN,'black', linewidth=2)
        plt.plot(TEcx,TEcy,'black', linewidth=2)
        plt.plot(TEcx2,TEcy2,'black', linewidth=2)
        plt.plot(XPIN,[i+ptc*cx for i in YPIN],'black',linewidth=2)
        plt.plot(XSIN,[i+ptc*cx for i in YSIN],'black', linewidth=2)
        plt.plot(TEcx,[i+ptc*cx for i in TEcy],'black', linewidth=2)
        plt.plot(TEcx2,[i+ptc*cx for i in TEcy2],'black', linewidth=2)
        plt.axis('off')
        plt.axis('equal')
        plt.margins(0.05,0.05)
    
    delta = 0.5*abs(XSIN[-1]-XPIN[-1])		### RECALCULATE TRAILING EDGE THICKNESS
    print XIN
    return(XIN,YIN,XPIN,YPIN,XSIN,YSIN,TEcx,TEcy,TEcx2,TEcy2)

k = Profile2(0,70,0.03)
#plt.show()

