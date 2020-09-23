import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os,sys,time,copy,math,scipy,numpy
from scipy import interpolate as sciiint
import ts_tstream_default
from ts import ts_tstream_grid, ts_tstream_type,ts_tstream_load_balance
import pylab
import scipy.optimize as SciOpt

#########################################
### MODULAR GEOMETRY AND MESHING TOOL ###
### WRITTEN BY CJC95 2017	      ###
### UPDATED BY J WALDREN 2018	      ###
### UPDATED BY N CLARK 2019	      ###
### CHECKED AGAIN CJC95 2019	      ###
#########################################






class ROW:		### EACH ROW IS DEFINED BY A "ROW" CLASS THAT STORES ITS DETAILS
    def __init__(self): 
	### THIS INIT sets some values to defaults - most are overwritten on generation
	self.name	=	''		### option to name row
	self.X1=[]				### INLET METAL ANGLES
	self.X2=[]				### EXIT METAL ANGLES
	self.Cx=[]				### AXIAL CHORDS
	self.RADIUS =[]				### RADIUS OF SECTIONS
	self.RPM = 0.				### RPM OF ROW
	self.Nblades = 1.0			### # of blades in row
	self.shroud_hub = False			### if there is a hub shroud
	self.shroud_cas = False			### if there is a casing shroud
	self.upfrac=[0.5,0.3,0.2]		### fraction of chord upstream for [mx,shroudstart,shroudend]
	self.downfrac=[0.5,0.2,0.3]		### fraction of chord upstream for [mx,shroudend,shroudstart]
	self.shroudradius = 0.005	### thickness of shroud [upstream,downstream]
	self.shroudclearance = 0.002	### shroud clearance at upstream and downstream edge
	self.Nfins = 0
	self.finwidth = 0.05		### fin width as fraction of shroud length
	self.finclearance = 0.0025	### fin clearance

def set_aspect_equal_3d(ax):	### THIS IS A QUICK FUNCTION TO SET THE ASPECT RATIO OF 3D PLOTS TO BE EQUAL
	xlim=ax.get_xlim3d()
	ylim=ax.get_ylim3d()
	zlim=ax.get_zlim3d()

	xmean=np.mean(xlim)
	ymean=np.mean(ylim)
	zmean=np.mean(zlim)
	
	plot_radius=max([abs(lim-mean_)
				for lims,mean_ in 	((xlim,xmean),
							(ylim,ymean),
							(zlim,zmean))
				for lim in lims])
	ax.set_xlim3d([xmean-plot_radius,xmean+plot_radius])
	ax.set_ylim3d([ymean-plot_radius,ymean+plot_radius])
	ax.set_zlim3d([zmean-plot_radius,zmean+plot_radius])


def vinokur_dist( xst, xen,n, s1, s2):	### vinokur distribution 
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
    #print trans(0.00001), b, math.
    left = 0.00001
    delta = SciOpt.brentq(trans, left, 100)
    #delta = scipy.optimize.fmin_powell(trans, trans(0.0001), disp=False)
    u = 0.5*(1.0 + np.tanh(delta*(eta/nm1 - 0.5))/(np.tanh(delta*0.5)))
    s = u/(a + (1-a)*u)
    return xst + (xen-xst)*s



def vinokur_dist_NC(x0,x1,n,spac1,spac2):   ### FUNCTION TO PRODUCE A VINOKUR DISTRIBUTION WRITTEN BY N CLARK 2018
    n = n - 1
    frac1 = spac1/(x1-x0)
    frac2 = spac2/(x1-x0)
    delS1 = frac1*n
    delS2 = frac2*n
    A = math.sqrt(delS2/delS1)              ## calculation of delta from requred end spacings
    B = 1/(math.sqrt(delS1*delS2))
    if B == 1:
        delta = 0.0000000000000001
    elif B <= 2.7829681:                   ## inverse calculation of y = sinhx/x
        BDash = B - 1
        delta = math.sqrt(6*BDash)*(1 - 0.15*BDash + 0.057321429*BDash**-2 - 0.024907295*BDash**-3 + 0.0077424461*BDash**-4 - 0.0010794123*BDash**-5)
    else:
        v = math.log(B)
        w = 1/B - 0.028527431
        delta = v+(1+1/v)*math.log(2*v)-0.02041793+0.24902722*w+1.9496443*w**2-2.6294547*w**3+8.56795911*w**4
    x = np.linspace(0,1,n+1)                               ## set of n equally spaced point between 0 and 1
    u = 0.5*(1+np.tanh(delta*(x-0.5))/np.tanh(delta/2))    ## adjusts points to vinokur distrubution based on input scaling factor
    s = u/(A+(1-A)*u)                                      ## adds different spacings at either end
    s = s*(x1-x0) + x0                                     ## scales to desired range
#    if abs(s[1]-s[0])>math.sqrt(1.3)*frac1*abs(x1-x0) or abs(s[1]-s[0])<math.sqrt(0.77)*frac1*abs(x1-x0):
    if abs(s[1]-s[0])>1.3*frac1*abs(x1-x0) or abs(s[1]-s[0])<0.77*frac1*abs(x1-x0):
        print('WARNING mesh size at upstream boundary, scale should be between 1.14 and 0.877')	### WARNINGS ON EXPANSION RATIOS (KEEP WITHIN 1.2)
        print(abs(s[1]-s[0])/(frac1*abs(x1-x0)))
#    if abs(s[-1]-s[-2])>math.sqrt(1.3)*frac2*abs(x1-x0) or abs(s[-1]-s[-2])<math.sqrt(0.77)*frac2*abs(x1-x0):
    if abs(s[-1]-s[-2])>1.3*frac2*abs(x1-x0) or abs(s[-1]-s[-2])<0.77*frac2*abs(x1-x0):
        print('WARNING mesh size at downstream boundary, scale should be between 1.14 and 0.877') ### WARNINGS ON EXPANSION RATIOS (KEEP WITHIN 1.2)
        print(abs(s[-1]-s[-2])/(frac2*abs(x1-x0)))
    for i in range(2,len(s)):
        if abs(s[i]-s[i-1])>1.3*abs(s[i-1]-s[i-2]) or abs(s[i]-s[i-1])<0.77*abs(s[i-1]-s[i-2]):
            print('WARNING scaling gradient of mesh too large, scale should be between 1.3 and 0.77') 	### WARNINGS ON EXPANSION RATIOS (KEEP WITHIN 1.2)
            print(abs(s[i]-s[i-1])/abs(s[i-1]-s[i-2]))
    return s



############################################
### NEW ADDITIONS FOR PROFILE GENERATION ###
############################################


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
    for i in xrange(n):
        t = i / float(n - 1)
        yield bezier(t, points)

    
def GEN_BZ(controls,plotting):
    """ Generate x,y of a bezier from control points """
    
    spacings = np.linspace(0,1,len(controls)+2)[1:-1]

    CPs = []
    CPs.append([0.0,0.0])   
    for i in range(len(controls)):
        CPs.append([spacings[i],controls[i]])
    CPs.append([1.0,1.0])

    steps = 100    
    controlPoints=CPs
   
    oldPoint = controlPoints[0]
    x=[]
    y=[]
    for point in bezier_curve_range(steps, controlPoints):
        #print point
        x.append(point[0])
        y.append(point[1])          
        
    if plotting == True:
        plt.figure('Geometry Definition')
        plt.subplot(221)
        plt.plot(x,y,'-b')
        oldPoint=controlPoints[0]
        for point in controlPoints:
            plt.plot([oldPoint[0],point[0]],[oldPoint[1],point[1]],'--k')
            plt.plot(point[0],point[1],'ok')
            oldPoint=point
        plt.axis('equal')
        plt.xlabel('surface fraction')
        plt.ylabel('turning fraction')
        plt.show()
    return(np.array(x),np.array(y),CPs)
    
def GEN_TURNING(Xi1,Xi2,controls,plotting):
    """ Generate x,y of a bezier from control points """
    
    spacings = np.linspace(0,1,len(controls)+2)[1:-1]
    print spacings

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
        #print point
        x.append(point[0])
        y.append(point[1])          
        
    if plotting == True:
        plt.figure('Geometry Definition')
        plt.subplot(221)
        plt.plot(x,y,'-b')
        oldPoint=controlPoints[0]
        for point in controlPoints:
            plt.plot([oldPoint[0],point[0]],[oldPoint[1],point[1]],'--k')
            plt.plot(point[0],point[1],'ok')
            oldPoint=point
        plt.axis('equal')
        plt.xlabel('surface fraction')
        plt.ylabel('turning fraction')
        #plt.axis([0,1,0,1])
       
        
        
    ### integrate profile
    xcam = [0.]
    ycam = [0.]
    
    DXi = Xi2-Xi1
    for i in range(1,len(x)):
        Xi = Xi1+DXi*(y[i]+y[i-1])*0.5
        xcam.append(xcam[-1]+(x[i]-x[i-1])*np.cos(Xi))
        ycam.append(ycam[-1]+(x[i]-x[i-1])*np.sin(Xi))
        
    if False == True:
        plt.figure('Geometry Definition')
        plt.subplot(224)
        plt.plot(xcam,ycam,'-b')

        plt.axis('equal')
        #plt.show()
        
    DX=xcam[-1]
    ycam=ycam/DX
    xcam=xcam/DX
        
    return(np.array(xcam),np.array(ycam))

def GEN_BZ_2D(controls_x,controls_y,plotting):
    """ Generate x,y of a bezier from control points """

    CPs = []
    for i in range(len(controls_x)):
        CPs.append([controls_x[i],controls_y[i]])

    steps = 100    
    controlPoints=CPs
    oldPoint = controlPoints[0]
    x=[]
    y=[]
    for point in bezier_curve_range(steps, controlPoints):
        #print point
        x.append(point[0])
        y.append(point[1])          
        
    if plotting == True:
        plt.figure('Geometry Definition')
        plt.subplot(222)
        plt.plot(x,y,'-r')
        oldPoint=controlPoints[0]
        for point in controlPoints:
            plt.plot([oldPoint[0],point[0]],[oldPoint[1],point[1]],'--m')
            plt.plot(point[0],point[1],'om')
            oldPoint=point
        #plt.axis('equal')
        #plt.xlabel('axial - camber')
        #plt.ylabel('tangential - camber')
    return(np.array(x),np.array(y),CPs)
    
def GEN_RMSE(controls,Yref):
    """ generate a RMSE between a reference Y and a bezier defined by control points"""
    (xn,yn,CPn)=GEN_BZ(controls)
    return(RMSE(yn,Yref))


    
def F_CF(psi): # function to generate class function of psi
    C = (psi**0.5)*(1.0-psi)
    return(C)

    
def F_SC(psi,i,Nm1): # function to generate shape component
    C = (psi**i)*((1.0-psi)**(Nm1-i))
    return(C)
    
def F_SF(psi,AIs,plotting): ### function generates Shape function based on reduced design parameters

    N = len(AIs)	### number of coefficients
    Nm1=N-1		### N-1
    npts = len(psi)	### number of points in distribution
    S = np.zeros((npts)) ### initilise shape space distribution

    #plotting = 1 ### optional plot out
    if plotting==1:
        plt.figure('Geometry Definition')
        plt.subplot(222)
        #plt.title('shape function')
        plt.xlabel('psi')
        plt.ylabel('Shape Function: S(psi)')

    for i in range(N):	### iterate over each contirbution
        Ai=AIs[i]
        Ki = binomial(i,Nm1)		### calculate bernstein polynomial
        SC = F_SC(psi,i,Nm1)		### calculate shape component
        S[:]=S[:]+Ai*Ki*SC[:]

        if plotting==1:### another optional plot out
            plt.plot(psi,Ai*Ki*SC[:],'-',label='component'+str(i))
        

    if plotting==1:### another optional plot out
        plt.plot(psi,S,'--k',label='total of components')
        plt.plot(psi,0.3*S[0]*(1.-psi)**18,'--r')
        
    S=S+0.3*S[0]*(1.-psi)**18
    if plotting==1:
            plt.plot(psi,S,'-k',label='total of components')
	#plt.legend(loc=0)
        #plt.show()  ### uncomment to plot imediatly     
    return(S)
    
def F_TF(psi,AIs,TE,plotting,spacing): ### function that returns thickness contributions due to parameters

    if spacing == None:
    	S=F_SF(psi,AIs,plotting)	### calculate shape function
    else:
    	S=F_SF(spacing,AIs,plotting)	### calculate shape function
    #plt.show()
    #S2 = np.interp(psi,psi2,S)
    C=F_CF(psi)		### calculate class function

    TF =S*C	+psi**0.5*TE/2.	### combine to form thickness funtion
    TFMAX=2.*TF.max()
    if plotting==1:	### optional plot out
        plt.figure('Geometry Definition')
        plt.subplot(223)
        #plt.title('thickness function')
        plt.xlabel('psi')
        plt.ylabel('thickness')
        plt.plot(psi,TF/TFMAX,'-r')
        plt.plot(psi,(psi**0.5*TE/2.)/TFMAX,'--r')
        plt.plot(psi,(S*C	+psi**0.5*TE/2)/TFMAX,'--b')
        #plt.plot(psi,-TF,'-b')
        plt.axis('equal')
    return(TF)
    
def F_TF_bez(controls_cam_x,controls_cam_y,TE,Oval_frac,plotting,spacing): ### function that returns thickness contributions due to parameters

    
    psi = dist_vino(200, 0, 1, 1./2000, 1./500)  
    X,S,cp = GEN_BZ_2D(controls_cam_x,controls_cam_y,plotting)
    
 
    
    S = np.interp(psi,X,S)
    #print S
    
    if plotting==1:	### optional plot out
        plt.figure('Geometry Definition')
        plt.subplot(222)
        #plt.title('thickness function')
        plt.plot(psi,S,'--m',label='total of components')
    
    S=S+Oval_frac*S[0]*(1.-psi)**18
    if plotting==1:	### optional plot out
        plt.figure('Geometry Definition')
        plt.subplot(222)
        #plt.title('thickness function')
        plt.xlabel('psi')
        plt.ylabel('SS')
        plt.plot(psi,S,'-m',label='total of components')
    #plt.show()
    #S2 = np.interp(psi,psi2,S)
    C=F_CF(psi)		### calculate class function

    TF =S*C	+psi*TE/2.	### combine to form thickness funtion
    TFMAX=2.*TF.max()
    if plotting==1:	### optional plot out
        plt.figure('Geometry Definition')
        plt.subplot(223)
        #plt.title('thickness function')
        plt.xlabel('psi')
        plt.ylabel('thickness')
        plt.plot(psi,TF/TFMAX,'-r')
        plt.plot(psi,(psi*TE/2.)/TFMAX,'--r')
        #plt.plot(psi,(S*C	+psi**0.5*TE/2)/TFMAX,'-b')
        #plt.plot(psi,-TF,'-b')
        plt.axis('equal')
        
        print 'Peak thickness:',2.*TF.max()*100,'% Cx'
        for i in range(len(TF)):
            if TF[i]==TF.max():
                print '@ ',100.*psi[i],' % camber line'
                break
            
        Bte = np.degrees(np.arctan((TF[-2]-TF[-1])/(psi[-1]-psi[-2])))
        print 'Boat tail:',Bte
    
    return(TF)
    
    
def calc_norm(x_cam,y_cam,plotting):
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
        #print A,portiondroop
	#plotting=1							### calculate normal 2 as 1/grad1
        if plotting==1:	### optional plot out
            plt.figure('Geometry Definition')
            plt.subplot(224)
            #plt.title('Camber line - showing normal vectors')
            plt.xlabel('psi')
            plt.ylabel('zeta')
            plt.plot(x_cam,y_cam,'-k',label = 'camber')
            plt.axis('equal')
            for i in range(len(x_cam))[::5]:
                plt.plot([x_cam[i],x_cam[i]+norm[i,0]*0.01],[y_cam[i],y_cam[i]+norm[i,1]*0.01],'k')	### plot normal vectors
                plt.plot([x_cam[i],x_cam[i]-norm[i,0]*0.01],[y_cam[i],y_cam[i]-norm[i,1]*0.01],'k')
	    #plt.show()
        return(norm)    
        
def F_Make(Xi1,Xi2,controls_cam,controls_thk_x,controls_thk_y,Tte,Oval_frac,plotting):        ### function that compiles a blade from reduced variables 


    ### turning distribution
    Xi1 = np.radians(Xi1)
    Xi2 = np.radians(Xi2)
    (x_cam,y_cam)=GEN_TURNING(Xi1,Xi2,controls_cam,plotting)      ### generate camber line from bezier curve
   

    Gamma = np.degrees(np.arctan((y_cam[-1])/x_cam[-1]))
    print 'Xi1:',Xi1,'Xi2:',Xi2,'Gamma:',Gamma

    s2=np.zeros((len(x_cam)))			### initialise camber length array
    for i in range(1,len(x_cam)):
            s2[i]=s2[i-1]+((x_cam[i]-x_cam[i-1])**2+(y_cam[i]-y_cam[i-1])**2)**0.5	
    psi = s2/s2.max()
    

    
    psi_new = dist_vino(200, 0, 1, 1./2000, 1./500)   
    
    x_cam = np.interp(psi_new,psi,x_cam)    
    y_cam = np.interp(psi_new,psi,y_cam)   
    psi =psi_new
    norm = calc_norm(x_cam,y_cam,plotting=True)  
    #print psi
    #thk = F_TF(psi,controls_thk,Tte,plotting=1,spacing=None)         ### calculate thickness distribution
    
    thk = F_TF_bez(controls_thk_x,controls_thk_y,Tte,Oval_frac,plotting=1,spacing=None)


    Z_U = y_cam+thk*norm[:,1]	### apply thickness onto camber line for upper
    Z_L = y_cam-thk*norm[:,1]	### apply thickness onto camber line for lower
    X_U = x_cam+thk*norm[:,0]	### repeat upper for x rather than y
    X_L = x_cam-thk*norm[:,0]

    if plotting==1:		### optional plot out of complete blade
	#print 'Gamma:',Gamma,' Xi1:',Xi1,' Xi2:', Xi2
        plt.figure('Geometry Definition')
        plt.subplot(224)
        #plt.title('Blade Profile - F_Make')
        plt.axis('equal')
        plt.plot(x_cam,y_cam,'-c',label = 'Construct line')
        plt.plot(X_U,Z_U,'-b', label = 'Upper')
        plt.plot(X_L,Z_L,'-r', label = 'Lower')
        #for i in range(len(psi)):				### optional plot of applied thickness vectors
         #   plt.plot([psi[i],X_U[i]],[y_cam[i],Z_U[i]],'-y')
         #   plt.plot([psi[i],X_L[i]],[y_cam[i],Z_L[i]],'-y')
	#plt.legend(loc=0)
        #plt.show()
    X = np.append(X_L[::-1],X_U[1:])	### combine axial distributions upper and lower
    Z = np.append(Z_L[::-1],Z_U[1:])	### combine tangential distributions upper and lower
    
    #(X,Z) = Rescale(X,Z)		### perform rescale to ensure unit chord
    return(X,Z)
    
### SHAPE SPACE


def dist_vino(n, xst, xen, s1, s2):	### vinokur distribution 
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
    #print trans(0.00001), b, math.
    left = 0.00001
    delta = SciOpt.brentq(trans, left, 100)
    #delta = scipy.optimize.fmin_powell(trans, trans(0.0001), disp=False)
    u = 0.5*(1.0 + np.tanh(delta*(eta/nm1 - 0.5))/(np.tanh(delta*0.5)))
    s = u/(a + (1-a)*u)
    return xst + (xen-xst)*s     

    
def F_CF(psi): # function to generate class function of psi
    C = (psi**0.5)*(1.0-psi)
    return(C)

    
def F_SC(psi,i,Nm1): # function to generate shape component
    C = (psi**i)*((1.0-psi)**(Nm1-i))
    return(C)
    
def F_SF(psi,AIs,plotting): ### function generates Shape function based on reduced design parameters

    N = len(AIs)	### number of coefficients
    Nm1=N-1		### N-1
    npts = len(psi)	### number of points in distribution
    S = np.zeros((npts)) ### initilise shape space distribution

    #plotting = 1 ### optional plot out
    if plotting==1:
        plt.figure('Geometry Definition')
        plt.subplot(222)
        #plt.title('shape function')
        plt.xlabel('psi')
        plt.ylabel('Shape Function: S(psi)')

    for i in range(N):	### iterate over each contirbution
        Ai=AIs[i]
        Ki = binomial(i,Nm1)		### calculate bernstein polynomial
        SC = F_SC(psi,i,Nm1)		### calculate shape component
        S[:]=S[:]+Ai*Ki*SC[:]

        if plotting==1:### another optional plot out
            plt.plot(psi,Ai*Ki*SC[:],'-',label='component'+str(i))
        

    if plotting==1:### another optional plot out
        plt.plot(psi,S,'--k',label='total of components')
        plt.plot(psi,0.3*S[0]*(1.-psi)**18,'--r')
        
    S=S+0.3*S[0]*(1.-psi)**18
    if plotting==1:
            plt.plot(psi,S,'-k',label='total of components')
	#plt.legend(loc=0)
        #plt.show()  ### uncomment to plot imediatly     
    return(S)
    
def F_TF(psi,AIs,TE,plotting,spacing): ### function that returns thickness contributions due to parameters

    if spacing == None:
    	S=F_SF(psi,AIs,plotting)	### calculate shape function
    else:
    	S=F_SF(spacing,AIs,plotting)	### calculate shape function
    #plt.show()
    #S2 = np.interp(psi,psi2,S)
    C=F_CF(psi)		### calculate class function

    TF =S*C	+psi**0.5*TE/2.	### combine to form thickness funtion
    TFMAX=2.*TF.max()
    if plotting==1:	### optional plot out
        plt.figure('Geometry Definition')
        plt.subplot(223)
        #plt.title('thickness function')
        plt.xlabel('psi')
        plt.ylabel('thickness')
        plt.plot(psi,TF/TFMAX,'-r')
        plt.plot(psi,(psi**0.5*TE/2.)/TFMAX,'--r')
        plt.plot(psi,(S*C	+psi**0.5*TE/2)/TFMAX,'--b')
        #plt.plot(psi,-TF,'-b')
        plt.axis('equal')
    return(TF)
    
def F_TF_bez(controls_cam_x,controls_cam_y,TE,Oval_frac,plotting,spacing): ### function that returns thickness contributions due to parameters

    
    psi = dist_vino(200, 0, 1, 1./2000, 1./500)  
    X,S,cp = GEN_BZ_2D(controls_cam_x,controls_cam_y,plotting)
    
 
    
    S = np.interp(psi,X,S)
    #print S
    
    if plotting==1:	### optional plot out
        plt.figure('Geometry Definition')
        plt.subplot(222)
        #plt.title('thickness function')
        plt.plot(psi,S,'--m',label='total of components')
    
    S=S+Oval_frac*S[0]*(1.-psi)**18
    if plotting==1:	### optional plot out
        plt.figure('Geometry Definition')
        plt.subplot(222)
        #plt.title('thickness function')
        plt.xlabel('psi')
        plt.ylabel('SS')
        plt.plot(psi,S,'-m',label='total of components')
    #plt.show()
    #S2 = np.interp(psi,psi2,S)
    C=F_CF(psi)		### calculate class function

    TF =S*C	+psi*TE/2.	### combine to form thickness funtion
    TFMAX=2.*TF.max()
    if plotting==1:	### optional plot out
        plt.figure('Geometry Definition')
        plt.subplot(223)
        #plt.title('thickness function')
        plt.xlabel('psi')
        plt.ylabel('thickness')
        plt.plot(psi,TF/TFMAX,'-r')
        plt.plot(psi,(psi*TE/2.)/TFMAX,'--r')
        #plt.plot(psi,(S*C	+psi**0.5*TE/2)/TFMAX,'-b')
        #plt.plot(psi,-TF,'-b')
        plt.axis('equal')
        
        print 'Peak thickness:',2.*TF.max()*100,'% Cx'
        for i in range(len(TF)):
            if TF[i]==TF.max():
                print '@ ',100.*psi[i],' % camber line'
                break
            
        Bte = np.degrees(np.arctan((TF[-2]-TF[-1])/(psi[-1]-psi[-2])))
        print 'Boat tail:',Bte
    
    return(TF)
    
    
def calc_norm(x_cam,y_cam,plotting):
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
        #print A,portiondroop
	#plotting=1							### calculate normal 2 as 1/grad1
        if plotting==1:	### optional plot out
            plt.figure('Geometry Definition')
            plt.subplot(224)
            #plt.title('Camber line - showing normal vectors')
            plt.xlabel('psi')
            plt.ylabel('zeta')
            plt.plot(x_cam,y_cam,'-k',label = 'camber')
            plt.axis('equal')
            for i in range(len(x_cam))[::5]:
                plt.plot([x_cam[i],x_cam[i]+norm[i,0]*0.01],[y_cam[i],y_cam[i]+norm[i,1]*0.01],'k')	### plot normal vectors
                plt.plot([x_cam[i],x_cam[i]-norm[i,0]*0.01],[y_cam[i],y_cam[i]-norm[i,1]*0.01],'k')
	    #plt.show()
        return(norm)    
        
def F_Make(Xi1,Xi2,controls_cam,controls_thk_x,controls_thk_y,Tte,Oval_frac,plotting):        ### function that compiles a blade from reduced variables 


    ### turning distribution
    Xi1 = np.radians(Xi1)
    Xi2 = np.radians(Xi2)
    (x_cam,y_cam)=GEN_TURNING(Xi1,Xi2,controls_cam,plotting)      ### generate camber line from bezier curve
   

    Gamma = np.degrees(np.arctan((y_cam[-1])/x_cam[-1]))
    print 'Xi1:',Xi1,'Xi2:',Xi2,'Gamma:',Gamma

    s2=np.zeros((len(x_cam)))			### initialise camber length array
    for i in range(1,len(x_cam)):
            s2[i]=s2[i-1]+((x_cam[i]-x_cam[i-1])**2+(y_cam[i]-y_cam[i-1])**2)**0.5	
    psi = s2/s2.max()
    

    
    psi_new = dist_vino(200, 0, 1, 1./2000, 1./500)   
    
    x_cam = np.interp(psi_new,psi,x_cam)    
    y_cam = np.interp(psi_new,psi,y_cam)   
    psi =psi_new
    norm = calc_norm(x_cam,y_cam,plotting=True)  
    #print psi
    #thk = F_TF(psi,controls_thk,Tte,plotting=1,spacing=None)         ### calculate thickness distribution
    
    thk = F_TF_bez(controls_thk_x,controls_thk_y,Tte,Oval_frac,plotting=1,spacing=None)


    Z_U = y_cam+thk*norm[:,1]	### apply thickness onto camber line for upper
    Z_L = y_cam-thk*norm[:,1]	### apply thickness onto camber line for lower
    X_U = x_cam+thk*norm[:,0]	### repeat upper for x rather than y
    X_L = x_cam-thk*norm[:,0]

    if plotting==1:		### optional plot out of complete blade
	#print 'Gamma:',Gamma,' Xi1:',Xi1,' Xi2:', Xi2
        plt.figure('Geometry Definition')
        plt.subplot(224)
        #plt.title('Blade Profile - F_Make')
        plt.axis('equal')
        plt.plot(x_cam,y_cam,'-c',label = 'Construct line')
        plt.plot(X_U,Z_U,'-b', label = 'Upper')
        plt.plot(X_L,Z_L,'-r', label = 'Lower')
        #for i in range(len(psi)):				### optional plot of applied thickness vectors
         #   plt.plot([psi[i],X_U[i]],[y_cam[i],Z_U[i]],'-y')
         #   plt.plot([psi[i],X_L[i]],[y_cam[i],Z_L[i]],'-y')
	#plt.legend(loc=0)
        #plt.show()
    X = np.append(X_L[::-1],X_U[1:])	### combine axial distributions upper and lower
    Z = np.append(Z_L[::-1],Z_U[1:])	### combine tangential distributions upper and lower
    
    #(X,Z) = Rescale(X,Z)		### perform rescale to ensure unit chord
    return(X_L,Z_L,X_U,Z_U)


def Profile2(X1,X2,TKTE):		### function based on profgen.f (JDD) - GENERATES SIMPLE BLADE PROFILES


    #### THIS FUNCTION SHOULD BE REPLACED AT SOME POINT TO USE A BETTER BLADE SECTION GENERATOR
    ### EXCUSE THE WEIRD FORMULATION, THIS IS CONVERSION FROM FORTRAN

    #X2 = X2#*1.07

    ### USER HARD CODED VALUES BELOW - COULD BE TAKEN OUTSIDE FUNCTION CJC
    
    NPTS = 500	        ### points on aerofoil
    TKMAX = 0.15        ### max thickness
    XTMAX = 0.4	        ### max thickness location
    TKLE = 0.1*TKMAX	### leading edge thickness
    #TKTE = 0.015	### trailing edge thicknes  CHANGE WAKEDIST CODE IF THIS VALUE CHANGES
    XIND = 3.0	        ### thickness profile shape

    controls_cam = [-0.0,-0.5,-0.8]
    Rle = 0.05
    Tte = TKTE
    Beta_te = 4.0
    Oval_frac = 0.3
    #controls_cam = [0,0.2,0.0]
    controls_thk_x = [0.0,0.5,1.0]
    controls_thk_y = [(2.*Rle)**0.5*(1.-Oval_frac),0.25,np.tan(np.radians(Beta_te))+Tte/2.]

    (XSIN,YSIN,XPIN,YPIN) = F_Make(X1,X2,controls_cam,controls_thk_x,controls_thk_y,Tte,Oval_frac,plotting=True)

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

    XPIN = np.interp(np.linspace(0,Sle,500),Snew,Xnew)[::-1]
    YPIN = np.interp(np.linspace(0,Sle,500),Snew,Ynew)[::-1]
    XSIN = np.interp(np.linspace(Sle,Snew[-1],500),Snew,Xnew)
    YSIN = np.interp(np.linspace(Sle,Snew[-1],500),Snew,Ynew)

    XScale = max(max(XPIN),max(XSIN))-Xle
    XPIN = (XPIN-Xle)/XScale
    XSIN = (XSIN-Xle)/XScale
    XIN = (XIN-Xle)/XScale
			
    YPIN=(YPIN-yoffset)/XScale	### center on peak thickness - DEFAULTS STACKING TO PEAK THICKNESS
    YSIN=(YSIN-yoffset)/XScale	### center on peak thickness - DEFAULTS STACKING TO PEAK THICKNESS
    YIN = (YIN - yoffset)/XScale ### center on peak thickness - DEFAULTS STACKING TO PEAK THICKNESS

    
    if 1==1:			### OPTIONAL PLOT OUT OF THE SECTION PROFILE
        plt.figure(num = 'PROFILE GENERATOR')
        plt.axis('equal')
        plt.plot(XIN,YIN,'-xb',label = 'Camber')
        plt.plot(XPIN,YPIN,'-xr',label = 'PS')
        plt.plot(XSIN,YSIN,'-xg',label = 'SS')
	plt.legend(loc=0)


	
	
        #plt.show()

    delta = 0.5*abs(XSIN[-1]-XPIN[-1])		### RECALCULATE TRAILING EDGE THICKNESS
        
    return(XIN,YPIN,YSIN,XPIN,XSIN,delta)




def Profile(X1,X2,TKTE):		### function based on profgen.f (JDD) - GENERATES SIMPLE BLADE PROFILES


    #### THIS FUNCTION SHOULD BE REPLACED AT SOME POINT TO USE A BETTER BLADE SECTION GENERATOR
    ### EXCUSE THE WEIRD FORMULATION, THIS IS CONVERSION FROM FORTRAN

    #X2 = X2#*1.07

    ### USER HARD CODED VALUES BELOW - COULD BE TAKEN OUTSIDE FUNCTION CJC
    
    NPTS = 500	        ### points on aerofoil
    TKMAX = 0.15        ### max thickness
    XTMAX = 0.4	        ### max thickness location
    TKLE = 0.1*TKMAX	### leading edge thickness
    #TKTE = 0.015	### trailing edge thicknes  CHANGE WAKEDIST CODE IF THIS VALUE CHANGES
    XIND = 3.0	        ### thickness profile shape

    ### GENERATE TURNING DICTRIBUTION
    XIN = np.linspace(0,1,NPTS)
    SLOPE = np.tan(np.radians(np.linspace(0,1,NPTS)**1.0*(X2-X1)+X1)) ##edit power to move max camber position - DEFAULT LINEAR 
    YIN=np.zeros((NPTS))
    YSIN=np.zeros((NPTS))
    YPIN=np.zeros((NPTS))
    XSIN=np.zeros((NPTS))
    XPIN=np.zeros((NPTS))
    for i in range(1,NPTS):
    	DX=XIN[i]-XIN[i-1]
        YIN[i]=YIN[i-1]+DX*SLOPE[i]		### INTEGRATE SLOPE TO GET CAMBER LINE
        
    XDIF=XIN-XIN[0]		### CALCULATE AXIAL DISTANCE
    YDIF=YIN-YIN[0]		### CALCULATE TANGENTIAL DISTANCE
    STAGGER = XDIF/((XDIF*XDIF+YDIF*YDIF)**0.5)			### CALCULATE STAGGER ANGLE ATAN(DY/DX)

    ### GENERATE THICKNESS DISTRIBUTION BASED 
    POWER = np.log(0.5)/np.log(XTMAX)
    XTRANS=XIN**POWER
    TKLIN = TKLE+XIN**0.5*(TKTE-TKLE)	
    FX=abs((XTRANS-0.5)*2)
    TKIN = (TKLIN+TKMAX*(1.0-FX**XIND))/STAGGER

    ### ASSEMBLE Pressure and Suction surfaces
    for i in range(1,NPTS):
        if XIN[i]<=TKLE:
            YPIN[i]=YIN[i]+0.5*TKIN[i]*(XIN[i]/TKLE)**0.2
            YSIN[i]=YIN[i]-0.5*TKIN[i]*(XIN[i]/TKLE)**0.2
            XPIN[i]=XIN[i]
            XSIN[i]=XIN[i]
       # elif XIN[i]>TKLE and SLOPE[i]*SLOPE[1]>0:
        #    YPIN[i]=YIN[i]+0.5*TKIN[i]
        #    YSIN[i]=YIN[i]-0.5*TKIN[i]
        #    XPIN[i]=XIN[i]
         #   XSIN[i]=XIN[i]
        else:
            YPIN[i]=YIN[i]+0.5*np.cos(np.arctan(SLOPE[i]))*TKIN[i]
            YSIN[i]=YIN[i]-0.5*np.cos(np.arctan(SLOPE[i]))*TKIN[i]
            XPIN[i]=XIN[i]-0.5*np.sin(np.arctan(SLOPE[i]))*TKIN[i]
            XSIN[i]=XIN[i]+0.5*np.sin(np.arctan(SLOPE[i]))*TKIN[i]
        if XIN[i-1]<XTMAX and XIN[i]>=XTMAX:
            yoffset = YIN[i]

 

    Xnew =np.zeros((len(XPIN)+len(XSIN)-1))
    Xnew[:500]=XPIN[::-1]
    Xnew[499:]=XSIN[:]

    Ynew =np.zeros((len(YPIN)+len(YSIN)-1))
    Ynew[:500]=YPIN[::-1]
    Ynew[499:]=YSIN[:]

    Snew = np.zeros((len(Ynew)))
    for i in range(1,len(Xnew)):
	Snew[i]=Snew[i-1]+((Xnew[i]-Xnew[i-1])**2.+(Ynew[i]-Ynew[i-1])**2.)**0.5

    for i in range(1,len(Xnew)):
	if Xnew[i]==Xnew.min():
		Sle=Snew[i]
		Xle=Xnew[i]
		break

    XPIN = np.interp(np.linspace(0,Sle,500),Snew,Xnew)[::-1]
    YPIN = np.interp(np.linspace(0,Sle,500),Snew,Ynew)[::-1]

    XScale = max(max(XPIN),max(XSIN))-Xle
    XPIN = (XPIN-Xle)/XScale
    XSIN = (XSIN-Xle)/XScale
    XIN = (XIN-Xle)/XScale
			
    YPIN=(YPIN-yoffset)/XScale	### center on peak thickness - DEFAULTS STACKING TO PEAK THICKNESS
    YSIN=(YSIN-yoffset)/XScale	### center on peak thickness - DEFAULTS STACKING TO PEAK THICKNESS
    YIN = (YIN - yoffset)/XScale ### center on peak thickness - DEFAULTS STACKING TO PEAK THICKNESS

    

    

    
    if 1==1:			### OPTIONAL PLOT OUT OF THE SECTION PROFILE
        plt.figure(num = 'PROFILE GENERATOR')
        plt.axis('equal')
        plt.plot(XIN,YIN,'-xb',label = 'Camber')
        plt.plot(XPIN,YPIN,'-xr',label = 'PS')
        plt.plot(XSIN,YSIN,'-xg',label = 'SS')
	plt.legend(loc=0)

	plt.figure(num = 'Thickness Distribution')
	plt.plot(XIN,TKLIN,'-k')
	plt.plot(XIN,TKIN,'-r')
	
        plt.show()

    delta = 0.5*abs(XSIN[NPTS-1]-XPIN[NPTS-1])		### RECALCULATE TRAILING EDGE THICKNESS
        
    return(XIN,YPIN,YSIN,XPIN,XSIN,delta)



def set_default(g,Cp,gam,Rgas,visc):		### FUNCTION THAT SETS DEFAULT VARIABLES WITHIN TURBOSTREAM - LOADED FROM ts_tstream_default.txt
    # SET APPLICATION VARIABLES - SOLVER WIDE VARIABLES
    for name in ts_tstream_default.av:
        val = ts_tstream_default.av[name]
	#print name, val
        if type(val) == type(1):
            g.set_av(name, ts_tstream_type.int, val)
        else:
            g.set_av(name, ts_tstream_type.float, val)
	    
    g.set_av('cp', ts_tstream_type.float, Cp)	   
    g.set_av('ga', ts_tstream_type.float, gam)
    g.set_av('viscosity', ts_tstream_type.float, visc)  
     
    # SET BLOCK VARIABLES - BLOCK ONLY VARIABLES
    for bid in g.get_block_ids():
        for name in ts_tstream_default.bv:
            val = ts_tstream_default.bv[name]
	    #print name, val
            if type(val) == type(1):
                g.set_bv(name, ts_tstream_type.int, bid, val)
            else:
                g.set_bv(name, ts_tstream_type.float, bid, val)
	

    ### SET VISCOSITY MODEL DEFAULT - SOLVER DOESNT DO IT PROPERLY
    trans_dyn_vis=(0.84*(2.0+4.6))*visc
    print 'voscosity:',visc
    for bid in g.get_block_ids():
        		b=g.get_block(bid)
			g.set_bp("trans_dyn_vis",ts_tstream_type.float,bid,np.zeros([b.nk,b.nj,b.ni],np.float32)+visc)


def linearspac(delx,Nv):	### A LINEAR SPACING FUNCTION 

	### RETURNS A LINEARLY VARYING SPACING FUNCTION THAT GOES FROM 0->1 WITH A CONSTANT EXPANSION RATIO
	### distribution spans 0-1 with fixed spacing @0

	exprat = 1.0
	for it in range(100):
		spac = np.ones((Nv-1))
		for i in range(Nv-1):
			spac[i]=spac[i-1]*exprat
		spac=spac*delx
		dx=sum(spac)
		if dx<1.0:exprat=exprat*1.01
		else:exprat=exprat*0.99

	#spac=((spac-spac.min())/(spac.max()-spac.min()))
	dist = np.zeros((Nv))
	for i in range(1,Nv):
		dist[i]=dist[i-1]+spac[i-1]
	dist=((dist-dist.min())/(dist.max()-dist.min()))

	#print 'expansion ratio used:',exprat
	if exprat>1.3: print 'WARNING- EXPANSION RATIO HIGH:',exprat
	return(dist)



###############################
### BELOW FUNCTION IS LARGE ###
### IT ADDS PATCHES LINKING ###
### ALL THE BLOCKS TOGETHER ###
###############################

### NATHAN STOPPED USING THIS APPARENTLY



def Poisson_smth(x,r,t,Nsmth):
	### perform poisson smoothing

	#plt.figure()
	#plt.plot(x,t,'-k')
	#plt.plot(x.T,t.T,'-k')

	for it in range(Nsmth):
			xnew=copy.deepcopy(x)
			tnew=copy.deepcopy(t)
			rnew=copy.deepcopy(r)

			xnew[1:-1,1:-1] = 0.25*(x[:-2,1:-1]+x[2:,1:-1]+x[1:-1,:-2]+x[1:-1,2:])
			tnew[1:-1,1:-1] = 0.25*(t[:-2,1:-1]+t[2:,1:-1]+t[1:-1,:-2]+t[1:-1,2:])
			rnew[1:-1,1:-1] = 0.25*(rnew[:-2,1:-1]+rnew[2:,1:-1]+rnew[1:-1,:-2]+rnew[1:-1,2:])
		
			x=copy.deepcopy(xnew)
			r=copy.deepcopy(rnew)

			t=copy.deepcopy(tnew)

			
			#xgrid1[i,NI_up-5:NI_up+5,:]=copy.deepcopy(xgridnew1[i,NI_up-5:NI_up+5,:])
			#rgrid1[i,NI_up-5:NI_up+5,:]=copy.deepcopy(rgridnew1[i,NI_up-5:NI_up+5,:])
			#tgrid1[i,NI_up-5:NI_up+5,:]=copy.deepcopy(tgridnew1[i,NI_up-5:NI_up+5,:])

	#plt.plot(x,t,'-r')
	#plt.plot(x.T,t.T,'-r')
	#plt.show()

	return(xnew,rnew,tnew)


#####################################
### MAIN ROUTINE TO GENERATE MESH ###
### ROWS: LIST OF ROW CLASS FILES ###
### VXIN: INLET AXIAL VELOCITY    ###
### ROIN: INLET DENSITY 	  ###
### PSHUBIN : INLET HUB PS        ###
#####################################

def MESH(ROWS,VXin,ROin,PSHUBin,gam,Rgas,Cp,visc):	### MAIN MESHING PROGRAM

	print 'STARTING MAIN MESHING ROUTINE...'

	rownumber = 0			### IDENTIFY NUMBER OF ROWS TO BE MESHED
	nrows = len(ROWS)
	print 'NUMBER OF ROWS:',nrows
	g=ts_tstream_grid.TstreamGrid()	### initialise grid


	################################
	### CREATE STREAMLINE SPLINES -- ENABLES CURVED ENDWALLS
	### CAN CAUSE PROBLEMS IN SHORTER MACHINES

	XH = np.zeros((nrows+2))	### INTIIALISE ARRAY FOR HUB AXIAL DISTANCE
	XC = np.zeros((nrows+2))	### INITIALISE ARRAY FOR CAS AXIAL DISTANCE
	RH= np.zeros((nrows+2))		### INITIALISE ARRAY FOR HUB RADIUS
	RC = np.zeros((nrows+2))	### INITIALISE ARRAY FOR CAS RADIUS

	### ITERATE THROUGH ROWS SETTING RADIUS AND AXIAL LOCATIONS AT CENTER POINTS OF ROW
	### FIRST ENTRY IS AT 0 AXIAL DISTANCE
	for i in range(len(ROWS)):
		row=ROWS[i]
		if i==0:
			XH[i+1]=XH[i]+row.Cx*(0.5+row.upfrac[0])	### FIRST ROW LOCATED AT HALF A CHORD AND AN UPSTREAM SPACING
		else:
			rowm1=ROWS[i-1]
			XH[i+1]=XH[i]+row.Cx*(0.5+row.downfrac[0])+rowm1.Cx*(0.5+row.upfrac[0]) ### ITERATE SPACINGS TO FINS SUBSEQUENT LOCATIONS
		RH[i+1]=row.RADIUS[0]			### SET RADII HUB
		RC[i+1]=row.RADIUS[-1]			### SET RADII CAS
	XH[-1]=XH[-2]+row.Cx*(0.5+row.downfrac[0])	### FINAL ENTRY IS HALF A CHORD AND A DOWNSTREAM SPACING 
	XC[:]=XH[:]					### SET ALL CASING AXIAL LOCATIONS THE SAME AS THE HUBS
	RH[0]=RH[1]					### SET INLET RADIUS TO SAME AS FIRST ROW
	RC[0]=RC[1]					### SET INLET RADIUS TO SAME AS FIRST ROW
	RH[-1]=RH[-2]					### SET EXIT RADIUS TO SAME AS LAST ROW
	RC[-1]=RC[-2]					### SET EXIT RADIUS TO SAME AS LAST ROW

	if nrows>1:
		h_spline = sciiint.splrep(XH, RH,k=3)	### generate splines for hub - DEFAULT ORDER K=3
		c_spline = sciiint.splrep(XC, RC,k=3)	### generate splines for cas - DEFAULT ORDER K=3
		XSL = np.linspace(XH[0],XH[-1],400)
		HSL = sciiint.splev( XSL,h_spline)
		CSL = sciiint.splev( XSL,c_spline)

	else: 
		XSL = XH
		HSL = RH
		CSL = RC

	if 1==1:
		plt.figure(num = 'ANNULUS SPLINES')
		plt.title('annulus spline fits')
		plt.hold('on')
		plt.plot(XH,RH,'ob',label = 'hub control')
		plt.plot(XC,RC,'or',label = 'cas control')
		plt.plot(XSL,HSL,'-b',label = 'hub spline')
		plt.plot(XSL,CSL,'-r',label = 'cas spline')
		plt.axis('equal')
		#plt.show()
	### END OF ENDWALL SPLINES


	### CALCULATE MASSFLOW OF WHOLE MACHINE
	A=np.pi*(RC[0]**2-RH[0]**2)		### INLET AXIAL AREA
	MDOT=A*VXin*ROin			### MASS FLOW
	P0IN = PSHUBin+0.5*ROin*VXin**2		### ESTIMATE OF INLET TOTAL PRESSURE - COULD BE IMPROVED
	T0IN = P0IN/(Rgas*ROin)			### ESTIMATE OF INLET TOTAL TEMPERATURE - COULD BE IMPROVED
	YAWIN = ROWS[0].X1[0]			### INLET FLOW ANGLE - STARTING WITH STATOR

	print 'MDOT:', MDOT, 'P0IN: ', P0IN, 'T0IN: ', T0IN, 'ROin: ', ROin, 'VXin:', VXin, 'YAWIN:', YAWIN

	#######################################################
	### ESTIMATE PERFORMANCE AND GENERATE INITIAL GUESS ###
	#######################################################
	p0I=np.zeros((len(ROWS)))	### INTIALISE ARRAY OF TOTAL PRESSURE @ ROW INLET
	p0E=np.zeros((len(ROWS)))	### INTIALISE ARRAY OF TOTAL PRESSURE @ ROW INLET
	p0rI=np.zeros((len(ROWS)))	### INTIALISE ARRAY OF TOTAL PRESSURE @ ROW INLET
	p0rE=np.zeros((len(ROWS)))	### INTIALISE ARRAY OF TOTAL PRESSURE @ ROW INLET
	pI=np.zeros((len(ROWS)))	### INITIALISE ARRAY OF STATIC PRESSURE
	pE=np.zeros((len(ROWS)))	### INITIALISE ARRAY OF STATIC PRESSURE
	t0I=np.zeros((len(ROWS)))	### INITIALISE ARRAY OF TOTAL TEMPERATURE
	t0E=np.zeros((len(ROWS)))	### INITIALISE ARRAY OF TOTAL TEMPERATURE
	t0rI=np.zeros((len(ROWS)))	### INITIALISE ARRAY OF TOTAL TEMPERATURE
	t0rE=np.zeros((len(ROWS)))	### INITIALISE ARRAY OF TOTAL TEMPERATURE
	vI=np.zeros((len(ROWS)))	### INITIALISE ARRAY OF INLET VELOCITY
	vE=np.zeros((len(ROWS)))	### INITIALISE ARRAY OF EXIT VELOCITY
	vtE=np.zeros((len(ROWS)))	### INITIALISE ARRAY OF EXIT VELOCITY
	vtI=np.zeros((len(ROWS)))	### INITIALISE ARRAY OF EXIT VELOCITY
	roI=np.zeros((len(ROWS)))	### INITIALISE ARRAY OF INLET DENSITY

	print 'Cp: ', Cp, 'Gam: ', gam, 'R: ', Rgas	

	fdev = 0.97

	### ITERATE THROUGH ROWS -note cjc currently uses hub angles -  midspan would be better
	### THIS IS SETTING INITIAL GUESS
	for i in range(len(ROWS)):
		row=ROWS[i]
		if i == 0:
			###first row
			p0I[i]=P0IN
			t0I[i]=T0IN
			vtI[i] = VXin*np.tan(np.radians(YAWIN))

		else:
			p0I[i]=p0E[i-1]	### carried from previous exit
			t0I[i]=t0E[i-1]	### carried from previous exit
			vtI[i]=vtE[i-1]

		rtip =row.RADIUS[-1]
		rhub=row.RADIUS[0]
		rmid=(rtip+rhub)/2.				### CALCULATE MIDSPAN RADIUS
		A=np.pi*(rtip**2-rhub**2)			### CALCULATE AREA
		U = row.RPM*(rmid*2.*np.pi)/60.			### BLADE SPEED OF ROW

		vtrI = vtI[i]-U

		t0rI[i]=t0I[i]-(vtI[i]**2.)/(2.*Cp)+(vtrI**2.)/(2.*Cp)
		p0rI[i]=p0I[i]*(t0rI[i]/t0I[i])**(gam/(gam-1))
		

		### INLET CALCULATION

			
		ro = p0rI[i]/(Rgas*t0rI[i])				### CALCULATE DENSITY
		c = (gam*Rgas*t0rI[i])**0.5			#BL## CALCULATE SPEED OF SOUND
		for it in range(10):
			vx = MDOT/(A*ro)			### CALCULATE AXIAL VELOCITY
			v=(vx**2.+vtrI**2)**0.5 ### CALCULATE RELATIVE INLET ABSOLUTE VELOCITY
			t=t0rI[i]*(1+(gam-1)*((v/c)**2)/2)**(-1.0)	### CALCULATE STATIC PRESSURE
			c = (gam*Rgas*t)**0.5			#BL## RECLAULATE SOS
			pI[i] = p0rI[i]*(1.+(gam-1.)*((v/c)**2.)/2)**(-gam/(gam-1))	### CALCULATE STATIC PRESSURE			
			ro = pI[i]/(Rgas*t)			### RECALCULATE DENSITY
		vI[i]=v
		roI[i] = ro
		print 'vx:', vx, 'ro:', ro

		#pI[i]=pI[i]-0.5*ro*vt**2+0.5*ro*vt_abs**2
		yawpredict = np.degrees(np.arctan(vtrI/vx))
		print 'predicted inlet angle:',yawpredict,row.X1[0]
		if abs(yawpredict-row.X1[0])>5.0:print 'warning -traingles dont make sense'

		### TO and P0 drop calculation
		v2est = vx/np.cos(np.radians(row.X2[0]*fdev))

		p0rE[i]=p0rI[i]-0.03*(0.5*ro*v2est**2) ### assumed 5% loss coefficient - P0 DROP
		t0rE[i]=t0rI[i]


		### EXIT CALCULATION
		ro = p0rE[i]/(Rgas*t0rE[i])				### CALCULATE DENSITY
		c = (gam*Rgas*t0rE[i])**0.5			#BL## CALCULATE SPEED OF SOUND
		for it in range(10):
			vx = MDOT/(A*ro)			### CALCULATE AXIAL VELOCITY
			v=vx/np.cos(np.radians(row.X2[0]*fdev))   ### CALCULATE RELATIVE INLET ABSOLUTE VELOCITY
			vtr = vx*np.tan(np.radians(row.X2[0]*fdev))
			vt = vtr+U
			t=t0rE[i]*(1+(gam-1)*((v/c)**2)/2)**(-1.0)	### CALCULATE STATIC PRESSURE
			c = (gam*Rgas*t)**0.5			#BL## RECLAULATE SOS
			pE[i] = p0rE[i]*(1.+(gam-1.)*((v/c)**2.)/2.)**(-gam/(gam-1.))	### CALCULATE STATIC PRESSURE			
			ro = pE[i]/(Rgas*t)			### RECALCULATE DENSITY

		vE[i]=v
		
		t0E[i]=t0rE[i]-(vtr**2)/(2.*Cp)+(vt**2)/(2.*Cp)
		p0E[i]=p0rE[i]*(t0E[i]/t0rE[i])**(gam/(gam-1.))
		
		Re = ro*v*row.Cx/visc
		print 'Reynolds number:',Re
		

		vtE[i]=vt
		print 'velocities: vx:',vx,'vexit:',vE[i],'vinlet:',vI[i],'vtexit,rel:',vtr,'vtexit:',vt




	if 1==1: # OPTIONAL PLOT OF INITIAL GUESS
		plt.figure(num = 'INITIAL GUESS')
		plt.subplot(311)
		plt.hold('on')
		plt.plot(p0I,'-ob',label='P0,in')
		plt.plot(p0E,'-oc',label = 'P0,out')
		plt.plot(p0rI,'--xb',label='P0r,in')
		plt.plot(p0rE,'--xc',label = 'P0r,out')
		plt.plot(pI,'-or',label = 'Ps,in')
		plt.plot(pE,'-om',label = 'Ps,out')
		plt.legend(loc=0)
		#plt.xlabel('row')
		plt.ylabel('Pressure(Pa)')
		plt.subplot(312)	
		plt.plot(t0I,'-ob')
		plt.plot(t0E,'-oc')
		#plt.xlabel('row')
		plt.ylabel('T0(K)')
		plt.subplot(313)	
		plt.plot(vI,'-xb')
		plt.plot(vE,'-xc')
		plt.xlabel('row')
		plt.ylabel('Velocity')
		#plt.show()

	#################################
	#### CODE FOR GENERATING MESH ###
	#################################

	xshift = 0.0	### cumulative shift for future rows AXIAL OFFSET
	irow = 0	### set counter
	bidsup=[]	### intilise list for upstream block ids
	bidsmid1=[]	### initialise list for mid1 blocks
	bidsmid2=[]	### initialise list for mid2 blocks
	bidsdown=[]	### initialise list for downstream blocks
	bidssrdblock = {}	### initialise shroud block list

	### shroud blocks added by nathan C
        no_blocks = 1+S1.Nfins  
	if S1.Nfins == 0:
		no_blocks = 2           
        for n in range(no_blocks+4):
                bidssrdblock[n] = []


	### ITERATE THROUGH EACH ROW GENERATING GRID
	for row in ROWS:
		print 'GENERATING ROW #:',rownumber,'name:',row.name,'RPM:',row.RPM

		### CALCULATE ROW FEATURES
		Nsections = len(row.X1)			# NUMBER OF SECTIONS
		rhub = min(row.RADIUS) 			### find hub radius
		if rhub<=0.0:print 'ERROR ####  HUB RADIUS IS <0'
		rtip = max(row.RADIUS) 			### find cas radius
		H=rtip-rhub				### span
		hubtiprat = rhub/rtip			### hubtip ratio

		### DEFINE AXIAL LOCATIONS OF FEATURES
		xblock0 = -0.2*row.Cx ### extra distance on inlet block - added by NATHAN C
		xstart = 0.0+xshift				### ROW START
		xle = xstart + row.upfrac[0]*row.Cx		### LEADING EDGE
		xsrdupstart = xle-row.upfrac[1]*row.Cx		### UPSTREAM SHROUD START
		xsrdupend = xle-row.upfrac[2]*row.Cx		### UPSTREAM SHROUD END
		xte = xle+row.Cx				### TRAILING EDGE
		xsrddnstart = xte+row.downfrac[1]*row.Cx	### DOWNSTREAM SHROUD START
		xsrddnend = xte+row.downfrac[2]*row.Cx		### DOWNSTREAM SHROUD END
		xend = xte+row.downfrac[0]*row.Cx		### DOMAIN END

		#################################################
		### MESH COUNTS AND SPACINGS	#################
		### USER UNPUTS BELOW           #################
		#################################################

		### MESH COUNTS
		N1 = 20		### UPSTREAM OF SHROUD
		N2 = 15		### INSIDE SHROUD
		N3 = 20 	### UPSTREAM OF LE     
		N4 = 70		### PASSAGE
		N5 = 25 	### DOWNSTREAM OF PASSAGE  
		N6 = 20		### DOWNSTREAM OF SHROUD
		N7 = 12		### DEPTH OF SHROUD
		N8 = 40		### AXIAL LENGTH OF SHROUD
		N9 = 20		### OVER KNIFE TIP --- TO REMOVE
	
		NJ = 50		### SPANWISE 
		NK = 50		### PITCHWISE
	
		#### GRID SPACINGS
		LEspacing = 0.001*row.Cx 	### leading edge clustering as fraction of AXIAL chord
		TEspacing = 0.001*row.Cx 	### trailing edge clustering as fraction of AXIAL chord
		EWspacing = 0.003 		### EW clustering as fraction of span
		BLspacing = 0.003 		### BL clustering as fraction of pitch
		#EWspacing = 0.005 		### EW clustering as fraction of span
		#BLspacing = 0.005 		### BL clustering as fraction of pitch
		SRDspacing = 0.05*row.Cx*(row.upfrac[1]-row.upfrac[2]) ### shroud clustering as fraction of shroud width


		##########################
		### end of user input ####
		##########################
		### calculate total i direction mesh for main block
		NI = N1+N2-1+N3-1+N4-1+N5-1+N2-1+N6-1		### SUMMING ALONG BLOCKS IGNORING OVERLAPPING POINTS

		### GENERATE SECTIONS
		X = np.zeros((Nsections,NI,NK))			### INIT. ARRAY FOR AXIAL NODES
		R = np.zeros((Nsections,NI,NK))			### INIT. ARRAY FOR RADIAL NODES
		RT = np.zeros((Nsections,NI,NK))		### INIT ARRAY FOR TANGENTIAL NODES

		colors = ['-b','-g','-r','-c','-m','-k','-y']	### SET COLORS FOR PLOTTING MANY SECTIONS AT ONCE

		### ITERATE  ACROSS SECTIONS GENERATING PROFILES
		for N in range(Nsections):		### ITERATE THROUGH SECTIONS
			print 'generating section:',N


			### CALCULATE THE PITCH/NUMBER OF BLADES
			if row.Nblades>=3:
				S=2.*np.pi*row.RADIUS[N]/row.Nblades			### PITCH IS FOUND DIVIDING 2*PI*R BY NBLADES
				print 'nblades specified, using pitch to chord:',S/row.Cx
			if row.Nblades<3:						### if number of blades is below 3 treat number as hub pitch to chord
				S=row.Nblades*row.Cx*row.RADIUS[N]/row.RADIUS[0]	### PITCH IS FOUND FROM HUB P2C -cjc could be better if mispan was used
				row.Nblades = 2.*(np.pi*row.RADIUS[N]/S)
				print 'pitch to chord specified, using nblades=',row.Nblades

			LEspacing = 0.5*S*BLspacing#/row.Cx
			TEspacing = 0.5*S*BLspacing#/row.Cx
			print 'spacings:',LEspacing,TEspacing


			### GENERATE PITCHWISE DISTRIBUTION TO BE USED THROUGHOUT
			pitchdist = vinokur_dist(0,1,NK,BLspacing,BLspacing) 	### fixed symetric spacing used throughout domain - could be not fixed -would be better
			pitchdist_U = vinokur_dist(0,1,NK,BLspacing*3.0,BLspacing*3.0) 	### fixed symetric spacing used throughout domain - could be not fixed -would be better

			### GENERATE AEROFOIL NEW
			[xin,ypin,ysin,xpin,xsin,delta] = Profile2(row.X1[N],row.X2[N],row.thkTE)	### load seciton geometry
			xin = xin*row.Cx+xle
			xpin = xpin*row.Cx+xle
			xsin = xsin*row.Cx+xle
			ypin = ypin*row.Cx
			ysin = ysin*row.Cx
			delta = delta*row.Cx
			yle = ypin[0]
			yte = ypin[-1]


			### GENERATE AXIAL DISTRIBUTION

			### UPSTREAM OF SHROUD	
			dist = linearspac(SRDspacing/(xsrdupstart-xstart),N1)		### use a linear spacing upstream of shroud
			dist = (abs(dist[::-1]-1))*(xsrdupstart-xstart)+xstart
			for k in range(NK):
				X[N,:N1,k]=dist[:]

			if rownumber == 0:			### ADDITIONAL SPACE IN DOMAIN INLET BLOCK - ADDED NATHAN C
				dist = linearspac(SRDspacing/(xsrdupstart-xblock0),N1)
				dist = (abs(dist[::-1]-1))*(xsrdupstart-xblock0)+xblock0
				for k in range(NK):
					X[N,:N1,k]=dist[:]

			### INSIDE SHROUD
			print('INSIDE UPSTREAM SHROUD')
			dist = vinokur_dist(xsrdupstart,xsrdupend,N2,SRDspacing,SRDspacing) 	### use a vinokur spacing within the shroud block
			for k in range(NK):
				X[N,N1-1:N1-1+N2,k]=dist[:]


			### BETWEEN SHROUD AND LE
			print('BETWEEN SHROUD AND LEADING EDGE')
			dist = vinokur_dist(xsrdupend,xle,N3,SRDspacing,LEspacing) 		### use a vinokur spacing between the shroud and the leading edge
			for k in range(NK):
				X[N,N1-1+N2-1:N1-1+N2-1+N3,k]=dist[:]		

			for i in range(N1+N2-1+N3-1):
				if i <N1+N2-1-1:
					axialfrac = 0.0
				else:
					axialfrac = (float((i-(N1-1+N2-1)))/(N3-1))**1
				#axialfrac = max(float(i-(N1+N2-1))/(N1+N2-1+N3-1.-1),0)
				ylow= yle+(X[N,i,0]-xle)*np.tan(np.radians(row.X1[N]))

				localdist = pitchdist_U+(pitchdist-pitchdist_U)*axialfrac
				RT[N,i,:]=ylow+S*localdist#pitchdist

			### BLADE PASSAGE		--- use two vinokur spacings within the blade passage -one for each surface
			if row.X2[N]>=0:		### assumption of exit flow angle for stator
				print('UPPER PASSAGE DIST')
				dist_top = vinokur_dist(xle,(xte),N4,LEspacing,TEspacing)
				print('LOWER PASSAGE DIST')
				dist_bottom = vinokur_dist(xle,(xte-2.*delta),N4,LEspacing,TEspacing)
			if row.X2[N] < 0:		### if blade row is flipped so are PS and SS in the generator
				print('UPPER PASSAGE DIST')
				dist_top = vinokur_dist(xle,(xte-2.*delta),N4,LEspacing,TEspacing)
				print('LOWER PASSAGE DIST')
				dist_bottom = vinokur_dist(xle,(xte),N4,LEspacing,TEspacing)
			y_top = np.interp(dist_top,xsin,ysin)+S		### interpolate upper surface - with pitch added on
			y_bottom = np.interp(dist_bottom,xpin,ypin)	### interpolate lower surface - no pitch added


			### iterate through the passage generating X RT grid
			for j in range(N4):
				length = math.sqrt((y_top[j]-y_bottom[j])**2+(dist_top[j]-dist_bottom[j])**2)
				cos_theta = (y_top[j]-y_bottom[j])/length
				sin_theta = (dist_top[j]-dist_bottom[j])/length
				for z in range(NK):
					X[N,N1-1+N2-1+N3-1+j,z] = dist_bottom[j]+sin_theta*length*pitchdist[z]
					RT[N,N1-1+N2-1+N3-1+j,z] = y_bottom[j]+cos_theta*length*pitchdist[z]


        		pitchdist_U=pitchdist
			### TE---> SHROUD  
			print('TE-> SHROUD DIST - UPPER')
			#print(xsin[-1],xsrddnstart,N5,TEspacing,SRDspacing)

			top_x = vinokur_dist(xsin[-1],xsrddnstart,N5,TEspacing,SRDspacing)	#vinokur between upper te and shroud
			print('TE->SHROUD DIST - LOWER')
			#print(xpin[-1],xsrddnstart,N5,TEspacing,SRDspacing)
			bottom_x = vinokur_dist(xpin[-1],xsrddnstart,N5,TEspacing,SRDspacing)	### vinokur between lower te and shroud
			top_y = ysin[-1]+S+(top_x-xsin[-1])*np.tan(np.radians(row.X2[N]*0.95))
			bottom_y = ypin[-1]+(bottom_x-xpin[-1])*np.tan(np.radians(row.X2[N]*0.95))

			for j in range(N5):
    
				length = math.sqrt((top_y[j]-bottom_y[j])**2+(top_x[j]-bottom_x[j])**2)
				cos_theta = (top_y[j]-bottom_y[j])/length
				sin_theta = (top_x[j]-bottom_x[j])/length
				axialfrac = (bottom_x[j]-bottom_x[0])/row.Cx
				if axialfrac<0.05:
					axialfrac = 0.0
				elif axialfrac>0.25:axialfrac = 1.0
				else:axialfrac = (axialfrac - 0.05)/0.2
				if axialfrac >1.0:axialfrac=1.0
				localdist = pitchdist+(pitchdist_U-pitchdist)*axialfrac
    
				for z in range(NK):
        
					X[N,N1-1+N2-1+N3-1+N4-1+j,z] = bottom_x[j]+sin_theta*length*localdist[z]
					RT[N,N1-1+N2-1+N3-1+N4-1+j,z] = bottom_y[j]+cos_theta*length*localdist[z]

			### INSIDE SHROUD	- vinokur distribution within shroud
			print('INSIDE DOWNSTREAM SHROUD')
			dist = vinokur_dist(xsrddnstart,xsrddnend,N2,SRDspacing,SRDspacing)

			for k in range(NK):
				X[N,N1-1+N2-1+N3-1+N4-1+N5-1:N1-1+N2-1+N3-1+N4-1+N5-1+N2,k]=dist[:]

			### DOWNSTREAM OF SHROUD	- constant expansion outside of shroud
			print('DOWNSTREAM OF SHROUD')
			dist = linearspac(SRDspacing/(xend-xsrddnend),N6)#[::-1]
			dist = dist*(xend-xsrddnend)+xsrddnend

			for k in range(NK):
				X[N,N1-1+N2-1+N3-1+N4-1+N5-1+N2-1:N1-1+N2-1+N3-1+N4-1+N5-1+N2-1+N6,k]=dist[:]

			### DOWNSTREAM GRID FOLLOWS EXIT ANGLE --- set RT of everything downstream of trailing edge
			for i in range(N1-1+N2-1+N3-1+N4-1+N5-1,N1-1+N2-1+N3-1+N4-1+N5-1+N2-1+N6):
				ylow = ypin[-1]+(X[N,i,0]-xpin[-1])*np.tan(np.radians(row.X2[N]*0.95))# 7% deviation
				RT[N,i,:]=ylow+length*pitchdist_U

			if N == 0:			### if section is hub initialise arrays
				NW = 0
				X2 = np.zeros((Nsections,N5-1+N2-1+N6,NK))
				R2 = np.zeros((Nsections,N5-1+N2-1+N6,NK))
				RT2 = np.zeros((Nsections,N5-1+N2-1+N6,NK))


			if 0==1:
				ile = N1-1+N2-1+N3
				iles = ile -5
				ilee=ile
				Nsmth = 20
				### smooth leading edge region
				for Ismth in range(N3-5):
					if Ismth == 0:
						Nsmth = 2
						(X[N,iles-Ismth:ilee+Ismth,Ismth:],Rn,RT[N,iles-Ismth:ilee+Ismth,Ismth:]) = Poisson_smth(X[N,iles-Ismth:ilee+Ismth,Ismth:],RT[N,iles-Ismth:ilee+Ismth,Ismth:],RT[N,iles-Ismth:ilee+Ismth,Ismth:],Nsmth)
					else:
						Nsmth=2
						(X[N,iles-Ismth:ilee+Ismth,Ismth:-Ismth],Rn,RT[N,iles-Ismth:ilee+Ismth,Ismth:-Ismth]) = Poisson_smth(X[N,iles-Ismth:ilee+Ismth,Ismth:-Ismth],RT[N,iles-Ismth:ilee+Ismth,Ismth:-Ismth],RT[N,iles-Ismth:ilee+Ismth,Ismth:-Ismth],Nsmth)

			if 1==0:
				ile = N1-1+N2-1+N3-1+N4
				iles = ile -10
				ilee=ile+10
				Nsmth = 20
				### smooth T edge region
				for Ismth in range(15):
					if Ismth == 0:
						Nsmth = 15
						(X[N,iles-Ismth:ilee+Ismth,Ismth:],Rn,RT[N,iles-Ismth:ilee+Ismth,Ismth:]) = Poisson_smth(X[N,iles-Ismth:ilee+Ismth,Ismth:],RT[N,iles-Ismth:ilee+Ismth,Ismth:],RT[N,iles-Ismth:ilee+Ismth,Ismth:],Nsmth)
					else:
						Nsmth=5
						(X[N,iles-Ismth:ilee+Ismth,Ismth:-Ismth],Rn,RT[N,iles-Ismth:ilee+Ismth,Ismth:-Ismth]) = Poisson_smth(X[N,iles-Ismth:ilee+Ismth,Ismth:-Ismth],RT[N,iles-Ismth:ilee+Ismth,Ismth:-Ismth],RT[N,iles-Ismth:ilee+Ismth,Ismth:-Ismth],Nsmth)


			X2[N,:,:NK] = X[N,N1-1+N2-1+N3-1+N4-1:N1-1+N2-1+N3-1+N4-1+N5-1+N2-1+N6,:NK]
			RT2[N,:,:NK] = RT[N,N1-1+N2-1+N3-1+N4-1:N1-1+N2-1+N3-1+N4-1+N5-1+N2-1+N6,:NK]



			### WAKE AREA BEHIND THICK TRAILING EDGE
			if xsin[-1] != xpin[-1] or ysin[-1] != ypin[-1]:	### check for thick trailing edge
				
			### TE---> SHROUD WAKE 
				top_x_wake = bottom_x
				bottom_x_wake = top_x
				top_y_wake = bottom_y+S
				bottom_y_wake = top_y

				if row.X2[N] < 0:		### switch if rotor
        				top_y_wake = bottom_y
        				bottom_y_wake = top_y-S

				if N == 0:			### IF HUB SECTION CALCULATE THE NUMBER OF POINTS NEEDED
					scal = (np.pi*(row.RADIUS[0]+row.RADIUS[1])/row.Nblades)/(row.thkTE*row.Cx) #CALCULATE RATIO OF PITCH TO TE THICKNESS

					edge_BL = S*BLspacing*0.5
					T_TE = (row.thkTE*row.Cx)

					NW = min(int(1.5/(scal*BLspacing)),30)	### ESTIMATE NUMBER OF CELS REQUIRED
					NW = min(int(0.5*T_TE/edge_BL),30)	### ESTIMATE NUMBER OF CELS REQUIRED
					print 'TE CELL COUNTS: TE THICKNESS:',T_TE*1000.,', BLcell size:',edge_BL*1000., ', Ncells=',NW, T_TE/edge_BL

					if NW <= 10:
						print 'LESS THAN spaced TE POINTS NEEDED -USING linear'
						NW = int(0.8*round(T_TE/edge_BL)+1)
						print 'replacing Nte:',NW
						#pitchwake = np.linspace(0,1,NW)

					if NW <= 5:
						print 'LESS THAN 5 TE POINTS NEEDED -USING 5'
						NW = 5
						pitchwake = np.linspace(0,1,NW)
					else:
						#print('wakedist')
						pitchwake = vinokur_dist(0,1,NW,edge_BL/T_TE,edge_BL/T_TE)
						#print('end')
				x_val2 = np.zeros([Nsections,N5-1+N2-1+N6,NW])
				y_val2 = np.zeros([Nsections,N5-1+N2-1+N6,NW])

				for j in range(N5):
    
				        length = math.sqrt((top_y_wake[j]-bottom_y_wake[j])**2+(top_x_wake[j]-bottom_x_wake[j])**2)
				        cos_theta = (top_y_wake[j]-bottom_y_wake[j])/length
				        sin_theta = (top_x_wake[j]-bottom_x_wake[j])/length
    
				        for z in range(NW):
        
					        x_val2[N,j,z] = bottom_x_wake[j]+sin_theta*length*pitchwake[z]
        					y_val2[N,j,z] = bottom_y_wake[j]+cos_theta*length*pitchwake[z]

				### INSIDE SHROUD WAKE
				dist = vinokur_dist(xsrddnstart,xsrddnend,N2,SRDspacing,SRDspacing)

				for k in range(NW):
				        x_val2[N,N5-1:N5-1+N2,k]=dist[:]

				### DOWNSTREAM OF SHROUD WAKE
				dist = linearspac(SRDspacing/(xend-xsrddnend),N6)#[::-1]
				dist = dist*(xend-xsrddnend)+xsrddnend

				for k in range(NW):
	        			x_val2[N,N5-1+N2-1:N5-1+N2-1+N6,k]=dist[:]
    
				### DOWNSTREAM GRID FOLLOWS EXIT ANGLE WAKE
				for i in range(N5-1,N5-1+N2-1+N6):
				        ylow = ysin[-1]+S+(x_val2[N,i,0]-xsin[-1])*np.tan(np.radians(row.X2[N]*0.95))### 7% deviation
				        y_val2[N,i,:]=ylow+length*pitchwake

				if N == 0:
					X2 = np.zeros((Nsections,N5-1+N2-1+N6,NK-1+NW))
					R2 = np.zeros((Nsections,N5-1+N2-1+N6,NK-1+NW))
					RT2 = np.zeros((Nsections,N5-1+N2-1+N6,NK-1+NW))

				X2[N,:,0:NK] = X[N,N1-1+N2-1+N3-1+N4-1:N1-1+N2-1+N3-1+N4-1+N5-1+N2-1+N6,:]
				RT2[N,:,0:NK] = RT[N,N1-1+N2-1+N3-1+N4-1:N1-1+N2-1+N3-1+N4-1+N5-1+N2-1+N6,:]

				X2[N,:,NK-1:NK-1+NW] = x_val2[N,:,:]
				RT2[N,:,NK-1:NK-1+NW] = y_val2[N,:,:]

				if row.X2[N] < 0:
					y_val2[N,N5-1:N5-1+N2-1+N6,:] = y_val2[N,N5-1:N5-1+N2-1+N6,:] - S

					X2[N,:,NW-1:NW-1+NK] = X[N,N1-1+N2-1+N3-1+N4-1:N1-1+N2-1+N3-1+N4-1+N5-1+N2-1+N6,:]
					RT2[N,:,NW-1:NW-1+NK] = RT[N,N1-1+N2-1+N3-1+N4-1:N1-1+N2-1+N3-1+N4-1+N5-1+N2-1+N6,:]

					X2[N,:,0:NW] = x_val2[N,:,:]
					RT2[N,:,0:NW] = y_val2[N,:,:]

			### Plot profile with mesh
			if 0==1:
			   if N == 0 and row.X2[N] > 0:
				plt.figure(70)
				plt.subplot(221)
				plt.plot(dist_top,y_top-S,'o')
				plt.plot(dist_bottom,y_bottom,'o')
				plt.xticks([])
				plt.yticks([])
				plt.xlim(0.003,0.014)
				plt.ylim(-0.003,0.008)
				#plt.axis('equal')

				plt.subplot(222)
				plt.plot(X[N,N1-1+N2-1+N3-1:N1-1+N2-1+N3-1+N4,:(NK-1)/2-8],RT[N,N1-1+N2-1+N3-1:N1-1+N2-1+N3-1+N4,:(NK-1)/2-8],'-g') 
				plt.plot(X[N,N1-1+N2-1+N3-1:N1-1+N2-1+N3-1+N4,:(NK-1)/2-8].T,RT[N,N1-1+N2-1+N3-1:N1-1+N2-1+N3-1+N4,:(NK-1)/2-8].T,'-g')
				plt.plot(X[N,N1-1+N2-1+N3-1:N1-1+N2-1+N3-1+N4,(NK-1)/2+8:],RT[N,N1-1+N2-1+N3-1:N1-1+N2-1+N3-1+N4,(NK-1)/2+8:]-S,'-g') 
				plt.plot(X[N,N1-1+N2-1+N3-1:N1-1+N2-1+N3-1+N4,(NK-1)/2+8:].T,RT[N,N1-1+N2-1+N3-1:N1-1+N2-1+N3-1+N4,(NK-1)/2+8:].T-S,'-g')
				if xsin[-1] != xpin[-1] or ysin[-1] != ypin[-1]:
					plt.plot(X2[N,:N5-21,:(NK-1)/2-8],RT2[N,:N5-21,:(NK-1)/2-8],'-m')
					plt.plot(X2[N,:N5-21,:(NK-1)/2-8].T,RT2[N,:N5-21,:(NK-1)/2-8].T,'-m')
					plt.plot(X2[N,:N5-21,(NK-1)/2+8:],RT2[N,:N5-21,(NK-1)/2+8:]-S,'-m')
					plt.plot(X2[N,:N5-21,(NK-1)/2+8:].T,RT2[N,:N5-21,(NK-1)/2+8:].T-S,'-m')
				plt.xticks([])
				plt.yticks([])
				plt.xlim(0.003,0.014)
				plt.ylim(-0.003,0.008)
				#plt.axis('equal')

			### SET RADIUS TO SPANWISE FRACTION OF SECTION 
			R[N,:,:]=(row.RADIUS[N]-row.RADIUS[0])/(row.RADIUS[-1]-row.RADIUS[0])
			R2[N,:,:]=(row.RADIUS[N]-row.RADIUS[0])/(row.RADIUS[-1]-row.RADIUS[0])




			### PLOT BLADE TO BLADE GRIDS
	
######
			if 1==1 and N == 0:
				plt.figure(num = 'Blade to Blade Mesh')
				#plt.subplot(200+rownumber+10*nrows+1)
				label = 'inlet:'+str(row.X1[N])+' exit:'+str(row.X2[N])
				plt.plot(X[N,:N1,:],RT[N,:N1,:],'-b')
				plt.plot(X[N,N1-1:N1-1+N2-1+N3-1+N4,:],RT[N,N1-1:N1-1+N2-1+N3-1+N4,:],'-g')
				plt.plot(X[N,:N1,:].T,RT[N,:N1,:].T,'-b') 
				plt.plot(X[N,N1-1:N1-1+N2-1+N3-1+N4,:].T,RT[N,N1-1:N1-1+N2-1+N3-1+N4,:].T,'-g')
				# plot new thickess wake
				if xsin[-1] != xpin[-1] or ysin[-1] != ypin[-1]:
					plt.plot(X2[N,:N5-1+N2,:],RT2[N,:N5-1+N2,:],'-m')
					plt.plot(X2[N,N5-1+N2-1:N5-1+N2-1+N6,:],RT2[N,N5-1+N2-1:N5-1+N2-1+N6,:],'-c')
					plt.plot(X2[N,:N5-1+N2,:].T,RT2[N,:N5-1+N2,:].T,'-m')
					plt.plot(X2[N,N5-1+N2-1:N5-1+N2-1+N6,:].T,RT2[N,N5-1+N2-1:N5-1+N2-1+N6,:].T,'-c')
					plt.axis('equal')
					plt.xlabel('axial distance(m)')
					plt.ylabel('tangential distance(m)')
					plt.title(row.name+'B2B')

######			
			#plt.show()

			for i in range(1,len(X[N,:,0])):
				if X[N,i,0]-X[N,i-1,0]==0.0:print 'error - duplicate'

		#####################################################################################################
		### HAVE NOW GENERATED 2D MESH ON ALL SECTIONS DESIGNATED - NOW NEED TO INTERPOLATE BETWEEN THEM  ###
		#####################################################################################################
		
		### GENERATE BLOCKS

		### generate spanwise distribution used throughout
		spandist = vinokur_dist(0,1,NJ,EWspacing,EWspacing) 

		### UPSTREAM OF SHROUD BLOCK
		b = ts_tstream_type.TstreamBlock()		### INIT BLOCK
		### ASSIGN BLOCK NUMBER
		if rownumber ==0:
			b.bid=0
		else:
			b.bid = g.get_block_ids()[-1]+1
		###APPEND LIST OF UPSTREAM BLOCKS - USED FOR PATCHING
		bidsup.append(b.bid)
		b.ni = N1	### NUMBER AXIAL
		b.nj = NJ	### NUMBER RADIAL
		b.nk = NK	### NUMBER TANGENTIAL
		g.add_block(b)	### ADD BLOCK TO GRID
		x= np.zeros((b.nk,b.nj,b.ni),np.float32)	### INITIALISE X ARRAY
		r= np.zeros((b.nk,b.nj,b.ni),np.float32)	### INITIALISE X ARRAY
		rt= np.zeros((b.nk,b.nj,b.ni),np.float32)	### INITIALISE X ARRAY


		### PERFORM INTERPOLATION FROM NSECTIONS TO NJ POINTS
		if Nsections>3.: ### for lots of sections use cubic
			for i in range(b.ni):
				for k in range(b.nk):

					x_spline = sciiint.splrep(R[:,i,k], X[:,i,k],k=3)	### generate splines between B2B planes
					rt_spline = sciiint.splrep(R[:,i,k], RT[:,i,k],k=3)	
					x[k,:,i] = sciiint.splev( spandist,x_spline)
					rt[k,:,i] = sciiint.splev( spandist,rt_spline)
					rcasing=np.interp(x[k,-1,i],XSL,CSL)			### INTERPOLATE FROM ANNULUS LINES
					rfoot=np.interp(x[k,0,i],XSL,HSL)			### INTERPOLATE FROM ANNULUS LINES
					r[k,:,i]=spandist*(rcasing-rfoot)+rfoot
		elif Nsections==3:	### for 3 sections use quadratic interp
			for i in range(b.ni):
				for k in range(b.nk):

					x_spline = sciiint.splrep(R[:,i,k], X[:,i,k],k=2)	### generate splines between B2B planes
					rt_spline = sciiint.splrep(R[:,i,k], RT[:,i,k],k=2)
					x[k,:,i] = sciiint.splev( spandist,x_spline)
					rt[k,:,i] = sciiint.splev( spandist,rt_spline)
					rcasing=np.interp(x[k,-1,i],XSL,CSL)			### INTERPOLATE FROM ANNULUS LINES
					rfoot=np.interp(x[k,0,i],XSL,HSL)			### INTERPOLATE FROM ANNULUS LINES
					r[k,:,i]=spandist*(rcasing-rfoot)+rfoot
		else:
			for i in range(b.ni):
				for k in range(b.nk):
					x[k,:,i] = X[0,i,k]+(X[-1,i,k]-X[0,i,k])*(spandist-spandist.min())/(spandist.max()-spandist.min())
					rt[k,:,i] =RT[0,i,k]+(RT[-1,i,k]-RT[0,i,k])*(spandist-spandist.min())/(spandist.max()-spandist.min())
					rcasing=np.interp(x[k,-1,i],XSL,CSL)			### INTERPOLATE FROM ANNULUS LINES
					rfoot=np.interp(x[k,0,i],XSL,HSL)			### INTERPOLATE FROM ANNULUS LINES
					r[k,:,i]=spandist*(rcasing-rfoot)+rfoot


		g.set_bp("x", ts_tstream_type.float, b.bid, x)			### SET X SPACINGS
		g.set_bp("r", ts_tstream_type.float, b.bid, r)			### SET R SPACINGS
		g.set_bp("rt", ts_tstream_type.float, b.bid, rt)		### SET RT SPACINGS



		### SET ROTATION OF BLOCK
		g.set_bv('rpm',ts_tstream_type.float,b.bid,row.RPM)	### BLOCK ALWAYS ROTATES WITH ATTACHED BLADE
		g.set_bv('rpmi1',ts_tstream_type.float,b.bid,0.0)		### AXIAL INLET FACE DOESNT ROTATE
		g.set_bv('rpmi2',ts_tstream_type.float,b.bid,row.RPM)
		if row.shroud_hub==True:
			if rownumber==0:
				g.set_bv('rpmj1',ts_tstream_type.float,b.bid,0.0)#ROWS[rownumber+1].RPM)	### set it so it rotates with subsequent row
			else:	g.set_bv('rpmj1',ts_tstream_type.float,b.bid,ROWS[rownumber-1].RPM)     ### set it so it rotates with previous rotor
		else:
			if rownumber==0:
				g.set_bv('rpmj1',ts_tstream_type.float,b.bid,0.0)#ROWS[rownumber+1].RPM)	### set it so it rotates with subsequent row
			else:	g.set_bv('rpmj1',ts_tstream_type.float,b.bid,row.RPM)#ROWS[rownumber-1].RPM)     ### set it so it rotates with previous rotor
		if row.shroud_cas==True:
			if rownumber==0:
				g.set_bv('rpmj2',ts_tstream_type.float,b.bid,0.0)#ROWS[rownumber+1].RPM)	### set it so it rotates with subsequent row
			else:	g.set_bv('rpmj2',ts_tstream_type.float,b.bid,ROWS[rownumber-1].RPM) ### ste it so it rotates with previous rotor
		else:
			if rownumber==0:
				g.set_bv('rpmj2',ts_tstream_type.float,b.bid,row.RPM)### casing spings with row
			else:	g.set_bv('rpmj2',ts_tstream_type.float,b.bid,row.RPM)#ROWS[rownumber-1].RPM) ### ste it so it rotates with previous rotor
			
		g.set_bv('rpmk1',ts_tstream_type.float,b.bid,row.RPM)	### K SURFACES SPIN WITH BLOCK
		g.set_bv('rpmk2',ts_tstream_type.float,b.bid,row.RPM)

		g.set_bv('xllim',ts_tstream_type.float,b.bid,0.03*S)
		g.set_bv('dampin_mul',ts_tstream_type.float,b.bid,1.0)	

		### SET BLOCK INITIAL GUESS
		g.set_bv('vgridin',ts_tstream_type.float,b.bid,vI[rownumber])
		g.set_bv('vgridout',ts_tstream_type.float,b.bid,vI[rownumber])
		g.set_bv('pstatin',ts_tstream_type.float,b.bid,pI[rownumber])
		g.set_bv('pstatout',ts_tstream_type.float,b.bid,pI[rownumber])
		g.set_bv('tstagin',ts_tstream_type.float,b.bid,t0I[rownumber])
		g.set_bv('tstagout',ts_tstream_type.float,b.bid,t0I[rownumber])
		if row.Nblades>3:
			g.set_bv('nblade',ts_tstream_type.int,b.bid,row.Nblades)
			g.set_bv('fblade',ts_tstream_type.float,b.bid,float(row.Nblades))
		else:
			g.set_bv('nblade',ts_tstream_type.int,b.bid,1)
			g.set_bv('fblade',ts_tstream_type.float,b.bid,1.0)




		### PASSAGE BLOCK  - Shroud to Trailing edge
		### ALL SAME AS BEFORE
		b = ts_tstream_type.TstreamBlock()
		b.bid = g.get_block_ids()[-1]+1
		b.ni = N2-1+N3-1+N4
		b.nj = NJ
		b.nk = NK
		g.add_block(b)
		x= np.zeros((b.nk,b.nj,b.ni),np.float32)
		r= np.zeros((b.nk,b.nj,b.ni),np.float32)
		rt= np.zeros((b.nk,b.nj,b.ni),np.float32)
		bidsmid1.append(b.bid)

		if Nsections>3.: ### for lots of sections use cubic
			for i in range(b.ni):
				for k in range(b.nk):

					x_spline = sciiint.splrep(R[:,N1-1+i,k], X[:,N1-1+i,k],k=3)	### generate splines between B2B planes
					rt_spline = sciiint.splrep(R[:,N1-1+i,k], RT[:,N1-1+i,k],k=3)	
					x[k,:,i] = sciiint.splev( spandist,x_spline)
					rt[k,:,i] = sciiint.splev( spandist,rt_spline)
					rcasing=np.interp(x[k,-1,i],XSL,CSL)
					rfoot=np.interp(x[k,0,i],XSL,HSL)
					r[k,:,i]=spandist*(rcasing-rfoot)+rfoot
		elif Nsections==3:	### for 3 sections use quadratic interp
			for i in range(b.ni):
				for k in range(b.nk):

					x_spline = sciiint.splrep(R[:,N1-1+i,k], X[:,N1-1+i,k],k=2)	### generate splines between B2B planes
					rt_spline = sciiint.splrep(R[:,N1-1+i,k], RT[:,N1-1+i,k],k=2)
					x[k,:,i] = sciiint.splev( spandist,x_spline)
					rt[k,:,i] = sciiint.splev( spandist,rt_spline)
					rcasing=np.interp(x[k,-1,i],XSL,CSL)
					rfoot=np.interp(x[k,0,i],XSL,HSL)
					r[k,:,i]=spandist*(rcasing-rfoot)+rfoot
		else:
			for i in range(b.ni):
				for k in range(b.nk):
					x[k,:,i] = X[0,N1-1+i,k]+(X[-1,N1-1+i,k]-X[0,N1-1+i,k])*(spandist-spandist.min())/(spandist.max()-spandist.min())
					rt[k,:,i] =RT[0,N1-1+i,k]+(RT[-1,N1-1+i,k]-RT[0,N1-1+i,k])*(spandist-spandist.min())/(spandist.max()-spandist.min())
					rcasing=np.interp(x[k,-1,i],XSL,CSL)
					rfoot=np.interp(x[k,0,i],XSL,HSL)
					r[k,:,i]=spandist*(rcasing-rfoot)+rfoot


		g.set_bp("x", ts_tstream_type.float, b.bid, x)
		g.set_bp("r", ts_tstream_type.float, b.bid, r)
		g.set_bp("rt", ts_tstream_type.float, b.bid, rt)



		### SET ROTATION OF BLOCK
		g.set_bv('rpm',ts_tstream_type.float,b.bid,row.RPM)	### BLOCK ALWAYS ROTATES WITH ATTACHED BLADE
		g.set_bv('rpmi1',ts_tstream_type.float,b.bid,row.RPM)	### hub spins with row 
		g.set_bv('rpmi2',ts_tstream_type.float,b.bid,row.RPM)   ### casing spings with row
		
		g.set_bv('rpmj1',ts_tstream_type.float,b.bid,row.RPM)	### ALL CASES IT SPINS WITH ROW
		g.set_bv('rpmj2',ts_tstream_type.float,b.bid,row.RPM)   ### ALL CASES IT SPINGS WITH ROW
	
		g.set_bv('rpmk1',ts_tstream_type.float,b.bid,row.RPM)	### K SURFACES SPIN WITH BLOCK
		g.set_bv('rpmk2',ts_tstream_type.float,b.bid,row.RPM)

		g.set_bv('xllim',ts_tstream_type.float,b.bid,0.03*S)
		g.set_bv('dampin_mul',ts_tstream_type.float,b.bid,1.0)	
		### SET INITIAL GUESS
		g.set_bv('vgridin',ts_tstream_type.float,b.bid,vI[rownumber])
		g.set_bv('vgridout',ts_tstream_type.float,b.bid,vE[rownumber])
		g.set_bv('pstatin',ts_tstream_type.float,b.bid,pI[rownumber])
		g.set_bv('pstatout',ts_tstream_type.float,b.bid,pE[rownumber])
		g.set_bv('tstagin',ts_tstream_type.float,b.bid,t0I[rownumber])
		g.set_bv('tstagout',ts_tstream_type.float,b.bid,t0E[rownumber])
		if row.Nblades>3:
			g.set_bv('nblade',ts_tstream_type.int,b.bid,row.Nblades)
			g.set_bv('fblade',ts_tstream_type.float,b.bid,float(row.Nblades))
		else:
			g.set_bv('nblade',ts_tstream_type.int,b.bid,1)
			g.set_bv('fblade',ts_tstream_type.float,b.bid,1.0)

		#####################################################

		### NOT SURE WHAT THIS IS DOING - PART OF NATHANS CODE
		k_old = [0]*(b.nk)   # changes nk in upstream block from NK to NK-1+NW so it corresponds to outlet shroud width
		for i in range(0,b.nk):
		        k_old[i] = i
		k_new = np.linspace(0,b.nk-1,NK-1+NW)
		rt_temp = np.zeros((NK-1+NW,b.nj,b.ni))

		for i in range(b.ni):
		        for j in range(b.nj):
		                rt_temp[:,j,i] = sciiint.interp1d(k_old, rt[k_old,j,i], kind='cubic')(k_new)
		
		rtpas2 = rt_temp
		xpas0 = x
		rpas0 = r
		rtpas0 = rt

		### end of nathans part



		### PASSAGE BLOCK 2 - TE to end of shroud
		b = ts_tstream_type.TstreamBlock()
		b.bid = g.get_block_ids()[-1]+1
		b.ni = N5-1+N2
		b.nj = NJ
		b.nk = NK-1+NW
		g.add_block(b)
		x= np.zeros((b.nk,b.nj,b.ni),np.float32)
		r= np.zeros((b.nk,b.nj,b.ni),np.float32)
		rt= np.zeros((b.nk,b.nj,b.ni),np.float32)
		bidsmid2.append(b.bid)

		if Nsections>3.: ### for lots of sections use cubic
			for i in range(b.ni):
				for k in range(b.nk):

					x_spline = sciiint.splrep(R2[:,i,k], X2[:,i,k],k=3)	### generate splines between B2B planes
					rt_spline = sciiint.splrep(R2[:,i,k], RT2[:,i,k],k=3)	
					x[k,:,i] = sciiint.splev( spandist,x_spline)
					rt[k,:,i] = sciiint.splev( spandist,rt_spline)
					rcasing=np.interp(x[k,-1,i],XSL,CSL)
					rfoot=np.interp(x[k,0,i],XSL,HSL)
					r[k,:,i]=spandist*(rcasing-rfoot)+rfoot
		elif Nsections==3:	### for 3 sections use quadratic interp
			for i in range(b.ni):
				for k in range(b.nk):

					x_spline = sciiint.splrep(R2[:,i,k], X2[:,i,k],k=2)	### generate splines between B2B planes
					rt_spline = sciiint.splrep(R2[:,i,k], RT2[:,i,k],k=2)
					x[k,:,i] = sciiint.splev( spandist,x_spline)
					rt[k,:,i] = sciiint.splev( spandist,rt_spline)
					rcasing=np.interp(x[k,-1,i],XSL,CSL)
					rfoot=np.interp(x[k,0,i],XSL,HSL)
					r[k,:,i]=spandist*(rcasing-rfoot)+rfoot
		else:
			for i in range(b.ni):
				for k in range(b.nk):
					x[k,:,i] = X2[0,i,k]+(X2[-1,i,k]-X2[0,i,k])*(spandist-spandist.min())/(spandist.max()-spandist.min())
					rt[k,:,i] =RT2[0,i,k]+(RT2[-1,i,k]-RT2[0,i,k])*(spandist-spandist.min())/(spandist.max()-spandist.min())
					rcasing=np.interp(x[k,-1,i],XSL,CSL)
					rfoot=np.interp(x[k,0,i],XSL,HSL)
					r[k,:,i]=spandist*(rcasing-rfoot)+rfoot



		g.set_bp("x", ts_tstream_type.float, b.bid, x)
		g.set_bp("r", ts_tstream_type.float, b.bid, r)
		g.set_bp("rt", ts_tstream_type.float, b.bid, rt)

		### SET ROTATION OF BLOCK
		g.set_bv('rpm',ts_tstream_type.float,b.bid,row.RPM)	### BLOCK ALWAYS ROTATES WITH ATTACHED BLADE
		g.set_bv('rpmi1',ts_tstream_type.float,b.bid,row.RPM)	### hub spins with row 
		g.set_bv('rpmi2',ts_tstream_type.float,b.bid,row.RPM)   ### casing spings with row
		g.set_bv('rpmj1',ts_tstream_type.float,b.bid,row.RPM)	### hub spins with row 
		g.set_bv('rpmj2',ts_tstream_type.float,b.bid,row.RPM)   ### casing spings with row
	
		g.set_bv('rpmk1',ts_tstream_type.float,b.bid,row.RPM)	### K SURFACES SPIN WITH BLOCK
		g.set_bv('rpmk2',ts_tstream_type.float,b.bid,row.RPM)

		g.set_bv('xllim',ts_tstream_type.float,b.bid,0.03*S)
		g.set_bv('dampin_mul',ts_tstream_type.float,b.bid,0.03)		
		g.set_bv('vgridin',ts_tstream_type.float,b.bid,vE[rownumber])
		g.set_bv('vgridout',ts_tstream_type.float,b.bid,vE[rownumber])
		g.set_bv('pstatin',ts_tstream_type.float,b.bid,pE[rownumber])
		g.set_bv('pstatout',ts_tstream_type.float,b.bid,pE[rownumber])
		g.set_bv('tstagin',ts_tstream_type.float,b.bid,t0E[rownumber])
		g.set_bv('tstagout',ts_tstream_type.float,b.bid,t0E[rownumber])
		g.set_bv('fmgrid',ts_tstream_type.float,b.bid,0.0)
		g.set_bv('poisson_fmgrid',ts_tstream_type.float,b.bid,0.0)
		if row.Nblades>3:
			g.set_bv('nblade',ts_tstream_type.int,b.bid,row.Nblades)
			g.set_bv('fblade',ts_tstream_type.float,b.bid,float(row.Nblades))
		else:
			g.set_bv('nblade',ts_tstream_type.int,b.bid,1)
			g.set_bv('fblade',ts_tstream_type.float,b.bid,1.0)


		#### NATHANS PART
		k_old = [0]*(b.nk)   # changes nk in downstream block from NK-1+NW to NK so it corresponds to inlet shroud width
		for i in range(0,b.nk):
		        k_old[i] = i
		k_new = np.linspace(0,b.nk-1,NK)
		rt_temp = np.zeros((NK,b.nj,b.ni))

		for i in range(b.ni):
		        for j in range(b.nj):
		                rt_temp[:,j,i] = sciiint.interp1d(k_old, rt[k_old,j,i], kind='cubic')(k_new)

		x[N,0,:] = X2[N,:b.ni,0]
		rtpas = rt_temp
		rtpas1 = rt
		xpas = x
		rpas = r
		### NATHANS PART


		### DOWNSTREAM BLOCK
		Nsofar = N5-1+N2-1
		b = ts_tstream_type.TstreamBlock()
		b.bid = g.get_block_ids()[-1]+1
		b.ni = N6
		b.nj = NJ
		b.nk = NK-1+NW
		g.add_block(b)
		x= np.zeros((b.nk,b.nj,b.ni),np.float32)
		r= np.zeros((b.nk,b.nj,b.ni),np.float32)
		rt= np.zeros((b.nk,b.nj,b.ni),np.float32)
		bidsdown.append(b.bid)
		if Nsections>3.: ### for lots of sections use cubic
			for i in range(b.ni):
				for k in range(b.nk):

					x_spline = sciiint.splrep(R2[:,Nsofar+i,k], X2[:,Nsofar+i,k],k=3)	### generate splines between B2B planes
					rt_spline = sciiint.splrep(R2[:,Nsofar+i,k], RT2[:,Nsofar+i,k],k=3)	
					x[k,:,i] = sciiint.splev( spandist,x_spline)
					rt[k,:,i] = sciiint.splev( spandist,rt_spline)
					rcasing=np.interp(x[k,-1,i],XSL,CSL)
					rfoot=np.interp(x[k,0,i],XSL,HSL)
					r[k,:,i]=spandist*(rcasing-rfoot)+rfoot
		elif Nsections==3:	### for 3 sections use quadratic interp
			for i in range(b.ni):
				for k in range(b.nk):

					x_spline = sciiint.splrep(R2[:,Nsofar+i,k], X2[:,Nsofar+i,k],k=2)	### generate splines between B2B planes
					rt_spline = sciiint.splrep(R2[:,Nsofar+i,k], RT2[:,Nsofar+i,k],k=2)
					x[k,:,i] = sciiint.splev( spandist,x_spline)
					rt[k,:,i] = sciiint.splev( spandist,rt_spline)
					rcasing=np.interp(x[k,-1,i],XSL,CSL)
					rfoot=np.interp(x[k,0,i],XSL,HSL)
					r[k,:,i]=spandist*(rcasing-rfoot)+rfoot
		else:
			for i in range(b.ni):
				for k in range(b.nk):
					x[k,:,i] = X2[0,Nsofar+i,k]+(X2[-1,Nsofar+i,k]-X2[0,Nsofar+i,k])*(spandist-spandist.min())/(spandist.max()-spandist.min())
					rt[k,:,i] =RT2[0,Nsofar+i,k]+(RT2[-1,Nsofar+i,k]-RT2[0,Nsofar+i,k])*(spandist-spandist.min())/(spandist.max()-spandist.min())
					rcasing=np.interp(x[k,-1,i],XSL,CSL)
					rfoot=np.interp(x[k,0,i],XSL,HSL)
					r[k,:,i]=spandist*(rcasing-rfoot)+rfoot

		g.set_bp("x", ts_tstream_type.float, b.bid, x)
		g.set_bp("r", ts_tstream_type.float, b.bid, r)
		g.set_bp("rt", ts_tstream_type.float, b.bid, rt)
		### SET ROTATION OF BLOCK
		g.set_bv('rpm',ts_tstream_type.float,b.bid,row.RPM)	### BLOCK ALWAYS ROTATES WITH ATTACHED BLADE
		g.set_bv('rpmi1',ts_tstream_type.float,b.bid,row.RPM)		### AXIAL INLET FACE DOESNT ROTATE
		g.set_bv('rpmi2',ts_tstream_type.float,b.bid,0.0)
		if row.shroud_hub==True:
			if rownumber==0:
				g.set_bv('rpmj1',ts_tstream_type.float,b.bid,ROWS[rownumber+1].RPM)	### set it so it rotates with subsequent row
			else:	g.set_bv('rpmj1',ts_tstream_type.float,b.bid,ROWS[rownumber-1].RPM) ### ste it so it rotates with previous rotor
		else:g.set_bv('rpmj1',ts_tstream_type.float,b.bid,row.RPM)### casing spings with row
		if row.shroud_cas==True:
			if rownumber==0:
				g.set_bv('rpmj2',ts_tstream_type.float,b.bid,ROWS[rownumber+1].RPM)	### set it so it rotates with subsequent row
			else:	g.set_bv('rpmj2',ts_tstream_type.float,b.bid,ROWS[rownumber-1].RPM) ### ste it so it rotates with previous rotor
		else:g.set_bv('rpmj2',ts_tstream_type.float,b.bid,row.RPM)### casing spings with row
	
		g.set_bv('rpmk1',ts_tstream_type.float,b.bid,row.RPM)	### K SURFACES SPIN WITH BLOCK
		g.set_bv('rpmk2',ts_tstream_type.float,b.bid,row.RPM)

		g.set_bv('xllim',ts_tstream_type.float,b.bid,0.03*S)
		g.set_bv('dampin_mul',ts_tstream_type.float,b.bid,0.03)	
		g.set_bv('vgridin',ts_tstream_type.float,b.bid,vE[rownumber])
		g.set_bv('vgridout',ts_tstream_type.float,b.bid,vE[rownumber])
		g.set_bv('pstatin',ts_tstream_type.float,b.bid,pE[rownumber])
		g.set_bv('pstatout',ts_tstream_type.float,b.bid,pE[rownumber])
		g.set_bv('tstagin',ts_tstream_type.float,b.bid,t0E[rownumber])
		g.set_bv('tstagout',ts_tstream_type.float,b.bid,t0E[rownumber])
		if row.Nblades>3:
			g.set_bv('nblade',ts_tstream_type.int,b.bid,row.Nblades)
			g.set_bv('fblade',ts_tstream_type.float,b.bid,float(row.Nblades))
		else:
			g.set_bv('nblade',ts_tstream_type.int,b.bid,1)
			g.set_bv('fblade',ts_tstream_type.float,b.bid,1.0)


		##########################
		#### ADD PATCHING ########
		##########################

		### PERIODIC UPSTREAM OF BLADE
		bid = bidsup[-1]		### just added upstream block
		bid0 = bid
		bid1 = bid
		b0 = g.get_block(bid0)
		b1 = g.get_block(bid1)
		p0 = ts_tstream_type.TstreamPatch()
		p0.bid = bid0
		p0.nxbid = bid1		### links to itself
		p0.ist = 0		### covers whole axial direction
		p0.ien = b0.ni
		p0.jst = 0		### covers whole radial dirction
		p0.jen = b0.nj
		p0.kst = 0
		p0.ken = 1
		p0.idir = 0		### i ->i
		p0.jdir = 1		### j->j
		p0.kdir = 2		### k->k
		p0.kind = 5	
		p1 = ts_tstream_type.TstreamPatch()
		p1.bid = bid1
		p1.nxbid = bid0
		p1.ist = 0
		p1.ien = b1.ni
		p1.jst = 0
		p1.jen = b1.nj
		p1.kst = b1.nk-1
		p1.ken = b1.nk
		p1.idir = 0
		p1.jdir = 1
		p1.kdir = 2
		p1.kind = 5

		p0.pid = g.add_patch(bid0, p0)
		p1.pid = g.add_patch(bid1, p1)
		p1.nxpid = p0.pid
		p0.nxpid = p1.pid

		### PERIODIC UPSTREAM OF BLADE 2
		bid0 = bidsmid1[-1]
		bid1 = bidsmid1[-1]
		b0 = g.get_block(bid0)
		b1 = g.get_block(bid1)

		p0 = ts_tstream_type.TstreamPatch()
		p0.bid = bid0
		p0.nxbid = bid1
		p0.ist = 0
		p0.ien = N2+N3-1
		p0.jst = 0
		p0.jen = b0.nj
		p0.kst = 0
		p0.ken = 1
		p0.idir = 0
		p0.jdir = 1
		p0.kdir = 2
		p0.kind = 5	
		p1 = ts_tstream_type.TstreamPatch()
		p1.bid = bid1
		p1.nxbid = bid0
		p1.ist = 0
		p1.ien = N2+N3-1
		p1.jst = 0
		p1.jen = b1.nj
		p1.kst = b1.nk-1
		p1.ken = b1.nk
		p1.idir = 0
		p1.jdir = 1
		p1.kdir = 2
		p1.kind = 5

		p0.pid = g.add_patch(bid0, p0)
		p1.pid = g.add_patch(bid1, p1)
		p1.nxpid = p0.pid
		p0.nxpid = p1.pid

		### PERIODIC AFTER TE
		bid0 = bidsmid2[-1]
		bid1 = bidsmid2[-1]
		b0 = g.get_block(bid0)
		b1 = g.get_block(bid1)

		p0 = ts_tstream_type.TstreamPatch()
		p0.bid = bid0
		p0.nxbid = bid1
		p0.ist = 0
		p0.ien = b0.ni
		p0.jst = 0
		p0.jen = b0.nj
		p0.kst = 0
		p0.ken = 1
		p0.idir = 0
		p0.jdir = 1
		p0.kdir = 2
		p0.kind = 5	

		p1 = ts_tstream_type.TstreamPatch()
		p1.bid = bid1
		p1.nxbid = bid0
		p1.ist = 0
		p1.ien = b1.ni
		p1.jst = 0
		p1.jen = b1.nj
		p1.kst = b1.nk-1
		p1.ken = b1.nk
		p1.idir = 0
		p1.jdir = 1
		p1.kdir = 2
		p1.kind = 5

		p0.pid = g.add_patch(bid0, p0)
		p1.pid = g.add_patch(bid1, p1)
		p1.nxpid = p0.pid
		p0.nxpid = p1.pid

		### DOWNSTREAM PERIODIC
		bid = bidsdown[-1]
		bid0 = bid
		bid1 = bid
		b0 = g.get_block(bid0)
		b1 = g.get_block(bid1)
		p0 = ts_tstream_type.TstreamPatch()
		p0.bid = bid0
		p0.nxbid = bid1
		p0.ist = 0
		p0.ien = b0.ni
		p0.jst = 0
		p0.jen = b0.nj
		p0.kst = 0
		p0.ken = 1
		p0.idir = 0
		p0.jdir = 1
		p0.kdir = 2
		p0.kind = 5	

		p1 = ts_tstream_type.TstreamPatch()
		p1.bid = bid1
		p1.nxbid = bid0
		p1.ist = 0
		p1.ien = b1.ni
		p1.jst = 0
		p1.jen = b1.nj
		p1.kst = b1.nk-1
		p1.ken = b1.nk
		p1.idir = 0
		p1.jdir = 1
		p1.kdir = 2
		p1.kind = 5

		p0.pid = g.add_patch(bid0, p0)
		p1.pid = g.add_patch(bid1, p1)
		p1.nxpid = p0.pid
		p0.nxpid = p1.pid

		### BLOCK LINKING PATCHES
		### UPSTREAM TO PASSAGE
		bid0 = bidsup[-1]
		bid1 = bidsmid1[-1]
		b0 = g.get_block(bid0)
		b1 = g.get_block(bid1)

		p0 = ts_tstream_type.TstreamPatch()
		p0.bid = bid0
		p0.nxbid = bid1
		p0.ist = b0.ni-1
		p0.ien = b0.ni
		p0.jst = 0
		p0.jen = b0.nj
		p0.kst = 0
		p0.ken = b0.nk
		p0.idir = 0
		p0.jdir = 1
		p0.kdir = 2
		p0.kind = 5	

		p1 = ts_tstream_type.TstreamPatch()
		p1.bid = bid1
		p1.nxbid = bid0
		p1.ist = 0
		p1.ien = 1
		p1.jst = 0
		p1.jen = b1.nj
		p1.kst = 0
		p1.ken = b1.nk
		p1.idir = 0
		p1.jdir = 1
		p1.kdir = 2
		p1.kind = 5

		p0.pid = g.add_patch(bid0, p0)
		p1.pid = g.add_patch(bid1, p1)
		p1.nxpid = p0.pid
		p0.nxpid = p1.pid

		### PASSAGE TO TE BLOCK
		bid0 = bidsmid1[-1]
		bid1 = bidsmid2[-1]
		b0 = g.get_block(bid0)
		b1 = g.get_block(bid1)

		p0 = ts_tstream_type.TstreamPatch()
		p0.bid = bid0
		p0.nxbid = bid1
		p0.ist = b0.ni-1
		p0.ien = b0.ni
		p0.jst = 0
		p0.jen = b0.nj
		p0.kst = 0
		p0.ken = b0.nk
		p0.idir = 0
		p0.jdir = 1
		p0.kdir = 2
		p0.kind = 5	

		p1 = ts_tstream_type.TstreamPatch()
		p1.bid = bid1
		p1.nxbid = bid0
		p1.ist = 0
		p1.ien = 1
		p1.jst = 0
		p1.jen = b1.nj
		### LOGIC BASED ON BLADE EXIT ANGLE
		if row.X2[N] >= 0:               ### i.e. for stators
			p1.kst = 0
			p1.ken = b0.nk
		else:		         		### for rotors
			p1.kst = b1.nk-b0.nk
			p1.ken = b1.nk
		p1.idir = 0
		p1.jdir = 1
		p1.kdir = 2
		p1.kind = 5
	
		p0.pid = g.add_patch(bid0, p0)
		p1.pid = g.add_patch(bid1, p1)
		p1.nxpid = p0.pid
		p0.nxpid = p1.pid

		### TE TO DOWNSTREAM
		bid0 = bidsmid2[-1]
		bid1 = bidsdown[-1]
		b0 = g.get_block(bid0)
		b1 = g.get_block(bid1)
		p0 = ts_tstream_type.TstreamPatch()
		p0.bid = bid0
		p0.nxbid = bid1
		p0.ist = b0.ni-1
		p0.ien = b0.ni
		p0.jst = 0
		p0.jen = b0.nj
		p0.kst = 0
		p0.ken = b0.nk
		p0.idir = 0
		p0.jdir = 1
		p0.kdir = 2
		p0.kind = 5	

		p1 = ts_tstream_type.TstreamPatch()
		p1.bid = bid1
		p1.nxbid = bid0
		p1.ist = 0
		p1.ien = 1
		p1.jst = 0
		p1.jen = b1.nj
		p1.kst = 0
		p1.ken = b1.nk
		p1.idir = 0
		p1.jdir = 1
		p1.kdir = 2
		p1.kind = 5

		p0.pid = g.add_patch(bid0, p0)
		p1.pid = g.add_patch(bid1, p1)
		p1.nxpid = p0.pid
		p0.nxpid = p1.pid

		### inlet or mixing plane
		############################################
		### INLET CONDITIONS USER INPUT REQUIRED ###
		############################################

		if rownumber == 0:
			### inlet patch
			bid0=0#bidsup[0]
			b0 = g.get_block(bid0)
			p0 = ts_tstream_type.TstreamPatch()
			p0.bid = bid0
			p0.ist = 0
			p0.ien = 1
			p0.jst = 0
			p0.jen = b0.nj
			p0.kst = 0
			p0.ken = b0.nk
			p0.idir = 0
			p0.jdir = 1
			p0.kdir = 2
			p0.kind = 0
			p0.pid = g.add_patch(bid0, p0)

			pstag = np.zeros((p0.ken-p0.kst, p0.jen-p0.jst, p0.ien-p0.ist), np.float32)
			tstag = np.zeros((p0.ken-p0.kst, p0.jen-p0.jst, p0.ien-p0.ist), np.float32)
			pitch = np.zeros((p0.ken-p0.kst, p0.jen-p0.jst, p0.ien-p0.ist), np.float32)
			yaw = np.zeros((p0.ken-p0.kst, p0.jen-p0.jst, p0.ien-p0.ist), np.float32)
			radius_in=(g.get_bp('r',0)[0,:,0])
			radius_in2=(g.get_bp('r',0)[:,:,0])
			radiust_in2=(g.get_bp('rt',0)[:,:,0])

			### generate inlet conditions -  default set up for inlet angle equal to X1  ######0 relative yaw angle######
			vxi = np.zeros((len(radius_in)))
			vti = np.zeros((len(radius_in)))
			vi = np.zeros((len(radius_in)))
			yawi = np.zeros((len(radius_in)))
			psi = np.zeros((len(radius_in)))
			tsi = np.zeros((len(radius_in)))
			t0i = np.zeros((len(radius_in)))
			p0i = np.zeros((len(radius_in)))
			p0_reli = np.zeros((len(radius_in)))
			omega = row.RPM*2.*np.pi/60.
			rhub_in=radius_in[0]
			span = radius_in[-1]-radius_in[0]
			walldist=np.zeros((len(radius_in)))
			for iw in range(len(walldist)):
				walldist[iw] = (min((radius_in[iw]-radius_in[0])/span,-(radius_in[iw]-radius_in[-1])/span))**0.5
				vxi[iw]=(vI[0]*0.9 * min(walldist[iw]/0.1**0.5,1)+vI[0]*0.1)*np.cos(np.radians(row.X1[0]))

			vti[:]=np.tan(np.radians(row.X1[0]))*vxi[:]
			vi[:]=(vxi[:]**2+vti[:]**2)**0.5
			yawi[:]=np.degrees(np.arctan(vti[:]/vxi[:]))
			psi[:]=PSHUBin#+ROin*0.5*(omega**2)*(radius_in[:]**2-rhub_in**2)
			p0i[:]=P0IN-0.5*ROin*(vi.max()**2-vi[:]**2)#+0.5*ROin*(omega**2)*radius_in[:]**2
			tsi[:]=psi[:]/(ROin*Rgas)
			t0i[:]=T0IN#-(vi.max()**2-VXin**2)/(2.*Cp)

			plt.figure(num = 'INLET CONDITIONS')
			plt.subplot(221)
			plt.plot(vxi,radius_in,'-xb')
			plt.xlabel('VX')
			plt.subplot(222)
			plt.plot(vti,radius_in,'-xb')
			plt.xlabel('VT')
			plt.subplot(223)
			plt.plot(p0i,radius_in,'-xb')
			plt.xlabel('P0')
			plt.subplot(224)
			plt.plot(psi,radius_in,'-xb')
			plt.xlabel('Ps')
			#plt.show()

			
			for j in range(len(pstag[:,0])):
			    for k in range(len(pstag[0,:])):
				pstag[j,k]=p0i[k]
				tstag[j,k]=t0i[k]
				yaw[j,k]=yawi[k]
				pitch[j,k]=0.0


			g.set_pv("rfin",ts_tstream_type.float,p0.bid,p0.pid,0.3)
			g.set_pv("sfinlet",ts_tstream_type.float,p0.bid,p0.pid,0.0)
		     	g.set_pp("pstag", ts_tstream_type.float, p0.bid, p0.pid, pstag)
		    	g.set_pp("tstag", ts_tstream_type.float, p0.bid, p0.pid, tstag)    	
		    	g.set_pp("pitch", ts_tstream_type.float, p0.bid, p0.pid, pitch)   
		    	g.set_pp("yaw", ts_tstream_type.float, p0.bid, p0.pid, yaw)
			g.set_bv("nimixl",ts_tstream_type.int,p0.bid,5)			### set mixing length in first block
			g.set_bv("sfin_mul",ts_tstream_type.float,p0.bid,1.)			### set mixing length in first block

		else:

			print 'ADDING MIXING PLANE'
			### mixing plane to upstream block
			bid0 = bidsdown[-2]
			bid1 = bidsup[-1]
			b0 = g.get_block(bid0)
			b1 = g.get_block(bid1)

			p0 = ts_tstream_type.TstreamPatch()
			p0.bid = bid0
			p0.nxbid = bid1
			p0.ist = b0.ni-1
			p0.ien = b0.ni
			p0.jst = 0
			p0.jen = b0.nj
			p0.kst = 0
			p0.ken = b0.nk
			p0.idir = 0
			p0.jdir = 1
			p0.kdir = 2
			p0.kind = 2	

			p1 = ts_tstream_type.TstreamPatch()
			p1.bid = bid1
			p1.nxbid = bid0
			p1.ist = 0
			p1.ien = 1
			p1.jst = 0
			p1.jen = b1.nj
			p1.kst = 0
			p1.ken = b1.nk
			p1.idir = 0
			p1.jdir = 1
			p1.kdir = 2
			p1.kind = 2

			p0.pid = g.add_patch(bid0, p0)
			p1.pid = g.add_patch(bid1, p1)
			p1.nxpid = p0.pid
			p0.nxpid = p1.pid

			#g.set_bv("nimixl",ts_tstream_type.int,p1.bid,5)			### set mixing length in first block
			
		#### exit block if final row ###

		#### exit patch
		mdot =MDOT
		if row.Nblades>3:
			Aestimate = np.pi*(rtip**2-rhub**2)
		else:
			mdot = mdot/nbladesfake
			Aestimate = (rtip-rhub)*S

		if rownumber == nrows-1:
			bid0=bidsdown[-1]
			b0 = g.get_block(bid0)
			p02 = ts_tstream_type.TstreamPatch()
			p02.bid = bid0
			p02.ist = b0.ni-1
			p02.ien = b0.ni
			p02.jst = 0
			p02.jen = b0.nj
			p02.kst = 0
			p02.ken = b0.nk
			p02.idir = 0
			p02.jdir = 1
			p02.kdir = 2
			p02.kind = 1
			p02.pid = g.add_patch(bid0, p02)
			g.set_pv("pout",ts_tstream_type.float,bid0,p02.pid,pE[-1])
			g.set_pv("ipout",ts_tstream_type.int,bid0,p02.pid,+3)
			g.set_pv("throttle_type",ts_tstream_type.int,bid0,p02.pid,1)
			g.set_pv("throttle_target",ts_tstream_type.float,bid0,p02.pid,mdot)
			g.set_pv("throttle_k0",ts_tstream_type.float,bid0,p02.pid,20.)
			g.set_pv("throttle_k1",ts_tstream_type.float,bid0,p02.pid,100.)
			g.set_pv("throttle_k2",ts_tstream_type.float,bid0,p02.pid,1600.0)


	
		rownumber=rownumber+1
		xshift=xend



	for bid in g.get_block_ids():
		print ''
		print 'Block Number: ', bid
		print 'RPM', g.get_bv('rpm',bid)
		print 'RPMi1', g.get_bv('rpmi1',bid)
		print 'RPMi2', g.get_bv('rpmi2',bid)
		print 'RPMj1', g.get_bv('rpmj1',bid)
		print 'RPMj2', g.get_bv('rpmj2',bid)
		print 'RPMk1', g.get_bv('rpmk1',bid)
		print 'RPMk2', g.get_bv('rpmk2',bid)
	
	set_default(g,Cp,gam,Rgas,visc)

	g.set_bv('nimixl', ts_tstream_type.int, 0, 5)
	#g.set_bv('fmgrid', ts_tstream_type.float, 0, 0.0)
	#g.set_bv('fmgrid', ts_tstream_type.float, 3, 0.0)
	#g.set_bv('poisson_fmgrid', ts_tstream_type.float, 3, 0.0)
	
	g.set_bv('dampin_mul', ts_tstream_type.float, 2,0.2)
	g.set_bv('dampin_mul', ts_tstream_type.float, 6,0.2)
	g.set_bv('dampin_mul', ts_tstream_type.float, 3,0.5)
	g.set_bv('dampin_mul', ts_tstream_type.float, 7,0.5)
	#g.set_bv('fmgrid', ts_tstream_type.float, 3,0.0)
	#g.set_bv('fmgrid', ts_tstream_type.float, 4,0.0)
	#g.set_bv('poisson_fmgrid', ts_tstream_type.float, 13, 0.0)

	for bid in g.get_block_ids():
		print ''
		print 'Block Number: ', bid
		print 'nimixl', g.get_bv('nimixl',bid)
		print 'fmgrid', g.get_bv('fmgrid',bid)
		print 'poisson_fmgrid', g.get_bv('poisson_fmgrid',bid)

	'''for av in g.get_av_ids():
		#print av, g.get_av(av)'''
    	ts_tstream_load_balance.load_balance(g, 1)
	fname='MESH'

    	g.write_hdf5(fname+'.hdf5')
    	g.write_xdmf(fname+".xdmf",'x','r','rt')

	print 'write successful'



#########################
####Start Main Routine###
#########################

#####################################
### DEFINE OPERATING CONDITIONS #####
#####################################

### INLET
V_correct = 1.0
Cp = 5187.
gam=1.6625
Rgas=Cp*(1.-1./gam)
visc = (4.18e-5)*V_correct

Mdot = 16.0*V_correct
P0in = 14500000.0
T0in = 950.0
RPM_shaft = 6782.*V_correct

ROin = P0in/(Rgas*T0in)

Rin = 0.404
Hin = 0.0084#226
Ain = np.pi*((Rin+Hin/2.)**2-(Rin-Hin/2.)**2)
for it in range(5):
	VXin = Mdot/(ROin*Ain)
	PSin=P0in-0.5*ROin*VXin**2


#################################
### GENERATE MACHINE GEOMETRY ###
#################################

ROWS = []	### create an array of rows
for irow in range(1):
	S1 = ROW()
	S1.name = 'S1'
	S1.X1=[-6.74,-6.68,-6.62]				### INLET METAL ANGLES: [Hub,Casing]
	S1.X2=[73,71,69]				### EXIT METAL ANGLES: [Hub,Casing]
	S1.Cx=0.0050				### AXIAL CHORDS
	S1.RADIUS =[0.404-0.0042,0.404,0.404+0.0042]#[0.214,0.2141]			### RADIUS OF SECTIONS
	S1.RPM = 0.0#-6782.0/2.			### RPM OF ROW
	S1.Nblades = 1.0			### # of blades in row/if below 3 this forms pitch to chord
	S1.shroud_hub = True			### if there is a hub shroud
	S1.shroud_cas = False			### if there is a casing shroud
	S1.upfrac=[0.5,0.3,0.1]		### fraction of chord upstream for [mx,shroudstart,shroudend]
	S1.downfrac=[0.5,0.1,0.3]		### fraction of chord upstream for [mx,shroudend,shroudend]
	S1.shroudthickness = 0.001	### thickness of shroud [upstream,downstream]
	S1.shroudclearance = 0.00025	### shroud clearance at upstream and downstream edge
	S1.Nfins = 0			### add fins
	S1.finclearance = 0.000125
	S1.finthickness = 0.05
	S1.thkTE = 0.3/(1000.*S1.Cx)

	ROWS.append(S1)
	#MESH(ROWS)
	R1 = ROW()
	R1.name = 'R1'
	R1.X1=[20.70,9.68,2.70]				### INLET METAL ANGLES: [Hub,Casing]
	R1.X2=[-71,-73,-71]				### EXIT METAL ANGLES: [Hub,Casing]
	R1.Cx=0.0052				### AXIAL CHORDS
	R1.RADIUS =[0.404-0.0044+0.0001,0.404,0.404+0.0044-0.0001]#[0.214,0.2141]	### RADIUS OF SECTIONS
	R1.RPM = RPM_shaft*V_correct			### RPM OF ROW
	R1.Nblades = 0.9			### # of blades in row/if below 3 this forms pitch to chord
	R1.shroud_hub = False			### if there is a hub shroud
	R1.shroud_cas = True	 		### if there is a casing shroud
	R1.upfrac=[0.5,0.3,0.1]		### fraction of chord upstream for [mx,shroudstart,shroudend]
	R1.downfrac=[0.5,0.1,0.3]		### fraction of chord upstream for [mx,shroudstart,shroudend]
	R1.shroudthickness = 0.001	### thickness of shroud [upstream,downstream]
	R1.shroudclearance = 0.00025	### shroud clearance at upstream and downstream edge
	R1.Nfins = 0			### add fins
	R1.finclearance = 0.000125
	R1.finthickness = 0.05
	R1.thkTE = 0.3/(1000.*R1.Cx)

	ROWS.append(R1)

S2 = ROW()
S2.name = 'S2'
S2.X1=[-5.16,-5.11,-5.06]				### INLET METAL ANGLES: [Hub,Casing]
S2.X2=[76.38+2.0,76.26+1.5,76.14+2.0]				### EXIT METAL ANGLES: [Hub,Casing]
S2.Cx=0.0053				### AXIAL CHORDS
S2.RADIUS =[0.484-0.0045+0.0001,0.484,0.484+0.0045-0.0001]#[0.214,0.2141]			### RADIUS OF SECTIONS
S2.RPM = 0.0#-6782.0/2.			### RPM OF ROW
S2.Nblades = 1.0			### # of blades in row/if below 3 this forms pitch to chord
S2.shroud_hub = True			### if there is a hub shroud
S2.shroud_cas = False			### if there is a casing shroud
S2.upfrac=[0.5,0.3,0.1]		### fraction of chord upstream for [mx,shroudstart,shroudend]
S2.downfrac=[0.5,0.1,0.3]		### fraction of chord upstream for [mx,shroudend,shroudend]
S2.shroudthickness = 0.001	### thickness of shroud [upstream,downstream]
S2.shroudclearance = 0.00025	### shroud clearance at upstream and downstream edge
S2.Nfins = 0			### add fins
S2.finclearance = 0.000125
S2.finthickness = 0.05
S2.thkTE = 0.3/(1000.*S2.Cx)

#ROWS.append(S2)

R2 = ROW()
R2.name = 'R2'
R2.X1=[9.43,5.11,0.81]				### INLET METAL ANGLES: [Hub,Casing]
R2.X2=[-76.14-1.0,-76.26,-76.38-1.0]				### EXIT METAL ANGLES: [Hub,Casing]
R2.Cx=0.0055				### AXIAL CHORDS
R2.RADIUS =[0.484-0.0047+0.00015,0.484,0.484+0.0047-0.00015]#[0.214,0.2141]	### RADIUS OF SECTIONS
R2.RPM = 6782.*V_correct			### RPM OF ROW
R2.Nblades = 0.9			### # of blades in row/if below 3 this forms pitch to chord
R2.shroud_hub = False			### if there is a hub shroud
R2.shroud_cas = True	 		### if there is a casing shroud
R2.upfrac=[0.5,0.3,0.1]		### fraction of chord upstream for [mx,shroudstart,shroudend]
R2.downfrac=[0.5,0.1,0.3]		### fraction of chord upstream for [mx,shroudstart,shroudend]
R2.shroudthickness = 0.001	### thickness of shroud [upstream,downstream]
R2.shroudclearance = 0.00025	### shroud clearance at upstream and downstream edge
R2.Nfins = 0			### add fins
R2.finclearance = 0.000125
R2.finthickness = 0.05
R2.thkTE = 0.3/(1000.*R2.Cx)

#ROWS.append(R2)

MESH(ROWS,VXin,ROin,PSin,gam,Rgas,Cp,visc)


plt.show()


















