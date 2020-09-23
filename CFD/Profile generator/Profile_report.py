import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import os,sys,time,copy,math,scipy,numpy
from scipy import interpolate as sciiint
import pylab
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
        ##printpoint
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
        plt.ylabel('fraction')
        plt.show()
    return(np.array(x),np.array(y),CPs)
    
def GEN_TURNING(Xi1,Xi2,controls,plotting):
    """ Generate x,y of a bezier from control points """
    
    spacings = np.linspace(0,1,len(controls)+2)[1:-1]
    #printspacings

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
        ##printpoint
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
    ##printS
    
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
        
        #print'Peak thickness:',2.*TF.max()*100,'% Cx'
        for i in range(len(TF)):
            if TF[i]==TF.max():
                #print'@ ',100.*psi[i],' % camber line'
                break
            
        Bte = np.degrees(np.arctan((TF[-2]-TF[-1])/(psi[-1]-psi[-2])))
        #print'Boat tail:',Bte
    
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
        ##printA,portiondroop
	#plotting=1							### calculate normal 2 as 1/grad1

        return(norm)    
        
def F_Make(Xi1,Xi2,controls_cam,controls_thk_x,controls_thk_y,Tte,Oval_frac,plotting):        ### function that compiles a blade from reduced variables 


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
    ##printpsi
    #thk = F_TF(psi,controls_thk,Tte,plotting=1,spacing=None)         ### calculate thickness distribution
    
    thk = F_TF_bez(controls_thk_x,controls_thk_y,Tte,Oval_frac,plotting=1,spacing=None)


    Z_U = y_cam+thk*norm[:,1]	### apply thickness onto camber line for upper
    Z_L = y_cam-thk*norm[:,1]	### apply thickness onto camber line for lower
    X_U = x_cam+thk*norm[:,0]	### repeat upper for x rather than y
    X_L = x_cam-thk*norm[:,0]

    if plotting==1:		### optional plot out of complete blade
	##print'Gamma:',Gamma,' Xi1:',Xi1,' Xi2:', Xi2
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
    ##printtrans(0.00001), b, math.
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
    ##printS
    
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
        
        #print'Peak thickness:',2.*TF.max()*100,'% Cx'
        for i in range(len(TF)):
            if TF[i]==TF.max():
                #print'@ ',100.*psi[i],' % camber line'
                break
            
        Bte = np.degrees(np.arctan((TF[-2]-TF[-1])/(psi[-1]-psi[-2])))
        #print'Boat tail:',Bte
    
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
        ##printA,portiondroop
	#plotting=1							### calculate normal 2 as 1/grad1
        return(norm)    
        
def F_Make(Xi1,Xi2,controls_cam,controls_thk_x,controls_thk_y,Tte,Oval_frac,plotting):        ### function that compiles a blade from reduced variables 


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
    ##printpsi
    #thk = F_TF(psi,controls_thk,Tte,plotting=1,spacing=None)         ### calculate thickness distribution
    
    thk = F_TF_bez(controls_thk_x,controls_thk_y,Tte,Oval_frac,plotting=0,spacing=None)


    Z_U = y_cam+thk*norm[:,1]	### apply thickness onto camber line for upper
    Z_L = y_cam-thk*norm[:,1]	### apply thickness onto camber line for lower
    X_U = x_cam+thk*norm[:,0]	### repeat upper for x rather than y
    X_L = x_cam-thk*norm[:,0]

    if plotting==1:		### optional plot out of complete blade
	##print'Gamma:',Gamma,' Xi1:',Xi1,' Xi2:', Xi2
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
    return(X_L,Z_L,X_U,Z_U,thk)


def Profile2(X1,X2,TKTE):		### function based on profgen.f (JDD) - GENERATES SIMPLE BLADE PROFILES


    #### THIS FUNCTION SHOULD BE REPLACED AT SOME POINT TO USE A BETTER BLADE SECTION GENERATOR
    ### EXCUSE THE WEIRD FORMULATION, THIS IS CONVERSION FROM FORTRAN

    #X2 = X2#*1.07

    ### USER HARD CODED VALUES BELOW - COULD BE TAKEN OUTSIDE FUNCTION CJC

    controls_cam = [0.0,-0.3,-0.8]
    Rle = 0.05
    Tte = TKTE
    Beta_te = 0.1
    Oval_frac = 0.7
    #controls_cam = [0,0.2,0.0]
    controls_thk_x = [0.0,0.5,1.0]
    controls_thk_y = [(2.*Rle)**0.5*(1.-Oval_frac),0.25,np.tan(np.radians(Beta_te))+Tte/2.]

    (XSIN,YSIN,XPIN,YPIN,thk) = F_Make(X1,X2,controls_cam,controls_thk_x,controls_thk_y,Tte,Oval_frac,plotting=False)
    
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
    points = 2000
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

    #thkcx = np.arange
    TE_circ_r = 0.5*((abs(XSIN[-1]-XPIN[-1]))**2+(abs(YSIN[-1]-YPIN[-1]))**2)**0.5
    TE_circ_cx = 0.5*(XSIN[-1]+XPIN[-1])
    TE_circ_cy = 0.5*(YSIN[-1]+YPIN[-1])
    TEcx = np.arange(XSIN[-1],TE_circ_cx+TE_circ_r,0.0000001)
    TEcx2 = np.arange(XPIN[-1],TE_circ_cx+TE_circ_r,0.0000001)
    TEcy = np.asarray([TE_circ_cy-(TE_circ_r**2-(i-TE_circ_cx)**2)**0.5 for i in TEcx])
    TEcy2 = np.asarray([TE_circ_cy+(TE_circ_r**2-(i-TE_circ_cx)**2)**0.5 for i in TEcx2])
    #print TE_circ_cy, TE_circ_cx,TE_circ_r, XSIN[-1],XPIN[-1]
    #print TE_circ_y, TE_circ_x
    #print YPIN.shape,XPIN.shape,YIN.shape,XIN.shape,YSIN.shape,XSIN.shape
    maxes = np.array([np.amax(XIN),np.amax(XPIN),np.amax(XSIN),np.amax(TEcx),np.amax(TEcx2)])
    mins = np.array([np.amin(XIN),np.amin(XPIN),np.amin(XSIN)])
    cx = np.amax(maxes)-np.amin(mins)
    ptc = 0.7
    sx = 1.7
    sy = 0.8
    plt.show()
    if 1==1:            ### Plot with lines
        i=25
        d=25
        arrows=0.015*cx
        n = 0.00001
        
        dist_min = 100
        throat = 0
        for j in range(len(YPIN)):
            if np.sqrt((XSIN[-1]-XPIN[j])**2+(YSIN[-1]+ptc*cx-YPIN[j])**2)<dist_min:
                dist_min = np.sqrt((XSIN[-1]-XPIN[j])**2+(YSIN[-1]+ptc*cx-YPIN[j])**2)
                throat = j
    
        cam1m = (YIN[1]-YIN[0])/(XIN[1]-XIN[0])
        lexs = np.arange(XSIN[points*0.04]-0.05*cx,XSIN[points*0.04],n)
        perp_x = np.arange(XIN[0],lexs[-1])
        perp_y = [YSIN[points*0.04]+(-1/cam1m)*(l-XSIN[points*0.04]) for l in perp_x]
        min_dis = 100
        po = 0
        for j in range(len(YPIN[:np.argmax(YPIN)])):
            for k in range(len(perp_x)):
                if np.sqrt((perp_x[k]-XPIN[j])**2+(perp_y[k]-YPIN[j])**2)<min_dis:
                    po = j
                    min_dis = np.sqrt((perp_x[k]-XPIN[j])**2+(perp_y[k]-YPIN[j])**2)
        lexp = np.arange(XPIN[po]-0.05*cx,XPIN[po],n)
        throat_bx = (XSIN[-1]+XPIN[throat])/2
        throat_by = (YSIN[-1]+cx*ptc+YPIN[throat])/2
        throat_dx = XSIN[-1]-XPIN[throat]
        throat_dy = YSIN[-1]+ptc*cx-YPIN[throat]
        throat_x = abs(throat_dx/2)-arrows*np.cos(np.arctan(throat_dy/throat_dx))
        throat_y = abs(throat_dy/2)-arrows*np.sin(np.arctan(throat_dy/throat_dx))
        guide1x = np.arange(XPIN[throat],XSIN[-1],0.001)
        guide1y = [YSIN[-1]+ptc*cx+(throat_dy/throat_dx)*(q-XSIN[-1]) for q in guide1x]
        ax1 = np.arange(XIN[0]-0.2*cx,XIN[0],0.001)
        cam1x = np.arange(XIN[0]-0.2*cx*np.cos(np.arctan(cam1m)),XIN[0],n)
        v1x = XIN[0]-0.2*cx*np.cos(np.arctan(cam1m)+np.radians(i))
        v1y = YIN[0]+ptc*cx-0.2*cx*np.sin(np.arctan(cam1m)+np.radians(i))
        v1dx = abs(v1x-XIN[0])-arrows*np.cos(np.arctan(cam1m)+np.radians(i))
        v1dy = abs(v1y-YIN[0]-ptc*cx)-arrows*np.sin(np.arctan(cam1m)+np.radians(i))
        guide2x = np.arange(XIN[0]-0.2*cx*np.cos(np.arctan(cam1m)+np.radians(i)),XIN[0],0.001)
        guide2y = [YIN[0]+ptc*cx+(v1dy/v1dx)*(h-XIN[0]) for h in guide2x]
        cam2m = (YIN[-1]-YIN[-2])/(XIN[-1]-XIN[-2])
        c0 = YIN[-1]-cam2m*XIN[-1]
        a = 1+cam2m**2
        b = 2*cam2m*c0-2*TE_circ_cy*cam2m-2*TE_circ_cx
        c = TE_circ_cx**2+TE_circ_cy**2+c0**2-2*TE_circ_cy*c0-TE_circ_r**2
        endx = (-b+np.sqrt(b**2-4*a*c))/(2*a)
        endy = cam2m*endx+c0
        cam2x = np.arange(endx,endx+0.2*cx*np.cos(np.arctan(cam2m)),n)
        v2x = endx-0.2*cx*np.cos(np.arctan(cam2m)+np.radians(d))
        v2y = endy+ptc*cx-0.2*cx*np.sin(np.arctan(cam2m)+np.radians(d))
        v2dx = abs(v2x-endx)-arrows*np.cos(np.arctan(cam2m)+np.radians(d))
        v2dy = abs(v2y-endy-ptc*cx)+arrows*np.sin(np.arctan(cam2m)+np.radians(d))
        guide3x = np.arange(endx,endx+0.2*cx*np.cos(np.arctan(cam2m)+np.radians(d)),0.001)
        guide3y = [endy+ptc*cx-(v2dy/v2dx)*(h-endx) for h in guide3x]
        ax2 = np.arange(endx,endx+(0.2+0.02)*cx,n)
        stagm = (YIN[0]-endy)/(XIN[0]-endx)
        stagx = np.arange(XIN[0]-0.2*cx*np.cos(np.arctan(stagm)),XIN[0],n)
        chord = np.sqrt((XIN[0]-endx)**2+(YIN[0]-endy)**2)
        ch_arrx = arrows*np.cos(np.arctan(abs(stagm)))
        ch_arry = arrows*np.sin(np.arctan(abs(stagm)))
        thkx = XIN[np.argmax(thk)]
        thky = YIN[np.argmax(thk)]
        thk_circ_r = 0.92*thk[np.argmax(thk)]
        thkcx = np.arange(thkx-0.99999999*thk_circ_r,thkx+0.99999999*thk_circ_r,0.000001)
        thkcy = np.asarray([thky-(thk_circ_r**2-(i-thkx)**2)**0.5 for i in thkcx])
        thkcy2 = np.asarray([thky+(thk_circ_r**2-(i-thkx)**2)**0.5 for i in thkcx])
        thkm =-(XIN[np.argmax(thk)+1]-XIN[np.argmax(thk)-1])/(YIN[np.argmax(thk)+1]-YIN[np.argmax(thk)-1])
        thk_delx = 0.8*thk_circ_r*np.cos(np.arctan(thkm))
        thk_dely = 0.8*thk_circ_r*np.sin(np.arctan(thkm))
        texs = np.arange(XSIN[-1],XSIN[-1]+0.05*cx,n)
        texp = np.arange(XPIN[-1],XPIN[-1]+0.05*cx,n)
        cam_max = np.argmax(YIN)
        cam_maxy = np.arange(YIN[0],YIN[cam_max],n)
        cx1 = np.arange(YPIN[np.argmin(XPIN)]+ptc*cx,np.amax(YPIN)+(ptc+0.06)*cx,n)
        cx2 = np.arange(TEcy[np.argmax(TEcx)]+ptc*cx,np.amax(YPIN)+(ptc+0.06)*cx,n)
        te_arr_bx = (texs[-1]+texp[-1])/2
        te_arr_by = ((YSIN[-1]+cam2m*(texs[-1]-XSIN[-1]))+(YPIN[-1]+cam2m*(texp[-1]-XPIN[-1])))/2
        te_arr_dx = texs[-1]-texp[-1]
        te_arr_dy = (YSIN[-1]+cam2m*(texs[-1]-XSIN[-1]))-(YPIN[-1]+cam2m*(texp[-1]-XPIN[-1]))
        te_arr_x = abs((texs[-1]-texp[-1])/2)-0.5*TE_circ_r*np.cos(np.arctan(te_arr_dy/te_arr_dx))
        te_arr_y = abs(((YSIN[-1]+cam2m*(texs[-1]-XSIN[-1]))-(YPIN[-1]+cam2m*(texp[-1]-XPIN[-1])))/2)-0.5*TE_circ_r*np.sin(np.arctan(te_arr_dy/te_arr_dx))
        le_arr_bx = (lexs[0]+lexp[0])/2
        le_arr_by = ((YSIN[points*0.04]+cam1m*(lexs[0]-XSIN[points*0.04]))+(YPIN[po]+cam1m*(lexp[0]-XPIN[po])))/2
        le_arr_dx = abs(lexs[0]-lexp[0])
        le_arr_dy = abs((YSIN[points*0.04]+cam1m*(lexs[0]-XSIN[points*0.04]))-(YPIN[po]+cam1m*(lexp[0]-XPIN[po])))
        le_arr_x = le_arr_dx/2-0.5*TE_circ_r*np.cos(np.arctan(le_arr_dy/le_arr_dx))
        le_arr_y = le_arr_dy/2-0.5*TE_circ_r*np.sin(np.arctan(le_arr_dy/le_arr_dx))
        throat_m = (YPIN[throat]-YPIN[throat+1])/(XPIN[throat]-XPIN[throat+1])
        uncov_max_x = XPIN[-1]+(YPIN[throat]+throat_m*(XPIN[-1]-XPIN[throat])-YPIN[-1])*np.cos(abs(np.arctan(throat_m)))*np.sin(abs(np.arctan(throat_m)))
        uncov_x = np.arange(XPIN[throat],uncov_max_x,n)
        
        plt.axis('equal')
        plt.margins(0.05,0.05)
        lw = 1
        plt.plot(ax1,[YIN[0]+ptc*cx for i in ax1],'black',linestyle='--',linewidth=lw)
        plt.plot(ax2,[endy+ptc*cx for i in ax2],'black',linestyle='--',linewidth=lw)
#       plt.arrow(endx+(0.2+0.02)*cx,endy+ptc*cx/2,0,(ptc/2-arrows/cx)*cx,ec=None,fc='black',head_length=arrows,head_width=arrows)
#       plt.arrow(endx+(0.2+0.02)*cx,endy+ptc*cx/2,0,-(ptc/2-arrows/cx)*cx,ec=None,fc='black',head_length=arrows,head_width=arrows)
        plt.plot(ax2,[endy for i in ax2],'black',linestyle='--',linewidth=lw)
        plt.plot(cam1x,[YIN[0]+ptc*cx+cam1m*(i-XIN[0]) for i in cam1x],'black',linestyle='--',dashes=(5,4),linewidth=lw)
        plt.plot(cam2x,[YIN[-1]+ptc*cx+cam2m*(i-XIN[-1]) for i in cam2x],'black',linestyle='--',dashes=(5,4),linewidth=lw)
        plt.plot(stagx,[YIN[0]+ptc*cx+stagm*(i-XIN[0]) for i in stagx],'black',linestyle='--',dashes=(5,4),linewidth=lw)
        #plt.arrow((endx+XIN[0])/2,(endy+YIN[0]+2*ptc*cx)/2,-(endx-XIN[0])/2+ch_arrx,-(endy-YIN[0])/2-ch_arry,ec=None,fc='black',head_length=arrows,head_width=arrows)
        #plt.arrow((endx+XIN[0])/2,(endy+YIN[0]+2*ptc*cx)/2,(endx-XIN[0])/2-ch_arrx,(endy-YIN[0])/2+ch_arry,ec=None,fc='black',head_length=arrows,head_width=arrows)
#       plt.arrow(thkx,thky,thk_delx,thk_dely,ec=None,fc='black',head_length=0.2*thk_circ_r,head_width=0.2*thk_circ_r)
#       plt.arrow(thkx,thky,-thk_delx,-thk_dely,ec=None,fc='black',head_length=0.2*thk_circ_r,head_width=0.2*thk_circ_r)
        plt.plot(texs,[YSIN[-1]+cam2m*(i-XSIN[-1]) for i in texs],'black',linestyle='--',dashes=(5,4),linewidth=lw)
        plt.plot(texp,[YPIN[-1]+cam2m*(i-XPIN[-1]) for i in texp],'black',linestyle='--',dashes=(5,4),linewidth=lw)
        plt.plot(lexs,[YSIN[points*0.04]+cam1m*(i-XSIN[points*0.04]) for i in lexs],'black',linestyle='--',dashes=(5,4),linewidth=lw)
        plt.plot(lexp,[YPIN[po]+cam1m*(i-XPIN[po]) for i in lexp],'black',linestyle='--',dashes=(5,4),linewidth=lw)
#       plt.arrow(te_arr_bx,te_arr_by,te_arr_x,te_arr_y,ec=None,fc='black',head_length=0.5*TE_circ_r,head_width=0.5*TE_circ_r)
#       plt.arrow(te_arr_bx,te_arr_by,-te_arr_x,-te_arr_y,ec=None,fc='black',head_length=0.5*TE_circ_r,head_width=0.5*TE_circ_r)
#       plt.arrow(le_arr_bx,le_arr_by,-le_arr_x,le_arr_y,ec=None,fc='black',head_length=0.5*TE_circ_r,head_width=0.5*TE_circ_r)
#        plt.arrow(le_arr_bx,le_arr_by,le_arr_x,-le_arr_y,ec=None,fc='black',head_length=0.5*TE_circ_r,head_width=0.5*TE_circ_r)
        plt.plot(thkcx,thkcy,'black', linewidth=0.5)
        plt.plot(thkcx,thkcy2,'black', linewidth=0.5)
        plt.plot([XIN[cam_max]for i in cam_maxy],cam_maxy,'black',linestyle='--',dashes=(5,4),linewidth=lw)
#       plt.arrow((XIN[0]+XIN[cam_max])/2,YIN[0],(XIN[0]-XIN[cam_max])/2+arrows,0,ec=None,fc='black',head_length=arrows,head_width=arrows)
#       plt.arrow((XIN[0]+XIN[cam_max])/2,YIN[0],-(XIN[0]-XIN[cam_max])/2-arrows,0,ec=None,fc='black',head_length=arrows,head_width=arrows)
        plt.plot([np.amin(mins) for i in cx1],cx1,'black',linestyle='--',dashes=(5, 4),linewidth=lw)
        plt.plot([np.amax(maxes) for i in cx2],cx2,'black',linestyle='--',dashes=(5,4),linewidth=lw)
#       plt.arrow((endx+XIN[0])/2,YPIN[np.argmax(YPIN)]+(ptc+0.06)*cx,(endx-XIN[0])/2-arrows,0,ec=None,fc='black',head_length=arrows,head_width=arrows)
#       plt.arrow((endx+XIN[0])/2,YPIN[np.argmax(YPIN)]+(ptc+0.06)*cx,-(endx-XIN[0])/2+arrows,0,ec=None,fc='black',head_length=arrows,head_width=arrows)
        #plt.arrow(v1x,v1y,v1dx,v1dy,ec=None,fc='black',head_length=arrows,head_width=arrows)
        #plt.arrow(endx,endy+ptc*cx,v2dx,-v2dy,ec=None,fc='black',head_length=arrows,head_width=arrows)
#       plt.arrow(throat_bx,throat_by,throat_x,throat_y,ec=None,fc='black',head_length=arrows,head_width=arrows)
#        plt.arrow(throat_bx,throat_by,-throat_x,-throat_y,ec=None,fc='black',head_length=arrows,head_width=arrows)
        plt.plot(guide1x,guide1y,'black',linewidth=0.1)
        plt.plot(guide2x,guide2y,'black',linewidth=0.1)
        plt.plot(guide3x,guide3y,'black',linewidth=0.1)
        plt.plot(uncov_x,[YPIN[throat]+throat_m*(i-XPIN[throat]) for i in uncov_x],'black',linestyle='--',linewidth=lw)
        
        plt.plot(XPIN,YPIN,'black',linewidth=2)
        plt.plot(XSIN,YSIN,'black', linewidth=2)
        plt.plot(TEcx,TEcy,'black', linewidth=2)
        plt.plot(TEcx2,TEcy2,'black', linewidth=2)
        plt.plot(XIN,YIN,'black',linestyle='--', linewidth=1)
        plt.plot(XPIN,[i+ptc*cx for i in YPIN],'black',linewidth=2)
        plt.plot(XSIN,[i+ptc*cx for i in YSIN],'black', linewidth=2)
        plt.plot(TEcx,[i+ptc*cx for i in TEcy],'black', linewidth=2)
        plt.plot(TEcx2,[i+ptc*cx for i in TEcy2],'black', linewidth=2)
        plt.plot(XIN,[i+ptc*cx for i in YIN],'black',linestyle='--', linewidth=1)
        plt.axis('off')
    if 2==1:			### Plot of a stage
        
        vx = 0.4*cx
        edge_space = 0.15
        vx1 = np.arange(np.amin(mins)-vx-edge_space*cx,np.amin(mins)-edge_space*cx,0.01)
        vx2 = np.arange(np.amax(maxes)+edge_space*cx,np.amax(maxes)+vx+edge_space*cx,0.01)
        vx2 = np.arange(np.amax(maxes)+edge_space*cx,np.amax(maxes)+vx+edge_space*cx,0.01)
        vx3 = np.arange(np.amax(maxes)+edge_space*cx+cx*sx,np.amax(maxes)+vx+edge_space*cx+cx*sx,0.01)
        y_pos = YIN[0]+0.7*ptc*cx
        v1y = [y_pos + np.tan(-np.radians(X1))*(i-vx1[0]) for i in vx1]
        v2y = [y_pos + np.tan(-np.radians(X2))*(i-vx2[0]) for i in vx2]
        v3y = [y_pos + np.tan(-np.radians(X1))*(i-vx3[0]) for i in vx3]
        w2y = [y_pos + np.tan(np.radians(X1))*(i-vx2[0]) for i in vx2]
        w3y = [y_pos + np.tan(np.radians(X2))*(i-vx3[0]) for i in vx3]
        uy = np.arange(np.amax(YPIN)+sy+ptc*cx+0.05*cx,np.amax(YPIN)+sy+ptc*cx+0.05*cx+0.4*vx,0.01)
        ux = [XPIN[np.argmax(YPIN)]+cx*sx for i in uy]
        anghx = np.arange(np.amin(mins)+0.2*cx,np.amin(mins)+0.55*cx,0.01)
        anghy = [YIN[0]+sy+ptc*cx for i in anghx]
        angpx = np.arange(np.amin(mins)+0.2*cx,np.amin(mins)+0.45*cx,0.01)
        angpy = [anghy[0]+np.tan(np.radians(X1))*(i-anghx[0]) for i in angpx]
        angnx = np.arange(np.amin(mins)+0.2*cx,np.amin(mins)+0.45*cx,0.01)
        angny = [anghy[0]+np.tan(np.radians(-X1))*(i-anghx[0]) for i in angnx]
        
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

        plt.plot(XPIN,-YPIN,'black',linewidth=2)
        plt.plot(XSIN,-YSIN,'black', linewidth=2)
        plt.plot(TEcx,-TEcy,'black', linewidth=2)
        plt.plot(TEcx2,-TEcy2,'black', linewidth=2)
        plt.plot(XPIN,[i+ptc*cx for i in -YPIN],'black',linewidth=2)
        plt.plot(XSIN,[i+ptc*cx for i in -YSIN],'black', linewidth=2)
        plt.plot(TEcx,[i+ptc*cx for i in -TEcy],'black', linewidth=2)
        plt.plot(TEcx2,[i+ptc*cx for i in -TEcy2],'black', linewidth=2)

        
        plt.plot([i+sx*cx for i in XPIN],[i+sy for i in YPIN],'black',linewidth=2)
        plt.plot([i+sx*cx for i in XSIN],[i+sy for i in YSIN],'black', linewidth=2)
        plt.plot([i+sx*cx for i in TEcx],[i+sy for i in TEcy],'black', linewidth=2)
        plt.plot([i+sx*cx for i in TEcx2],[i+sy for i in TEcy2],'black', linewidth=2)
        plt.plot([i+sx*cx for i in XPIN],[i+ptc*cx+sy for i in YPIN],'black',linewidth=2)
        plt.plot([i+sx*cx for i in XSIN],[i+ptc*cx+sy for i in YSIN],'black', linewidth=2)
        plt.plot([i+sx*cx for i in TEcx],[i+ptc*cx+sy for i in TEcy],'black', linewidth=2)
        plt.plot([i+sx*cx for i in TEcx2],[i+ptc*cx+sy for i in TEcy2],'black', linewidth=2)
        plt.axis('off')

    delta = 0.5*abs(XSIN[-1]-XPIN[-1])		### RECALCULATE TRAILING EDGE THICKNESS
    return(XIN,YIN,XPIN,YPIN,XSIN,YSIN,TEcx,TEcy,TEcx2,TEcy2)

k = Profile2(40,-53,0.03)
plt.savefig('Stage.eps',format='eps',bbox_inches='tight',transparent=True)
plt.show()
a = np.transpose(np.asarray([k[6],k[7]]))
#np.savetxt("pr.csv", a, delimiter=",")

