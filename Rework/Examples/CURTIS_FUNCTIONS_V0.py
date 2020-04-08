#### code to perform an inverse design
import numpy as np
import time,math,numpy,copy
import random
import scipy.optimize as SciOpt
import matplotlib.pyplot as plt
import numpy.linalg as linalg
import os
import subprocess
from numpy import (atleast_1d, eye, mgrid, argmin, zeros, shape, squeeze,
vectorize, asarray, sqrt, Inf, asfarray, isinf)
#################################
### CLASSES USED ###
#####################

class DESIGN:				### INITIALISE A CLASS TO STORE GEOMETRY AND FLOW CLASSES
	def __init__(self):
		self.GEO=[]
class CLASS_GEO:			### GEOMETRY CLASS - CONTAINS GEOMETRIC INFORMATION OF PROFILE
    def __init__(self): 
        self.X = []
        self.Z =[]
        self.P2C = []
	#return(self)

class CLASS_FLOW:			### FLOW CLASS - CONTAINS AERODYNAMIC INFORMATION OF PROFILE
    def __init__(self): 
        self.S = []	### series of surface samples
        self.CM=[]	### pressure surface CP samples
	self.A1 =[]	### inlet angle -set as inlet condition
        self.A2 = []	### exit angle - set as cost function component
	self.M2 = []	### exit  mach number 
	self.M1S = []	### inlet isen mach number -used as mises input
	self.Gamma = 1.4	### flow gamma used for isentropic relations default to air
	self.Re2 = []	### prescribed exit mach number
	self.Re1 =[]	### inlet reynolds number used as mises inlet condition

######################################
### CUSTOM MISES IDAT READER - CJC ###
######################################

class IDAT():
	data=[]		### CLASS FILE TO STORE IDAT DATA

def fread(fid, nelements, dtype):
     if dtype is np.str:
         dt = np.uint8  # WARNING: assuming 8-bit ASCII for np.str!
     else:
         dt = dtype

     data_array = np.fromfile(fid, dt, nelements)
     data_array.shape = (nelements, 1)

     return data_array

def fortread(f,n,precision):
	head = fread(f,2,np.int16)
	data = fread(f,n,precision)
	foot = fread(f,2,np.int16)
	#if head!=foot:print 'header/footer mismatch'
	return(data)

def read_idat(filename):

	### initialise class file
	Idat =IDAT()

	### open file
	f= open(filename,'rb')

	### read headers
	head =fread(f,2,np.int16)
	head =fread(f,32,np.str)	### doesnt give text gives int values for strings-go figure
	head =fread(f,2,np.int16)
	### load idat
	[Idat.nstati,Idat.nstatr]=fortread(f,2,np.int32)
	Idat.istate = fortread(f,Idat.nstati,np.int32)
	Idat.rstate = fortread(f,Idat.nstatr,np.float64)
	
	### unpack useful data
	#print Idat.rstate
	Idat.xinl = Idat.rstate[48]
	Idat.xout = Idat.rstate[49]
	Idat.pitch = Idat.rstate[50]
	Idat.sinl = Idat.rstate[22]
	Idat.binl = np.arctan(Idat.sinl)*180/np.pi

	# Unpack useful indicies
	Idat.ii = Idat.istate[0]
	Idat.jj = Idat.istate[1]
	Idat.nbl = Idat.istate[2]
	#print Idat.nbl
	ns = Idat.nbl*2
	#print ns
	iih = Idat.istate[4]
	iip = Idat.istate[5]
	nbitn = Idat.istate[7]
	ngmode = Idat.istate[12]

	# Read mesh indicies
	store = fortread(f,4*ns,np.int32)
	store = np.reshape(store,[4,ns])
	Idat.jbld = store[0,:]
	Idat.ninl = store[1,:]
	Idat.nbld = store[2,:]
	Idat.nout = store[3,:]
	store = fortread(f,3.*Idat.nbl,np.int32)
	store = np.reshape(store,[3,Idat.nbl])
	Idat.iib = store[0,:]
	Idat.ible = store[1,:]
	Idat.nwak = store[2,:]


	# Read blade coordinates and inlet/ outlet grid spacing arrays
	Idat.xb=np.zeros((Idat.nbl,Idat.iib))
	Idat.yb=np.zeros((Idat.nbl,Idat.iib))
	Idat.xpb=np.zeros((Idat.nbl,Idat.iib))
	Idat.ypb=np.zeros((Idat.nbl,Idat.iib))
	Idat.sb=np.zeros((Idat.nbl,Idat.iib))
	Idat.sginl=np.zeros((Idat.nbl,Idat.ii))
	Idat.sgout=np.zeros((Idat.nbl,Idat.ii))
	Idat.xw=np.zeros((Idat.nbl,Idat.ii))
	Idat.yw=np.zeros((Idat.nbl,Idat.ii))
	Idat.wgap=np.zeros((Idat.nbl,Idat.ii))
	

	for ib in range(Idat.nbl):
		store = fortread(f,5*Idat.iib[ib],np.float64)
		store = np.reshape(store,[5,Idat.iib[ib]])
		Idat.xb[ib,:]  = store[0,:]
		Idat.yb[ib,:]  = store[1,:]
		Idat.xpb[ib,:] = store[2,:]
		Idat.ypb[ib,:] = store[3,:]
		Idat.sb[ib,:]  = store[4,:]
		store = fortread(f,5*Idat.ii,np.float64)
		store = np.reshape(store,[5,Idat.ii])
		Idat.sginl[ib,:] = store[0,:]
		Idat.sgout[ib,:] = store[1,:]
		Idat.xw[ib,:]    = store[2,:]
		Idat.yw[ib,:]    = store[3,:]
		Idat.wgap[ib,:]  = store[4,:]

	#print Idat.xb

	# Read streamline/streamtube/blade adjustment modes
	store = fortread(f,4,np.int32)
	Idat.nbvrx = store[0]
	Idat.nmodx = store[1]
	Idat.nparx = store[2]
	Idat.nmovx = store[3]
	if Idat.nbvrx > 0: Idat.bvrn = fortread(f,Idat.nbvrx,np.float64)
	if Idat.nmodx > 0: Idat.modn = fortread(f,Idat.nmodx,np.float64)
	if Idat.nparx > 0: Idat.parn = fortread(f,Idat.nparx,np.float64)
	if Idat.nmovx > 0: Idat.movn = fortread(f,Idat.nmovx,np.float64)

	# Read grid coordinates and mass fraction between streamlines
	Idat.mfrac=np.zeros((Idat.jj))
	Idat.x=np.zeros((Idat.ii,Idat.jj))
	Idat.y=np.zeros((Idat.ii,Idat.jj))
	Idat.r=np.zeros((Idat.ii,Idat.jj))
	for jj in range(Idat.jj):
		# Mass fraction between streamlines
		Idat.mfrac[jj] = fortread(f,1,np.float64)
		# Grid coordinates
		store = fortread(f,3*Idat.ii,np.float64)
		store = np.reshape(store,[Idat.ii,3])
		Idat.x[:,jj] = store[:,0]
		Idat.y[:,jj] = store[:,1]
		Idat.r[:,jj] = store[:,2]

	# Calculate xp for each surface (percent chord on blade surfaces), based on
	# the geometry of blade 1
	Idat.xp = np.zeros((Idat.ii,2))
	Idat.xp[:,0]=100.*((Idat.x[:,0])/(max(Idat.xb[0,:])-min(Idat.xb[0,:])))
	Idat.xp[:,1]=100.*((Idat.x[:,-1])/(max(Idat.xb[0,:])-min(Idat.xb[0,:])))


	# Read boundary layer and wall parameters
	Idat.sg=np.zeros((2,Idat.ii))
	Idat.disp=np.zeros((2,Idat.ii))
	Idat.sdisp=np.zeros((2,Idat.ii))
	Idat.pspec=np.zeros((2,Idat.ii))
	Idat.dhsdhb=np.zeros((2,Idat.ii))
	Idat.th=np.zeros((2,Idat.ii))
	Idat.dstr=np.zeros((2,Idat.ii))
	Idat.uedg=np.zeros((2,Idat.ii))
	Idat.ctau=np.zeros((2,Idat.ii))
	Idat.tau=np.zeros((2,Idat.ii))
	Idat.hwall=np.zeros((2,Idat.ii))
	Idat.qwall=np.zeros((2,Idat.ii))
	Idat.mwall=np.zeros((2,Idat.ii))
	Idat.psiw=np.zeros((2,Idat.ii))
	for isc in range(ns):
		store=fortread(f,5*Idat.ii,np.float64)
		store=np.reshape(store,[Idat.ii,5])
		Idat.sg[isc,:]=store[:,0]
		Idat.disp[isc,:]=store[:,1]
		Idat.sdisp[isc,:]=store[:,2]
		Idat.pspec[isc,:]=store[:,3]
		Idat.dhsdhb[isc,:]=store[:,4]
		store=fortread(f,5*Idat.ii,np.float64)
		store=np.reshape(store,[Idat.ii,5])
		Idat.th[isc,:]=store[:,0]
		Idat.dstr[isc,:]=store[:,1]
		Idat.uedg[isc,:]=store[:,2]
		Idat.ctau[isc,:]=store[:,3]
		Idat.tau[isc,:]=store[:,4]
		store=fortread(f,4*Idat.ii,np.float64).T
		store=np.reshape(store,[Idat.ii,4])
		Idat.hwall[isc,:]=store[:,0]
		Idat.qwall[isc,:]=store[:,1]
		Idat.mwall[isc,:]=store[:,2]
		Idat.psiw[isc,:]=store[:,3]

	#print Idat.sg
	# Read meridional and streamtube details
	npl=(3 + 2.*(Idat.nbvrx+1)); # number of elements "per line"
	store=fortread(f,npl*iih,np.float64)
	store=np.reshape(store,[npl,iih])
	Idat.xh=store[:,0]
	Idat.rh=store[:,1]
	#Idat.rph=store[:,2]
	iist=4
	iend=iist+Idat.nbvrx
	Idat.bh=store[:,iist:-1]
	iist=iend+1
	iend=iist+Idat.nbvrx
	Idat.bph=store[:,iist:iend]
	Idat.ibhdef=fortread(f,Idat.nbvrx,np.int32)

	# Read specified pressure loss
	store=fortread(f,3*iip,np.float64)
	store=np.reshape(store,[3,iip]).T
	Idat.xrl=store[:,0]
	Idat.prl=store[:,1]
	Idat.prlx=store[:,2]

	#Read forces and moments
	store=fortread(f,5*Idat.nbl,np.float64)
	store=np.reshape(store,[5,Idat.nbl]).T
	Idat.bldfx=store[:,0]
	Idat.bldfy=store[:,1]
	Idat.bldmz=store[:,2]
	Idat.bldvl=store[:,3]
	Idat.ppsgw=store[:,4]

	# Get spline properties
	store=fortread(f,2*Idat.nbl,np.float64)
	store=np.reshape(store,[2,Idat.nbl]).T
	Idat.sble=store[:,0]
	Idat.sblold=store[:,1]
	store=fortread(f,5*Idat.nbl,np.float64)
	store=np.reshape(store,[5,Idat.nbl]).T
	Idat.sblegn=store[:,0]
	Idat.xblegn=store[:,1]
	Idat.yblegn=store[:,2]
	Idat.sstg=store[:,3]
	Idat.swak=store[:,4]

	# Read transition and 2nd derivatives of pressures (left hand and right hand respectivley)
	store = fortread(f,3*isc,np.float64)
	store = np.reshape(store,[3,isc]).T
	Idat.xtr  = store[:,0]
	Idat.pxx0 = store[:,1]
	Idat.pxx1 = store[:,2]
	store = fortread(f,2*isc,np.int32)
	store = np.reshape(store,[2,isc]).T
	Idat.itran = store[:,0]
	Idat.ktran = store[:,1]




	diff_x = Idat.x[1:,:] - Idat.x[:-1,:]
	diff_y = Idat.y[1:,:] - Idat.y[:-1,:]


	Idat.alpha = (180./np.pi)*np.arctan(diff_y/diff_x)


	ile =0

	for i in range(len(Idat.x[:,0])):			### search for leading edge using grid perioidicty
		if Idat.y[i,0]==Idat.y[i,-1]-Idat.pitch:
			ile = i
		if Idat.y[i,0]!=Idat.y[i,-1]-Idat.pitch:break
	Idat.xle = Idat.x[ile,0]
	Idat.yle = Idat.y[ile,0]+Idat.pitch

	Idat.xStag = Idat.x[:ile+1,0]
	Idat.yStag = Idat.y[:ile+1,0]+Idat.pitch
	Idat.Ms = ((Idat.r**(-0.4)-1.)*(2./0.4))**0.5		### calculate isentropic mach numbers
	### decompose into surface distributions



	return(Idat)

def calc_secondaries(Idat):

	### calculate blade geometry and surface lengths
	Idat.xps = Idat.x[:,0]
	Idat.xss = Idat.x[:,-1]
	Idat.yps = Idat.y[:,0]
	Idat.yss = Idat.y[:,-1]-Idat.pitch

	plt.figure()
	plt.axis('equal')
	plt.plot(Idat.xss,Idat.yss,'-b')
	plt.plot(Idat.xps,Idat.yps,'-r')
	dist = 1.0
	for i in range(len(Idat.xps)-1):
		for i2 in range(len(Idat.xss)-1):
			if Idat.xps[i]==Idat.xss[i2] and Idat.yps[i]==Idat.yps[i2] and Idat.xps[i]<0.5:
				ile_ps=i
				ile_ss=i2
			
			if Idat.xps[i+1]>1.0 and Idat.xps[i]<1.0   and Idat.xss[i2+1]>1.0 and Idat.xss[i2]<1.0  :
				ite_ps=i
				ite_ss=i2
	
	plt.plot(Idat.xss[ile_ss],Idat.yss[ile_ss],'og')
	plt.plot(Idat.xss[ite_ss],Idat.yss[ite_ss],'og')			
	plt.plot(Idat.xps[ile_ps],Idat.yps[ile_ps],'og')
	plt.plot(Idat.xps[ite_ps],Idat.yps[ite_ps],'og')
	### calculate secondary parameters

	### based on euler solver node and cell based
	Idat.x_node = Idat.x
	Idat.y_node = Idat.y
	Idat.x_cell = 0.25*(Idat.x[:-1,:-1]+Idat.x[1:,:-1]+Idat.x[1:,1:]+Idat.x[:-1,1:])
	Idat.y_cell = 0.25*(Idat.y[:-1,:-1]+Idat.y[1:,:-1]+Idat.y[1:,1:]+Idat.y[:-1,1:])

	### this is ro/ro01 (euler). note this is NOT isentropic
	Idat.r_cell = Idat.r[1:-1,1:-1]
	Idat.Ro_Ro01_cell = Idat.r_cell

	### calculate SW direction
	diff_x = Idat.x[1:,:]-Idat.x[:-1,:]
	diff_y = Idat.y[1:,:]-Idat.y[:-1,:]
	Idat.alpha_x = 0.5*(Idat.x[1:,:]+Idat.x[:-1,:])
	Idat.alpha_y = 0.5*(Idat.y[1:,:]+Idat.y[:-1,:])
	Idat.alpha = np.degrees(np.arctan(diff_y/diff_x))
	Idat.alpha_cell = 0.5*(Idat.alpha[:,:-1]+Idat.alpha[:,1:])
	Idat.alpha_node = Idat.alpha

	### I-Planes -calculate angles at each height
	Idat.I_plane_grid_plane_angle = np.degrees(np.arctan((Idat.x[:,1:]-Idat.x[:,:-1])/(Idat.y[:,1:]-Idat.y[:,:-1])))



	return()

#######################################
###  GEOMETRY GEN FUNCTIONS	#######
#######################################

def F_Ki(Nm1,i):      ## function to generate bernstein polynomial
    Ki = math.factorial(Nm1)/(math.factorial(Nm1-i)*math.factorial(i))
    return(Ki)
    
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
        plt.figure(1)
        plt.title('shape function')
        plt.xlabel('psi')
        plt.ylabel('Shape Function: S(psi)')

    for i in range(N):	### iterate over each contirbution
        Ai=AIs[i]
        Ki = F_Ki(Nm1,i)		### calculate bernstein polynomial
        SC = F_SC(psi,i,Nm1)		### calculate shape component
        S[:]=S[:]+Ai*Ki*SC[:]

        if plotting==1:### another optional plot out
            plt.plot(psi,Ai*Ki*SC[:],'-',label='component'+str(i))
        
    if plotting==1:### another optional plot out
        plt.plot(psi,S,'-k',label='total of components')
	plt.legend(loc=0)
        #plt.show()  ### uncomment to plot imediatly     
    return(S)
    
def F_CONVERT_N(AIs,Nnew):### function solves linear algebra transform to increase/decrease order of parameter representation
    if Nnew<len(AIs):print 'WARNING convert doesnt work to reduce variables - not always possible'
    Nm1=Nnew-1
    psi_temp = np.linspace(0,1,Nnew) 	### only needs be as long as the parameters we are moving to
    S = F_SF(psi_temp,AIs,0)   		### generate the shape space at those points
    b=S					### output of transform
    A = np.zeros((Nnew,Nnew))
    for i in range(Nnew):
            A[i,:]=F_Ki(Nm1,i)*F_SC(psi_temp,i,Nm1)	### matrix of transform
            
    x = linalg.solve(A.T,b)		### calculate input of matrix using linalg
    #print AIs,x
    return(x)
    
def F_TF(psi,AIs,plotting,spacing): ### function that returns thickness contributions due to parameters

    if spacing == None:
    	S=F_SF(psi,AIs,0)	### calculate shape function
    else:
    	S=F_SF(spacing,AIs,0)	### calculate shape function
    C=F_CF(psi)		### calculate class function
    TF =S*C		### combine to form thickness funtion
    if plotting==1:	### optional plot out
        plt.figure(2)
        plt.title('thickness function')
        plt.xlabel('psi')
        plt.ylabel('thickness')
        plt.plot(psi,TF,'-r')
        plt.axis('equal')
    return(TF)
    
def F_NVT(psi,plotting): ### function that dictates distribution of thickness that is normal and component that is tangential

    portion = 2		### number of components of which the first only will be normal
    ANs = np.zeros(portion)	### initialise components as 0's
    ANs[0]=1.0			### set first component to 1
    #ANs[1]=1.0			### set first component to 1
    FN = F_SF(psi,ANs,0)  	### generate distribution (uses shape space to do so-why not!)
    #FN[:]=0.0			### option to force whole thickness to tangential
    #FN[:]=1.0			### option to force whole thickness to normal
    FT=1.0-FN			### tangential component =1-normal thickness
    if plotting==1:		### optional plot out
        plt.figure(3)
        plt.plot(psi,FN,'-b',label='normal')
        plt.plot(psi,FT,'-r',label='tangential')
        plt.legend(loc=0)
        plt.title('direction of thickness addition 1-normal,0-tangential')
        plt.xlabel('psi')
        plt.ylabel('fraction in applied in each direction')
    
    return(FN,FT)

def F_NVT2(psi,plotting): ### function that dictates distribution of thickness that is normal and component that is tangential

    portion = 3		### number of components of which the first only will be normal
    ANs = np.zeros(portion)	### initialise components as 0's
    ANs[0]=1.0			### set first component to 1
    #ANs[1]=1.0			### set first component to 1
    FN = F_SF(psi,ANs,0)  	### generate distribution (uses shape space to do so-why not!)
    #FN[:]=0.0			### option to force whole thickness to tangential
    #FN[:]=1.0			### option to force whole thickness to normal
    FT=1.0-FN			### tangential component =1-normal thickness
    if plotting==1:		### optional plot out
        plt.figure(3)
        plt.plot(psi,FN,'-b',label='normal')
        plt.plot(psi,FT,'-r',label='tangential')
        plt.legend(loc=0)
        plt.title('direction of thickness addition 1-normal,0-tangential')
        plt.xlabel('psi')
        plt.ylabel('fraction in applied in each direction')
    
    return(FN,FT)
    



def F_CAMBER3(psi,Xi1,Xi2,Gamma,plotting): ### NEW CAMBER LINE DEFINITION - provides smoother LE DROOP
	global FLOW_Target
	n = (np.tan(np.radians(Xi2))+np.tan(np.radians(Xi1)))/np.tan(np.radians(Gamma))
	a = np.tan(np.radians(Xi2))/n
	b= -np.tan(np.radians(Xi1))/n

        y_cam = a*psi**n+b*(1-psi)**n		### add camber comonent as a straight line at stagger angle -simplest version


        grad = np.zeros((len(psi),2))		### initialise camber gradient array
        norm = np.zeros((len(psi),2))		### initialise camber normals array
        s=np.zeros((len(psi)))			### initialise camber length array
        for i in range(1,len(psi)):
            s[i]=s[i-1]+((psi[i]-psi[i-1])**2+(y_cam[i]-y_cam[i-1])**2)**0.5		### calculate cumsum length of camber
        grad[:-1,0]=(psi[1:]-psi[:-1])/(s[1:]-s[:-1])					### calculate dx/ds
        grad[:-1,1]=(y_cam[1:]-y_cam[:-1])/(s[1:]-s[:-1])				### calculate dy/ds
        grad[-1,:]=grad[-2,:]								### deal with final point
        norm[:,0]=-grad[:,1]								### calculate normal 1 as -1/grad2
        norm[:,1]=grad[:,0]	
        if 0==1:	### optional plot out
            plt.figure(4)
            plt.title('Camber line - showing normal vectors')
            plt.xlabel('psi')
            plt.ylabel('zeta')
            plt.plot(psi,y_cam,'-b',label = 'camber')
            plt.axis('equal')
            for i in range(len(psi)):
                plt.plot([psi[i],psi[i]+norm[i,0]*0.01],[y_cam[i],y_cam[i]+norm[i,1]*0.01],'r')	### plot normal vectors
                plt.plot([psi[i],psi[i]-norm[i,0]*0.01],[y_cam[i],y_cam[i]-norm[i,1]*0.01],'r')
	    plt.show()
        return(y_cam,norm)


        
def Rescale(X,Z):  ### function restores true Axial Chord to 1 and leading edge offset to 0 (could be modified for max thickness offset)
    Xmin = X.min()
    Xmax = X.max()
    for i in range(len(X)):
        if X[i]==Xmin:
            #print 'ile found'   #- debug output         
            ile = i
        if X[i]==Xmax:
            #print 'ite found' 	#- debug output           
            ite = i
    Cx = Xmax-Xmin
    Z= Z/Cx
    X=(X-Xmin)/Cx
    YS = F_sampleSURF(X,Z,20)			### perform 20 point equispaced sample for average Y
    Zshift=np.mean(YS)				### calc average Y
    Z=Z-Zshift					### align based on average Y
    return(X,Z)
    
def Realign(X,Z): ### function reorientates the aerofoil so it starts at TE SS and goes around LE to PS (always assumes profile is clockwise )

    for i in range(len(X)):		### find ile and ite
        if X[i]==X.min():ile=i
        if X[i]==X.max():ite = i

    X_new = np.append(X[ite:],X[:ite])	### reorder to start at TE
    Z_new = np.append(Z[ite:],Z[:ite])
    #plt.figure(5)			### debug plotting to check working
    #plt.title('debug plot - realign')
    #plt.plot(X,Z,'-b')
    #plt.plot(X[0],Z[0],'go')
    #plt.plot(X_new[0],Z_new[0],'or')
    return(X_new,Z_new)
               
        
def F_Make(psi,AIs_upper,AIs_lower,Gamma,Xi1,Xi2,Tte,plotting):        ### function that compiles a blade from reduced variables 

    ### REDUCED VARIABLE DESCRIPTION:  #########
    ### AIs_upper - upper shape space parameters
    ### AIs_lower - lower shape space parameters
    ### Gamma - stagger angle (degrees)
    ### Xi1 - inlet metal angle
    ### Xi2 - exit metal angle
    ### Tte - trailing edge thickness
   
    Z_U=np.zeros((len(psi)))	### initialise X,Z for both upper and lower surfaces
    Z_L=np.zeros((len(psi)))
    X_U=np.zeros((len(psi)))
    X_L=np.zeros((len(psi)))

    (y_cam,norm) = F_CAMBER3(psi,Xi1,Xi2,Gamma,plotting)       ### generate camber line
    l_cam = np.zeros((len(y_cam)))
    yfake = np.zeros((len(y_cam)))
    nfactor = 4.0
    for i in range(len(y_cam)):
	if psi[i]<0.5:
		yfake[i] = (1./nfactor)*(2.*(psi[i]-0.5))**nfactor
    for i in range(1,len(y_cam)):
	l_cam[i]=l_cam[i-1]+((psi[i]-psi[i-1])**2.+(yfake[i]-yfake[i-1])**2)**0.5
    l_cam=l_cam/l_cam.max()
    

    (FN,FT) =F_NVT(psi,0)           ### generate application directions -last option enables plotting when =1
    (FNU,FTU) =F_NVT2(psi,0)           ### generate application directions -last option enables plotting when =1
    (TF_U) = F_TF(psi,AIs_upper,0,l_cam)  ### generate upper thickness
    (TF_L) = F_TF(psi,AIs_lower,0,None)  ### generate lower thickness
    
    Z_U = y_cam+FTU*TF_U+FNU*norm[0,1]*TF_U+np.cos(np.radians(Xi2))*Tte*0.5*(psi**0.5)	### apply thickness onto camber line for upper
    Z_L = y_cam+FT*TF_L+FN*norm[:,1]*TF_L-np.cos(np.radians(Xi2))*Tte*0.5*(psi**0.5)	### apply thickness onto camber line for lower
    X_U = psi+FNU*norm[0,0]*TF_U+np.sin(np.radians(-Xi2))*Tte*0.5*psi**(0.5)		### repeat upper for x rather than y
    X_L = psi+FN*norm[:,0]*TF_L-np.sin(np.radians(-Xi2))*Tte*0.5*psi**(0.5)		### repeat lower for x rather than y

    if plotting==1:		### optional plot out of complete blade
	print 'Gamma:',Gamma,' Xi1:',Xi1,' Xi2:', Xi2
        plt.figure(7)
        plt.title('Blade Profile - F_Make')
        plt.axis('equal')
        plt.plot(psi,y_cam,'-c',label = 'Construct line')
        plt.plot(X_U,Z_U,'-b', label = 'Upper')
        plt.plot(X_L,Z_L,'-r', label = 'Lower')
        #for i in range(len(psi)):				### optional plot of applied thickness vectors
         #   plt.plot([psi[i],X_U[i]],[y_cam[i],Z_U[i]],'-y')
         #   plt.plot([psi[i],X_L[i]],[y_cam[i],Z_L[i]],'-y')
	plt.legend(loc=0)
        plt.show()
    X = np.append(X_L[::-1],X_U[1:])	### combine axial distributions upper and lower
    Z = np.append(Z_L[::-1],Z_U[1:])	### combine tangential distributions upper and lower
    (X,Z) = Rescale(X,Z)		### perform rescale to ensure unit chord

    return(X,Z)

    
      
def Generate_DEF(definition):  		### Generate a profile from a definition vector (vector of optimisation variables) 

    #print definition			### optional print of definiton vector
    NAIs = (len(definition)-2)/2	### calculate number of upper and lower surface variables (assumed they are equal)
    psi=dist_vino(180,0,1,0.002,0.003)
    for i in range(40):
	psi[i] = psi[i]*(i+1)/40.
    Xi1 = definition[0]     		### INLET METAL	(degrees)
    Xi2 =  definition[1]     		### EXIT METAL	(degrees)
    Gamma=definition[2]     		### STAGGER ANGLE	(degrees)
    P2C = definition[3]     		### Pitch-to-chord (based on axial chord)
    Zte=0.025#definition[4]		### trailing edge thickness - HARDCODED IN V0
    Rle = definition[4]   		### LEADING EDGE RADIUS (as fraction of axial chord)
    Beta = definition[5]            	### TRAILING EDGE WEDGE ANGLE (degrees)

    ### PERFORM REDUCTION OF DEFINITION TO RAW PARAMETERS

    ### CC edit to use more parameters on SS than PS - configured for double the number on the SS as the PS
    NAIs = (len(definition)-2)/2

    NAI_PS = (len(definition)-6)/3+2
    NAI_SS = (len(definition)-6)*2/3+2
    #print definition
    #print NAI_PS,NAI_SS


    AIs_upper=np.zeros((NAI_PS))		### initialise upper shape space parameters
    AIs_lower=np.zeros((NAI_SS))		### intialise lower shape space paracmeters
    AIs_upper[0] = (2.*Rle)**0.5	### first variable based on leading edge radius (see kulfan)
    AIs_lower[0] = -(2.*Rle)**0.5
    if NAIs>2:
        AIs_upper[1:-1] = definition[6:6+NAI_PS-2]
        AIs_lower[1:-1] = definition[6+NAI_PS-2:6+(NAI_PS-2)+(NAI_SS-2)]

    ### set final variable based on wedge and metal angles - tangential addition

    AIs_upper[-1] = np.tan(np.radians(Gamma))-np.tan(np.radians(Xi2-Beta/2.))
    AIs_lower[-1] = np.tan(np.radians(Gamma))-np.tan(np.radians(Xi2+Beta/2.))

    AIs_upper[-1] = np.tan(np.radians(Xi2))-np.tan(np.radians(Xi2-Beta/2.))
    AIs_lower[-1] = np.tan(np.radians(Xi2))-np.tan(np.radians(Xi2+Beta/2.))	### replacement for new camber
    #AIs_lower[-1] = np.tan(np.radians(Gamma)-np.radians(Xi2+Beta/2.))#-np.tan(np.radians(Xi2+Beta/2.))   
    ### set final variables based on wedge and metal angles - normal addition
    #AIs_upper[-1] = np.tan(np.radians(Gamma-(Xi2-Beta/2.)))-np.tan(np.radians(Xi2-Beta/2.))	### final variable to set wedge and exit metal
    #AIs_lower[-1] = np.tan(np.radians(Gamma-(Xi2+Beta/2.)))-np.tan(np.radians(Xi2+Beta/2.))	### ditto


    (X,Z)=F_Make(psi,AIs_upper,AIs_lower,Gamma,Xi1,Xi2,Zte,0)	### PERFORM MAKE OF PROFILE
    
    #F_Check_Deriv(X,Z,1)		### OPTIONAL CHECK FOR DERIVATIVE DATA
    return(X,Z,P2C)


###############################
### GENERAL USE FUNCTIONS #####
###############################

    
def F_sampleSURF(X,Z,Nsamples):				### equispaced sampling of a curve
    Xskip = 0.00#5					### initial skip included for leading edge pressure distribution
    Xsample = np.linspace(Xskip,1,Nsamples+2)[1:-1]	### generate N+2 points and neglect first and last
    ile = len(X)/2					### initial guess at ile incase it fails
    for i in range(len(X)):				### find ile
        if X[i]==X.min():ile=i
    YU = np.interp(Xsample,X[ile:],Z[ile:])		### split into upper and lower surfaces and interp onto equispaced
    YL = np.interp(Xsample,X[:ile][::-1],Z[:ile][::-1])

    #plt.figure(6)					### optional plot out for debug
    #plt.title('debug plot - F_samplesurf')
    #plt.plot(X,Z,'-b',label='X-Z')
    #plt.plot(X[ile],Z[ile],'-om',label='LE')
    #plt.plot(X[0],Z[0],'-oc',label = 'X[0],Z[0]')
    #plt.plot(X[ile+1],Z[ile+1],'-oy',label = 'X]ile+1],Z[ile+1]')
    #plt.plot(Xsample,YU,'or',label='sample-upper')
    #plt.plot(Xsample,YL,'og',label='sample-lower')
    #plt.show()
    YS = np.append(YU,YL)		### recombine to two sampled surfaces for return
    
    return(np.append(Xsample,Xsample),YS)	### return with recombined X samples
    
def F_Check_Deriv(X,Z,plotting): ### function that checks derivatives of profile to ensure continuous curvature

    ### curvature of curve found based on K(s) = (x'(s)y''(s)-x''(s)y'(s))/(x'(s)^2+y'(s)^2)

    S=np.zeros((len(X)))
    Vectors = np.zeros((len(X)-1,2))
    S2 = np.zeros((len(X)-2))
    GRADX=np.zeros((len(X)-2))		### x'(s)
    GRADY=np.zeros((len(X)-2))		### y'(s)
    CURVEX=np.zeros((len(X)-2))		### x''(s)
    CURVEY = np.zeros((len(X)-2))	### y''(s)
    CURVE = np.zeros((len(X)-2))	

    for i in range(1,len(X)):
        S[i]=S[i-1]+((X[i]-X[i-1])**2+(Z[i]-Z[i-1])**2)**0.5			### old method of curvature using vectors
        if X[i]==X.min():ile=i

    ### second curvature definition --- better method

    for i in range(len(X)-2):
	if X[i]==X[i+2]:print 'error - repeat x[i]==x[i+2]'
	if X[i]==X[i+1]:print 'error - repeat x[i]=x[i+1]'
	if Z[i]==Z[i+2]:print 'error - repeat z[i]=z[i+1]'
	S2[i]=(S[i+2]-S[i])/2.			### arclength
	GRADX[i] = 0.5*(X[i+2]-X[i])/S2[i]	### x gradient with s
	GRADY[i] = 0.5*(Z[i+2]-Z[i])/S2[i]	### y gradient with s
	CURVEX[i] = (X[i+2]+X[i]-2*X[i+1])/(S2[i]**2)	### x CURVATURE with s
	CURVEY[i] = (Z[i+2]+Z[i]-2*Z[i+1])/(S2[i]**2)	### y CURVATURE with s

    TRUECURVE = (GRADX*CURVEY-CURVEX*GRADY)/(GRADX**2+GRADY**2)**0.5	### calculate true intrinsic curvature

    if plotting==1:	### plot derivatives --- leading edge anomaly is seen in calculation, based in small x' at leading edge
        plt.figure(8)
        plt.title('derivative data')
	plt.plot(X[1:],Z[1:],'-k',label = 'profile')
	plt.plot(X[ile:-1],TRUECURVE[ile-1:]/max(abs(TRUECURVE)),'-b',label = 'Curvature/Curvature.max()-PS')
	plt.plot(X[1:ile+1],TRUECURVE[:ile]/max(abs(TRUECURVE)),'-r',label = 'Curvature/Curvature.max()-SS')
        #plt.plot(X[1:-1],ICURVE,'--r',label = 'd2x/dy2')
        plt.legend()
	plt.axis('equal')

        plt.figure(9)
        plt.title('derivative data')
        plt.plot(S[1:-1]-S[ile],TRUECURVE,'-xb',label = 'Curvature')
	plt.plot(S[1:-1]-S[ile],X[1:-1],'-g',label = 'x')
	plt.plot(S[1:-1]-S[ile],GRADX[:],'-r',label = "x'(s)")
	plt.plot(S[1:-1]-S[ile],GRADY[:],'-c',label = "y'(s)")
	plt.plot(S[1:-1]-S[ile],CURVEX[:],'--r',label = "x''(s)")
	plt.plot(S[1:-1]-S[ile],CURVEY[:],'--c',label = "y''(s)")
	plt.plot([S[0],S[-1]]-S[ile],[0,0],'-k',label='0')
        plt.legend(loc=0)
        
	plt.show()
    return(S[1:-1:],TRUECURVE)
        
def F_CONVERT_PARAM(definition,Nincrease):	### wraps around linalg routine to scale up control parameters of surfaces
						### EDITED BY CC TO DO PERFORM INCREASE TO PS AND DOUBLE THAT TO SS
    NAIs = (len(definition)-2)/2
    psi=np.linspace(0,1,500)
    Xi1 = definition[0]     ### INLET METAL
    Xi2 = definition[1]     ### EXIT METAL
    Gamma= definition[2]     ### STAGGER ANGLE
    P2C = definition[3]     ### P2C
    Zte = 0.025#definition[4] 
    Rle = max(definition[4],0.01)   ### LEADING EDGE RADIUS
    Beta = definition[5]            ### TRAILING EDGE WEDGE ANGLE

    ### CC edit for lop sided parameter countes x2 on SS
    NAI_PS = (len(definition)-6)/3+2
    NAI_SS = (len(definition)-6)*2/3+2


    AIs_upper=np.zeros((NAI_PS))
    AIs_lower=np.zeros((NAI_SS))
    AIs_upper[0] = (2.*Rle)**0.5
    AIs_lower[0] = -(2.*Rle)**0.5  
    if NAIs>2:
        AIs_upper[1:-1] = definition[6:6+NAI_PS-2]
        AIs_lower[1:-1] = definition[6+NAI_PS-2:6+(NAI_PS-2)+(NAI_SS-2)]
    AIs_upper[-1] = np.tan(np.radians(Gamma))-np.tan(np.radians(Xi2-Beta/2.))
    AIs_lower[-1] = np.tan(np.radians(Gamma))-np.tan(np.radians(Xi2+Beta/2.))

    AIs_upper[-1] = np.tan(np.radians(Xi2))-np.tan(np.radians(Xi2-Beta/2.))
    AIs_lower[-1] = np.tan(np.radians(Xi2))-np.tan(np.radians(Xi2+Beta/2.))	### replacement for new camber



    Nnew=NAIs+Nincrease			### new number of parameters defining thickness on each surface
    NPSnew = NAI_PS+Nincrease
    NSSnew = NAI_SS+Nincrease*2

    ### UPPER SURFACE
    Nm1=NPSnew-1				### number minus one for bernstein poly.
    psi_temp = np.linspace(0,1,NPSnew) 	### only needs be as long as the parameters we are moving to
    S = F_SF(psi_temp,AIs_upper,0)   	### generate the shape space at those points
    b=S					### output of transform
    A = np.zeros((NPSnew,NPSnew))		### initialise Amatrix
    for i in range(NPSnew):
            A[i,:]=F_Ki(Nm1,i)*F_SC(psi_temp,i,Nm1)	### gen Amatrix
            
    x_upper = linalg.solve(A.T,b)	### solve linalg transform to match shape space at those points
    ### LOWER SURFACE
    Nm1=NSSnew-1				### number minus one for bernstein poly.
    psi_temp = np.linspace(0,1,NSSnew) 	### only needs be as long as the parameters we are moving to
    S = F_SF(psi_temp,AIs_lower,0)   	### generate the shape space at those points
    b=S					### output of transform
    A = np.zeros((NSSnew,NSSnew))
    for i in range(NSSnew):
            A[i,:]=F_Ki(Nm1,i)*F_SC(psi_temp,i,Nm1)
            
    x_lower = linalg.solve(A.T,b)
    
    PARAM = np.append(np.append(definition[:6],x_upper[1:-1]),x_lower[1:-1])	### APPEND TOGETHER

    return(PARAM)

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
    left = 0.00001
    delta = SciOpt.brentq(trans, left, 100)
    u = 0.5*(1.0 + np.tanh(delta*(eta/nm1 - 0.5))/(np.tanh(delta*0.5)))
    s = u/(a + (1-a)*u)
    return xst + (xen-xst)*s


########################################################
### DEFINITION OF ANALYTIC LOADING DISTRIBUTION -CJC ###
########################################################

def OvalARC(X,args):
	### FUNCTION FOR FITTING OVAL COMPONENTS OF LOADING DISTRIBUTION
	xc = X[0]
	yc = X[1]
	Rx = X[2]
	Ry = X[3]


	xl = args[0]
	yl = args[1]
	gradl = args[2]
	xr = args[3]
	yr = args[4]
	gradr = args[5]

	thiL = np.degrees(np.arccos((xc-xl)/Rx))
	thiR = np.degrees(np.arccos((xc-xr)/Rx))

	xl_temp = xc-Rx*np.cos(np.radians(thiL))
	yl_temp = yc+Ry*np.sin(np.radians(thiL))
	xr_temp = xc-Rx*np.cos(np.radians(thiR))
	yr_temp = yc+Ry*np.sin(np.radians(thiR))
	gradl_temp = min((Ry/Rx)/np.tan(np.radians(thiL)),10000)	
	gradr_temp = min((Ry/Rx)/np.tan(np.radians(thiR)),10000)

	#print 'temp:',xl_temp,xr_temp,yl_temp,yr_temp,gradl_temp,gradr_temp
	#print 'target:',xl,xr,yl,yr,gradl,gradr
	#print X
	cost = 0
	cost = cost + ((xl-xl_temp))**2.0
	cost = cost + ((yl-yl_temp))**2.0
	
	cost = cost + 0.1*((gradl-gradl_temp))**2.0
	cost = cost + 10.*((xr-xr_temp))**2.0
	cost = cost + 10.*((yr-yr_temp))**2.0
	cost = cost + 0.1*((gradr-gradr_temp))**2.0
	cost = cost*(1.+abs(Rx-Ry)/abs(Rx+Ry))


	return(cost*1.)



def Analytic_CM(Mpeak,Mless,Mps,lpeak,lsss,lsps,leps,gte):	
	### DEFINITION OF ANALYTIC LOADING PROFILE USING 8 CONTROL VARIABLES
	### -- MAGNITUDES --
	### MPEAK - M/MTE AT PEAK (BACK SURFACE DIFFUSION)
	### MLESS - EXTRAPOLATED LE LOADING
	### MPS   - FLAT PS MACH (STRATFORD DIFFUSION)
	### -- LOCATIONS --
	### LPEAK - PEAK MACH SURFACE LENGTH LOCATION
	### LSSS - START OF LINEAR SS RAMP
	### LSPS - START OF PRESSURE SURFACE HOLD
	### LEPS - END OF PRESSURE SURFACE HOLD
	### -- GRADIENT --
	### GTE - GRADIENT AT TRAILING EDGE OF SS


	xp_SS = lpeak
	Cmp_SS = Mpeak
	### INPUT XPSS - PEAK AXIAL LOCATION SS
	### CMP_SS -PEAK ND MACH NUMBER SS

	GRAD_le_SS = 1.1	### LE SS gradient
	CURVE_le_SS=0.0		### LE SS Curvature
	CURVE_peak1_SS =10.0	### Peak->LE SS Curvature
	CURVE_peak2_SS =10.0	### Peak->TE SS Curvature
	GRAD_te_SS = gte#2.0	### TE SS gradient
	CURVE_te_SS=3.0		### TE SS curvature

	xp_PS = leps#0.25		### peak axial location PS
	Cmp_PS = Mps#0.18	### peak ND MACH PS
	GRAD_le_PS = 5.0	### LE PS gradient
	CURVE_le_PS=18.0	### LE PS curvature
	CURVE_peak1_PS =0.0	### Peak->LE PS curvature
	CURVE_peak2_PS =0.0	### Peak->TE PS curvature
	GRAD_te_PS = 3.5 #3.0	### TE PS gradient
	CURVE_te_PS=12.5	
	    
	### suction surface LE->peak
	xp = xp_SS
	Cmp=Cmp_SS
	GRAD_le = GRAD_le_SS
	CURVE_le=CURVE_le_SS
	CURVE_peak =CURVE_peak1_SS

	ns = [2.0,2.5,3.0,3.5]
	MLE = Mless/Cmp
	b=np.zeros((4))
	b[0]=1.-MLE
	b[1]=GRAD_le*(1-MLE)
	b[2]=CURVE_le
	b[3]=CURVE_peak

	A = np.zeros((4,4))
	for i in range(4):
	    A[0,i]=1.0
	    A[1,i]=ns[i]
	    A[2,i]=ns[i]*(ns[i]-1)
	A[3,0]=ns[0]*(ns[0]-1)
		    
	x = linalg.solve(A,b)

	a1 = x[0]
	a2=x[1]
	a3=x[2]
	a4 = x[3]
	n1=ns[0]
	n2=ns[1]
	n3=ns[2]
	n4=ns[3]

	x1 = np.linspace(0,xp,500)
	Cm_lin1 = np.linspace(0,Cmp,len(x1))
	Cm_pow1 = MLE*Cmp+(a1*f_power((1-x1/xp),n1,Cmp,0)+ a2*f_power((1-x1/xp),n2,Cmp,0)+ a3*f_power((1-x1/xp),n3,Cmp,0)+ a4*f_power((1-x1/xp),n4,Cmp,0))


	### suction surface Peak->TE
	xp = xp_SS
	Cmp=Cmp_SS
	GRAD_le = GRAD_te_SS
	CURVE_le=CURVE_te_SS
	CURVE_peak =CURVE_peak2_SS

	ns = [2.0,2.5,3.0,3.5]

	b=np.zeros((4))
	b[0]=1.
	b[1]=GRAD_le
	b[2]=CURVE_le
	b[3]=CURVE_peak

	A = np.zeros((4,4))
	for i in range(4):
	    A[0,i]=1.0
	    A[1,i]=ns[i]
	    A[2,i]=ns[i]*(ns[i]-1)
	A[3,0]=ns[0]*(ns[0]-1)
		    
	x = linalg.solve(A,b)

	a1 = x[0]
	a2=x[1]
	a3=x[2]
	a4 = x[3]
	n1=ns[0]
	n2=ns[1]
	n3=ns[2]
	n4=ns[3]
	x2 = np.linspace(xp,1,500)
	Cm_lin2 = np.linspace(Cmp,1.0,len(x2))
	Cm_pow2 = a1*f_power((x2-xp)/(1-xp),n1,(Cmp-1),1.)+a2*f_power((x2-xp)/(1-xp),n2,(Cmp-1),1.)+a3*f_power((x2-xp)/(1-xp),n3,(Cmp-1),1.)+a4*f_power((x2-xp)/(1-xp),n4,(Cmp-1),1.)



	### pressure surface LE->peak
	xp = xp_PS
	Cmp=Cmp_PS
	GRAD_le = GRAD_le_PS
	CURVE_le=CURVE_le_PS
	CURVE_peak =CURVE_peak1_PS

	ns = [2.0,2.5,3.0,3.5]

	b=np.zeros((4))
	b[0]=1.
	b[1]=GRAD_le
	b[2]=CURVE_le
	b[3]=CURVE_peak

	A = np.zeros((4,4))
	for i in range(4):
	    A[0,i]=1.0
	    A[1,i]=ns[i]
	    A[2,i]=ns[i]*(ns[i]-1)
	A[3,0]=ns[0]*(ns[0]-1)
		    
	x = linalg.solve(A,b)

	a1 = x[0]
	a2=x[1]
	a3=x[2]
	a4 = x[3]
	n1=ns[0]
	n2=ns[1]
	n3=ns[2]
	n4=ns[3]

	x3 = np.linspace(0,xp,500)
	Cm_lin3 = np.linspace(0,Cmp,len(x3))
	Cm_pow3 = (a1*f_power((1-x3/xp),n1,Cmp,0)+ a2*f_power((1-x3/xp),n2,Cmp,0)+ a3*f_power((1-x3/xp),n3,Cmp,0)+ a4*f_power((1-x3/xp),n4,Cmp,0))
	Cm_pow3 = np.linspace(Cmp,Cmp,len(x3))

	xp = xp_PS
	Cmp=Cmp_PS
	GRAD_le = GRAD_te_PS
	CURVE_le=CURVE_te_PS
	CURVE_peak =CURVE_peak2_PS

	ns = [2.0,2.5,3.0,3.5]

	b=np.zeros((4))
	b[0]=1.
	b[1]=GRAD_le
	b[2]=CURVE_le
	b[3]=CURVE_peak

	A = np.zeros((4,4))
	for i in range(4):
	    A[0,i]=1.0
	    A[1,i]=ns[i]
	    A[2,i]=ns[i]*(ns[i]-1)
	A[3,0]=ns[0]*(ns[0]-1)
		    
	x = linalg.solve(A,b)

	a1 = x[0]
	a2=x[1]
	a3=x[2]
	a4 = x[3]
	n1=ns[0]
	n2=ns[1]
	n3=ns[2]
	n4=ns[3]
	x4 = np.linspace(xp,1,500)
	Cm_lin4 = np.linspace(Cmp,1.0,len(x4))
	Cm_pow4 = a1*f_power((x4-xp)/(1-xp),n1,(Cmp-1),1.)+a2*f_power((x4-xp)/(1-xp),n2,(Cmp-1),1.)+a3*f_power((x4-xp)/(1-xp),n3,(Cmp-1),1.)+a4*f_power((x4-xp)/(1-xp),n4,(Cmp-1),1.)

	if 1==0:
		plt.figure(1)
		plt.plot(x1,Cm_lin1,'--b')
		plt.plot(x2,Cm_lin2,'--b')
		plt.plot(x3,Cm_lin3,'--b')
		plt.plot(x4,Cm_lin4,'--b')
		plt.plot(x1,Cm_pow1,'--r')
		plt.plot(x2,Cm_pow2,'--r')
		plt.plot(x3,Cm_pow3,'--r')
		plt.plot(x4,Cm_pow4,'--r')


    	XSS = np.append(x1,x2[1:])
    	XPS = np.append(x3,x4[1:])
    	CMSS = np.append(Cm_pow1,Cm_pow2[1:])
    	CMPS = np.append(Cm_pow3,Cm_pow4[1:])

	xskipPS= lsps
	xskipSS= lsss
  	for i in range(len(XSS)):
			if XSS[i]<xskipSS:continue
			else:
				yleloadSS=CMSS[i]
				break
	### oval arc fit
	xl = 0
	yl = 0

	xr = xskipSS
	yr = np.interp(xskipSS,XSS,CMSS)
	gradr = (np.interp(xskipSS+0.001,XSS,CMSS)-yr)/0.001
	gradl = 10.*(yr-yl)/(xr-xl)
	X0 = [xr,yl,xr,yr]
	ARGS = [xl,yl,gradl,xr,yr,gradr]

    	Result = SciOpt.fmin(OvalARC,X0,args=[ARGS],maxiter=50000000,disp=False,ftol=0.000000001)

	xc=Result[0]
	yc=Result[1]
	Rx=Result[2]
	Ry=Result[3]
	thil = np.degrees(np.arccos((xc-xl)/Rx))
	thir = np.degrees(np.arccos((xc-xr)/Rx))
	Tarc = np.linspace(thil+0.00001,thir-0.0001,100)
	Xarc = xc-Rx*np.cos(np.radians(Tarc))
	Yarc = yc+Ry*np.sin(np.radians(Tarc))


	for i in range(len(XSS)):
			if XSS[i]<xskipSS:
				CMSS[i] = np.interp(XSS[i],Xarc,Yarc)



  	for i in range(len(XPS)):
			if XPS[i]<xskipPS:continue
			else:
				yleloadPS=CMPS[i]
				break

	### oval arc fit
	xl = 0
	yl = 0

	xr = xskipPS
	yr = np.interp(xskipPS,XPS,CMPS)
	gradr = (np.interp(xskipPS+0.001,XPS,CMPS)-yr)/0.001

	gradl = 10.*(yr-yl)/(xr-xl)

	X0 = [xr*1.1,yl-yr*0.1,xr*1.1,yr*1.1]
	ARGS = [xl,yl,gradl,xr,yr,gradr]
    	Result = SciOpt.fmin(OvalARC,X0,args=[ARGS],maxiter=500000,disp=False,ftol=0.000000001)
	xc=Result[0]
	yc=Result[1]
	Rx=Result[2]
	Ry=Result[3]
	thil = np.degrees(np.arccos((xc-xl)/Rx))
	thir = np.degrees(np.arccos((xc-xr)/Rx))
	Tarc = np.linspace(thil+0.00001,thir-0.0001,100)
	Xarc = xc-Rx*np.cos(np.radians(Tarc))
	Yarc = yc+Ry*np.sin(np.radians(Tarc))


	for i in range(len(XPS)):
			if XPS[i]<xskipPS:
				CMPS[i]=np.interp(XPS[i],Xarc,Yarc)


    	X = np.append(XSS[::-1] ,XPS[1:])   
    	CM = np.append(CMSS[::-1] ,CMPS[1:]) 

    	return(X,CM)


def f_power(x,n,a,o):	### POWER LAW FUNCTION
    
    return(a*(1-x**n)+o)


def M2_to_M1(Gamma,M2,A1,A2):	### gives isentropic inlet mach based on flow angles and exit mach
	### input description
	### gamma - air gamma
	### M2 - exit mach
	### M1 - inlet mach
	### A1 - Inlet angle
	### A2 - Exit angle

	M1 = 0.1	### initial guess

	A2_A1 = np.cos(np.radians(A2))/np.cos(np.radians(A1))	### area ratio
	FM2 =  M2*(1.+0.5*(Gamma-1)*M2**2)**(-0.5*(Gamma+1)/(Gamma-1))
	
	### iterate with shrinking step size
	for i in range(10):
		FM1 =  M1*(1.+0.5*(Gamma-1)*M1**2)**(-0.5*(Gamma+1)/(Gamma-1))
		if FM1>A2_A1*FM2:
			M1=M1-0.1
			break
		else:	M1=M1+0.1
	for i in range(10):
		FM1 =  M1*(1.+0.5*(Gamma-1)*M1**2)**(-0.5*(Gamma+1)/(Gamma-1))
		if FM1>A2_A1*FM2:
			M1=M1-0.01
			break
		else:	M1=M1+0.01
	for i in range(10):
		FM1 =  M1*(1.+0.5*(Gamma-1)*M1**2)**(-0.5*(Gamma+1)/(Gamma-1))
		if FM1>A2_A1*FM2:
			M1=M1-0.001
			break
		else:	M1=M1+0.001

	for i in range(10):
		FM1 =  M1*(1.+0.5*(Gamma-1)*M1**2)**(-0.5*(Gamma+1)/(Gamma-1))
		if FM1>A2_A1*FM2:
			M1=M1-0.0001
			break
		else:	M1=M1+0.0001

	return(M1)

def RE2_to_RE1(RE2,A1,A2):	### assume constant viscosity and use axial chord
				### uses continuity to relate ro1v1 to ro2v2
	A2_A1 = np.cos(np.radians(A2))/np.cos(np.radians(A1))	### area ratio
	RE1 = RE2*A2_A1	
	### needs updating with viscosity variation due to temperature 
	
	return(RE1)

def Read_Mach(fname,Flow_Target):	### reads mach number distribution from output of edp (to be updated with binary read data at some point)

	f=open(fname)
	lines =f.readlines()
	read = 0
	S_SS=[]
	S_PS=[]
	M_SS=[]
	M_PS=[]
	for line in lines:
		if '! sig(upper)  Mach' in line: 
			read=1
			continue
		if '! sig(lower)  Mach' in line: 
			read=2
			continue
		if 'End' in line:break
		if read==1:
			S_PS.append(float(line.split()[0]))
			M_PS.append(float(line.split()[1]))
		if read==2:
			S_SS.append(float(line.split()[0]))
			M_SS.append(float(line.split()[1]))
	S_SS = np.array(S_SS)
	S_PS = np.array(S_PS)
	M_SS = np.array(M_SS)
	M_PS = np.array(M_PS)
	M_TE = 0.5*(M_SS[-1]+M_PS[-1])

	CM_SS=M_SS/M_SS[-1]#Flow_Target.M2#M_SS[-1]		### non-dimensionalise with final value on surface
	CM_PS=M_PS/M_PS[-1]#Flow_Target.M2#M_PS[-1]		### ditto
	CM_SS=CM_SS[1:]
	CM_PS=CM_PS[1:]	
	S_SS=S_SS[1:]
	S_PS=S_PS[1:]
	S_SS = (S_SS-S_SS.min())/(S_SS.max()-S_SS.min())
	S_PS = (S_PS-S_PS.min())/(S_PS.max()-S_PS.min())
	#print CM_PS
	S = np.append(S_SS[::-1] ,S_PS[1:])   			### combine surface distance
    	CM = np.append(CM_SS[::-1] ,CM_PS[1:])   		### combine mach coefficients
	M = np.append(M_SS[::-1] ,M_PS[1:]) 
	return(S,CM,M,M_TE)

def Read_YAW(fname):			### function that reads yaw angle from output file produced by iplot
	f=open(fname)
	lines=f.readlines()
	YP=0.0000001
	A2=0.0
	M2=0.0
	for line in lines:
		if 'S2' in line:
			if line.split()[4]=='NaN':
				A2='NaN'
			else:
				if line.split()[4] == 'deg.':
					A2 = float(line.split()[3])
				else:
					A2 = float(line.split()[4])
		if 'M1(ref)' in line:
			if line.split()[2]=='NaN':
				M1='NaN'
			else:

				M1 = float(line.split()[2])
		if 'M2(ref)' in line:
			if line.split()[2]=='NaN':
				#print line.split()
				M2='NaN'
			else:

				M2 = float(line.split()[2])

		#YP=0.0000001
		if 'Zeta' in line:
			#print 'viscous loss'
			if line.split()[2]=='NaN':
				YP='NaN'
			else:

				YP = float(line.split()[2])
		#print YP
	#print M2
	return(A2,M2,YP,M1)


def load(fname): 			### load external geometry listed as x,z , must be clockwise TE->TE

    f=open(fname,'r')
    lines=f.readlines()
    XPROF=np.zeros((len(lines)))
    YPROF=np.zeros((len(lines)))
    
    for i in range(len(lines)):
    		XPROF[i]=float(lines[i].split()[0])
    		YPROF[i]=float(lines[i].split()[1])
    for i in range(len(lines)):
    		if XPROF[i]==XPROF.min():
    			ile = i
    Cx = XPROF.max()-XPROF.min()
    YPROF2=(YPROF-YPROF[ile])/Cx
    XPROF2=(XPROF-XPROF[ile])/Cx
    (X,Z)=Realign(XPROF2,YPROF2)		### realign to center 
    (X,Z)=Rescale(X[::-1],-Z[::-1])  		### rescale to unit chord
    
    return(X,Z)

def run_command(command):	### function to run command in command line
    p = subprocess.Popen(command,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT)
    return iter(p.stdout.readline, b'')

################################################################
### FUNCTIONS USED TO GENERATE COST FUNCTION AND CONSTRAINTS ###
################################################################


def GAUS(MAG,CENT,WIDTH,x):				### GAUSSIAN FUNCTION
	f=MAG*np.exp(-(x-CENT)**2/(2*WIDTH**2))
	return(f)







def SURFLENGTHS(GEO,FLOW)	:
	### FUNCTION THAT PERFORMS SOME ANALYSIS OF THE PROFILE LENGTHS AND PEAK LOCATIONS ETC
	mini = 1000.
	for i in range(len(GEO.X)):
		if (GEO.X[i]-FLOW.IDAT.xle)**2+(GEO.Z[i]-FLOW.IDAT.yle)**2<mini:
			mini=(GEO.X[i]-FLOW.IDAT.xle)**2+(GEO.Z[i]-FLOW.IDAT.yle)**2
			GEO.ILE_FLOW=i

	for i in range(len(FLOW.S)):
		if FLOW.S[i]==0.0:
			iSLE=i

			FLOW.CM[i]=0.0
			break

	for i in range(len(FLOW.S)):
		if FLOW.CM[i]==FLOW.CM.max():
			FLOW.IPEAK = i
			break

	FLOW.iSLE=iSLE
	FLOW.ILE_FLOW = iSLE			### streamline based leading edge
	nSS = GEO.ILE_FLOW			### streamline based leading edge
	nPS = len(GEO.X)-GEO.ILE_FLOW
	LPS = 0
	LSS = 0
	for i in range(nSS-1):
		LSS=LSS+((GEO.X[i+1]-GEO.X[i])**2+(GEO.Z[i+1]-GEO.Z[i])**2)**0.5
	for i in range(nPS-1):
		LPS=LPS+((GEO.X[nSS+i+1]-GEO.X[nSS+i])**2+(GEO.Z[nSS+i+1]-GEO.Z[nSS+i])**2)**0.5


	GEO.LSS = LSS			### SS length
	GEO.LPS = LPS			### PS length
	GEO.LTOT = LSS+LPS		### total length
	### form surface mach and pressure distributions
	xps = FLOW.IDAT.x[:,0]
	Mps = FLOW.IDAT.Ms[:,0]
	xss = FLOW.IDAT.x[:,-2]
	Mss = FLOW.IDAT.Ms[:,-2]
	FLOW.MS = np.zeros((len(GEO.X)))
	FLOW.S2  = np.zeros((len(GEO.X)))	#### to put in distance
	FLOW.X  = np.zeros((len(FLOW.S)))
	FLOW.Z  = np.zeros((len(FLOW.S)))

	for i in range(GEO.ILE_FLOW):
		FLOW.S2[GEO.ILE_FLOW-i-1]=FLOW.S2[GEO.ILE_FLOW-i]+((GEO.X[GEO.ILE_FLOW-i-1]-GEO.X[GEO.ILE_FLOW-i])**2+(GEO.Z[GEO.ILE_FLOW-i-1]-GEO.Z[GEO.ILE_FLOW-i])**2)**0.5

	FLOW.S2[:GEO.ILE_FLOW]=FLOW.S2[:GEO.ILE_FLOW]/FLOW.S2[:GEO.ILE_FLOW].max()

	for i in range(GEO.ILE_FLOW,len(GEO.X)-1):
		FLOW.S2[i+1]=FLOW.S2[i]+((GEO.X[i+1]-GEO.X[i])**2+(GEO.Z[i+1]-GEO.Z[i])**2)**0.5

	FLOW.S2[GEO.ILE_FLOW:]=FLOW.S2[GEO.ILE_FLOW:]/FLOW.S2[GEO.ILE_FLOW:].max()

	GEO.S = FLOW.S2
	side='SS'
        for i in range(len(FLOW.S)):
		if side == 'SS' and i<iSLE:
			FLOW.X[i]=np.interp(FLOW.S[i],GEO.S[:GEO.ILE_FLOW][::-1],GEO.X[:GEO.ILE_FLOW][::-1])
			FLOW.Z[i]=np.interp(FLOW.S[i],GEO.S[:GEO.ILE_FLOW][::-1],GEO.Z[:GEO.ILE_FLOW][::-1])
		else:
			side='PS'
			FLOW.X[i]=np.interp(FLOW.S[i],GEO.S[GEO.ILE_FLOW:],GEO.X[GEO.ILE_FLOW:])
			FLOW.Z[i]=np.interp(FLOW.S[i],GEO.S[GEO.ILE_FLOW:],GEO.Z[GEO.ILE_FLOW:])
		

	side = 'SS'
	for i in range(len(GEO.X)):
		if side == 'SS' and i<GEO.ILE_FLOW:
			FLOW.MS[i]=np.interp(FLOW.S2[i],FLOW.S[:iSLE][::-1],FLOW.CM[:iSLE][::-1])
		else:
			side='PS'
			FLOW.MS[i]=np.interp(FLOW.S2[i],FLOW.S[FLOW.ILE_FLOW:],FLOW.CM[iSLE:])

	FLOW.RP2 = (1.+(0.4/2.0)*(FLOW.MS)**2)**(-1.4/0.2)
	FLOW.CM2=FLOW.MS/FLOW.MS[0]
	FLOW.CP2 = (1.-FLOW.RP2)/(1.-FLOW.RP2[0])
	#plt.figure()
	#plt.plot(GEO.X,FLOW.MS,'-xb')
	#plt.plot(GEO.X[0],FLOW.MS[0],'-ob')
	#plt.plot(FLOW.S[:],FLOW.CM[:],'-r')

	#plt.show()

	#print 'PSlength:',LPS,'SSlength:',LSS

	#plt.plot(GEO.X,GEO.Z,'-k')
	#plt.plot(GEO.X[GEO.ILE_FLOW],GEO.Z[GEO.ILE_FLOW],'or')
	#plt.plot(FLOW.IDAT.xle,FLOW.IDAT.yle,'og')
	#plt.show()
	return(GEO,FLOW)




######################################
### MODULES FOR INVERSE BUILD UP #####
######################################

### MACH NUMBER DISTRIBUTION FIT

### PARAMETERS TO GEOMETRY      ->      PARAM_TO_GEO
### GEOMETRY TO FLOW 		->	GEO_TO_FLOW
### FLOW TO COST FUNCTION	-> 	FLOW_TO_COST
### WRAPPER			-> 	PARAM_TO_COST_VIA_FLOW

### GEOMETRY FIT

### PARAMETERS TO GEOMETRY      ->      PARAM_TO_GEO
### GEOMETRY TO COST 		->	GEO_TO_FLOW
### WRAPPER			-> 	PARAM_TO_COST_VIA_GEO

def PARAM_TO_GEO(PARAM):	### MODULE - LINK PARAMETERS TO GEOMETRY DEFINITION - PARAMS ARE THOSE VARIED BY OPTIMIZER

    GEO = CLASS_GEO()					### Geometry class
    (GEO.X,GEO.Z,GEO.P2C)= Generate_DEF(PARAM)		### PARAM->Geometry module (default is CCfoils)
    
    return(GEO)			


#######################################
### FUNCTIONS TO SET UP MISES #########
#######################################
    
def WRITE_MISES_BLADE(GEO,FLOW,fout): ### WRITE MISES BLADE.CCFOIL

    f=open(fout,'w')
    f.write('CCfoil\n')
    f.write(str(np.tan(np.radians(FLOW.A1)))+' '+str(np.tan(np.radians(FLOW.A1)))+' 0.5 0.5 '+str(GEO.P2C)+'\n')		### WRITE ARBITARY FLOW ANGLES AND CORRECT PITCH TO CHORD
    for i in range(len(GEO.X)-1):
        f.write(str(GEO.X[i])+' '+str(GEO.Z[i])+' \n')			### WRITE SURFACE PROFILE IN X-Z
    f.write(str(GEO.X[-1])+' '+str(GEO.Z[-1]))
    f.close()
    #plt.figure()
    #plt.plot(GEO.X,GEO.Z,'-xk')
    #plt.show()

def WRITE_MISES_ISES(FLOW,fout): ### WRITE ISES FILE FOR MISES - CONTAINS VISCOUS INFORMATION AND TRANSITIONS

    #print 'written m1s:',FLOW.M1S
    f=open(fout,'w')
    f.write('1 2 5 6 !15\n')	### MISES VARIABLES - DEFAULT SETUP 1 2 5 6 -> TRANSONIC FIXED INLET MACH NUMBER AND INLET ANGLE
    f.write('1 4 3 6 !17 \n')
    #f.write('1 2 5 15 \n')	### MISES VARIABLES - DEFAULT SETUP 1 2 5 6 -> TRANSONIC FIXED INLET MACH NUMBER AND INLET ANGLE
    #f.write('1 4 3 17\n')	### CONSTRAINTS -DEFAULT SETUP 1 4 3 6 -> TRANSONIC FIXED INLET MACH NUMBER AND INLET ANGLE (SEE MISES MANUAL FOR ALTERNATIVES)
    f.write(str(FLOW.M1S)+' '+'0.0'+' '+str(np.tan(np.radians(FLOW.A1)))+' '+'-0.4 \n')	### INLET MACH, FLOW ANGLE AND INLET PLANE
    f.write(str(FLOW.M2)+' '+'0.0'+' '+str(np.tan(np.radians(FLOW.A2)))+' '+'1.4 \n')	### EXIT MACH FLOW ANGLE AND EXIT PLANE
    f.write('0.0'+' '+'0.0'+' '+str(FLOW.Gamma)+' '+'0.4 \n'	)			### GAS PROPETIES
    f.write(str(FLOW.Re1)+' '+'-0.5'+' '+'0.1'+' '+'0.1 \n'	)			### REYNOLDS NUMBER BASED ON INLET MACH
    f.write('0.3 0.3 \n')								### TRANSITION LOCATIONS MUST BE >0
    f.write('4    0.95    1.0    0.0000 \n')
    f.write('0  0')
    f.close()

def WRITE_MISES_ISES_INV(FLOW,fout): ### WRITE ISES FILE FOR MISES BUT SET INVISCID AND TO A FIXED EXIT MAHC NUMBER

    f=open(fout,'w')
    f.write('1 2 5 6 15 \n')	### VARIABLES -DEFUALT(1 2 5 6 15) FIXES INLET ANGLE AND EXIT MACH
    f.write('1 4 3 6 17\n')	### CONSTRAINTS - DEFAULT(1 4 3 6 17) FIXES INLET ANGLE AND EXIT MACH
    f.write(str(FLOW.M1S)+' '+'1.0'+' '+str(np.tan(np.radians(FLOW.A1)))+' '+'-0.5 \n')	### inlet
    f.write(str(FLOW.M2)+' '+'1.0'+' '+str(np.tan(np.radians(FLOW.A2)))+' '+'1.0 \n')	### exit
    f.write('0.0'+' '+'0.0'+' '+str(FLOW.Gamma)+' '+'0.4 \n'	)
    f.write(str(0.0)+' '+'-0.5'+' '+'0.1'+' '+'0.1 \n'	)			### DFAULTS REYNOLDS TO 0 HENCE INVISCID
    f.write('1.1 1.1 \n')
    f.write('4    0.95    1.0    0.0000 \n')
    f.write('0  0')
    f.close()


def RUN_MISES_TO_M2(FLOW_Target,GEO):
	#########################################################################################
	### WORK AROUND FOR OLD VERSION OF MISES THAT IS NOT STABLE TO FIXED EXIT MACH NUMBER ###
	#########################################################################################
	global PARAM_OLD  ,M1S_GLOBAL 

	### USER INPUT - DELAY FOR READ WRITE DELAYS IF THEY EXIST
	dT = 0.01	### time gaps for read write

	### insert code to solve geometry and resolve flow	#####
	FLOW=CLASS_FLOW()
	FLOW_temp=copy.deepcopy(FLOW_Target)

	WRITE_MISES_BLADE(GEO,FLOW_temp,'blade.CCfoil')	### write blade to MISES input
	fout='ises.CCfoil'
	WRITE_MISES_ISES(FLOW_temp,fout)		### WRITE ISES FILE SETTING UP VISCID SETUP

	time.sleep(dT)	### pause after file write
	os.system('sh run_process.txt > tmp_all')	### RUN initial mises run using shell script
	os.system('sh run_edp_proc.txt >OUTPUT')	### RUN edp process using shell script

	### iterate inlet mach to achieve exit mach
	error = 1.0	### default error to 1
	it = 0		### set counter to 0
	fname='OUTPUT'
	(FLOW.A2,FLOW.M2,FLOW.YP,FLOW.M1S)=Read_YAW(fname)	### read flow condition

	while error > 0.005 and it<30:		### iterate inlet mach number
			it=it+1
			fname='OUTPUT'
			
			if FLOW.A2 == 'NaN':			### break loop in diverging case
				print 'yaw angle: NAN'
				break
			if FLOW.M2 == 'NaN':break		### break loop in diverging case

			FLOW_temp.M1S = FLOW_temp.M1S*(FLOW_Target.M2/FLOW.M2)**0.8	### update inlet mach #
			
			WRITE_MISES_ISES(FLOW_temp,fout)		### WRITE ISES FILE SETTING UP VISCID SETUP
			time.sleep(dT)
			os.system('sh run_ises_single.txt > tmp_all')		### run at new mach with BL iteration
			time.sleep(dT)
			
			os.system('sh run_edp_proc.txt >OUTPUT')		### RUN edp process using shell script
			time.sleep(dT)
			fname='OUTPUT'

			(FLOW.A2,FLOW.M2,FLOW.YP,FLOW.M1S)=Read_YAW(fname)	### read inlet mach 
			if FLOW.A2 == 'NaN':break	

			if FLOW.M2 == 'NaN':error=0.0
			else:error = abs(1-FLOW.M2/FLOW_Target.M2)		
			#print 'error:',error,' M2:',FLOW.M2,'M1S:',FLOW.M1S	### print convergence to screen


        if error > 0.01:	### if convergence fails return diverged case
		print 'error setting M2 -> NAN1'
		FLOW.A2='NaN'
		FLOW.M2='NaN'


	os.system('sh eval_process.txt > tmp_all')	### RUN EDP AND IPLOT, DITTO
	fname='Mach.txt'
	(FLOW.S,FLOW.CM,FLOW.M,FLOW.MTE)=Read_Mach(fname,FLOW_Target)	### READ MACH NUMBER DISTRIBUTION

	filename = 'idat.CCfoil'
	FLOW.IDAT = read_idat(filename)	### read binary file

	fname='OUTPUT'
	(FLOW.A2,FLOW.M2,FLOW.YP,FLOW.M1S)=Read_YAW(fname)	### READ YAW, YP, INLET MACH
	FLOW.A1=FLOW_Target.A1

        if error > 0.01:
		FLOW.A2='NaN'
		FLOW.M2='NaN'
	return(FLOW)



def RUN_MISES_TO_M2_INC(FLOW_Target,GEO,INC):
	### ADDITIONAL STABILITY BY STARTING FROM NEGATIVE INCIDENCE
	global PARAM_OLD  ,M1S_GLOBAL 
	### USER INPUT - DELAY FOR READ WRITE DELAYS IF THEY EXIST
	dT = 0.01	### time gaps for read write

	### insert code to solve geometry and resolve flow	#####
	FLOW=CLASS_FLOW()
	FLOW_temp=copy.deepcopy(FLOW_Target)
	if FLOW_Target.M2*FLOW_Target.CM.max()>1.1:
		FLOW_temp.M1S = FLOW_temp.M1S*(0.75/FLOW_Target.M2)**0.25
	WRITE_MISES_BLADE(GEO,FLOW_temp,'blade.CCfoil')	### write blade to MISES input
	fout='ises.CCfoil'
	WRITE_MISES_ISES(FLOW_temp,fout)		### WRITE ISES FILE SETTING UP VISCID SETUP

	time.sleep(dT)
	os.system('sh run_process.txt > tmp_all')	### RUN initial mises run
	os.system('sh run_edp_proc.txt >OUTPUT')
	error = 1.0
	#FLOW_temp=copy.deepcopy(FLOW_Target)
	it = 0
	fname='OUTPUT'
	(FLOW.A2,FLOW.M2,FLOW.YP,FLOW.M1S)=Read_YAW(fname)
	#FLOW_temp.M1S = FLOW_Target.M1S
        for inc_var in np.linspace(0,INC,int(abs(INC)/0.5)):
		error=1
		it=0
		FLOW_temp.A1=FLOW_Target.A1+inc_var
		#if inc_var > 0 :
			#FLOW_temp.M1S=FLOW_temp.M1S*0.95	### prevents overspeed from reduced incidence
		print 'incidence shift:',inc_var#,error,FLOW_temp.A1
		while error > 0.0002 and it<20:		### iterate inlet mach number
			it=it+1
			fname='OUTPUT'
				### read inlet mach 

	
			#print FLOW.A2,FLOW.M2,FLOW.A1
			if FLOW.A2 == 'NaN':break
			if FLOW.M2 == 'NaN':break
			if inc_var != INC and inc_var > 0:
				FLOW_temp.M1S = FLOW_temp.M1S*(FLOW_Target.M2*0.8/FLOW.M2)**0.2
			else:
				FLOW_temp.M1S = FLOW_temp.M1S*(FLOW_Target.M2/FLOW.M2)**0.2
			WRITE_MISES_ISES(FLOW_temp,fout)		### WRITE ISES FILE SETTING UP VISCID SETUP
			time.sleep(dT)
			os.system('sh run_ises_single.txt > tmp_all')		### run at new mach with BL iteration
			time.sleep(dT)
			
			os.system('sh run_edp_proc.txt >OUTPUT')
			time.sleep(dT)
			fname='OUTPUT'
			#FLOW_temp.M1S = FLOW_temp.M1S*(FLOW_Target.M2/FLOW.M2)**0.2
			(FLOW.A2,FLOW.M2,FLOW.YP,FLOW.M1S)=Read_YAW(fname)	### read inlet mach 
			if FLOW.A2 == 'NaN':break	

			if FLOW.M2 == 'NaN':error=0.0
			else:error = abs(1-FLOW.M2/FLOW_Target.M2)		
			print 'error:',error,' M2:',FLOW.M2,'M1S:',FLOW.M1S
			#if FLOW.M2 == 'NaN':error=0.0
			#else:error = abs(1-FLOW.M2/FLOW_Target.M2)
			#print 'error:',error,' M2:',FLOW.M2,'M1S:',FLOW.M1S

        if error > 0.005:
		FLOW.A2='NaN'
		FLOW.M2='NaN'

 	#M1S_GLOBAL =M1S_GLOBAL*0.95+0.05*FLOW_temp.M1S
	os.system('sh eval_process.txt > tmp_all')	### RUN EDP AND IPLOT, DITTO
	fname='Mach.txt'
	(FLOW.S,FLOW.CM,FLOW.M,FLOW.MTE)=Read_Mach(fname,FLOW_Target)	### READ MACH NUMBER DISTRIBUTION
	#FLOW.OD = 100.*(1.-(FLOW.MTE/FLOW.M2))
	filename = 'idat.CCfoil'
	FLOW.IDAT = read_idat(filename)
	#(GEO,FLOW) = SURFLENGTHS(GEO,FLOW)
	#FLOW.S=FLOW.S2
	#FLOW.CM=FLOW.CM2
	fname='OUTPUT'
	(FLOW.A2,FLOW.M2,FLOW.YP,FLOW.M1S)=Read_YAW(fname)	### READ YAW, YP, INLET MACH
	FLOW.A1=FLOW_Target.A1

        if error > 0.005:
		FLOW.A2='NaN'
		FLOW.M2='NaN'
	return(FLOW)



def SAVE_PARAM(PARAM,savename): ### save parameters to a text file

	f = open(savename+'_param.txt','w')
	for i in range(len(PARAM)-1):
		f.write(str(PARAM[i])+'\n')
	f.write(str(PARAM[-1]))
	#print 'PARAMETERS SAVED TO:',savename
	return()

def LOAD_PARAM(savename): ### LOAD parameters FROM a text file

	f = open(savename+'_param.txt','r')
	lines=f.readlines()
	PARAM=np.zeros((len(lines)))
	for i in range(len(PARAM)):
		PARAM[i] = float(lines[i].split()[0])

	return(PARAM)

def LOAD_PARAM2(savename): ### LOAD parameters FROM a text file - slightly different naming convention

	f = open(savename,'r')
	lines=f.readlines()
	PARAM=np.zeros((len(lines)))
	for i in range(len(PARAM)):
		PARAM[i] = float(lines[i].split()[0])

	return(PARAM)

if __name__ == "__main__":
	
	global M1S_GLOBAL	### DITTO
	global cost_record,PS_record,SS_record,DUTY_record
	
	print 'Just a file full of functions'


