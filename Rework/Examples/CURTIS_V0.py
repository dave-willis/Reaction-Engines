import numpy as np			### numpy library for general matrix operations
import os,sys,time,numpy,copy			### os,sys libraries for talking to command line on linux
import pickle,cPickle, random		### pickle libraries used for storing trained RBF networks
import matplotlib.pyplot as plt		### matplotlib used for some plotting
from matplotlib.patches import Ellipse		### lots of widgets used for gui - to be replaced soon
from matplotlib import lines
from matplotlib.widgets import Button, Slider
from scipy.interpolate import interp1d	### 1d interpolation
import scipy.interpolate as sciint	### import interpolation library from scipy
from scipy.integrate import cumtrapz	### import trapezoidal integration
import scipy.optimize as SciOpt		### import optimization library from scipy
import numpy.linalg as linalg		### import numpys linear algebra library
import subprocess			### import subprocess library
from numpy import (atleast_1d, eye, mgrid, argmin, zeros, shape, squeeze,
vectorize, asarray, sqrt, Inf, asfarray, isinf)		### load lots of default functions
cwd = (os.getcwd()).strip()
import warnings



from CURTIS_FUNCTIONS_V0 import *


###############################################################################################
### THIS IS A DEMO VERSION OF THE CAMBRIDGE UNIVSERISTY RAPID TURBOMACHINERY INVERSE SYSTEM ###
### WRITTEN BY: CHRIS CLARK - CONTACTABLE ON: cjc95@cam.ac.uk #################################
############################################################################################### 



### INITIALIZE THE GUI TRHOUGH CLASSES/STRUCTURE
### MAIN WINDOW
class Struct:
    def __init__(self, **entries): 
        self.__dict__.update(entries)

### DEFINE DRAGGABLE POINTS
class DraggablePoint:
    def __init__(self, point):
        self.point = point
        self.press = None

    def connect(self):
        'connect to all the events we need'
        self.cidpress = self.point.figure.canvas.mpl_connect(
            'button_press_event', self.on_press)
        self.cidrelease = self.point.figure.canvas.mpl_connect(
            'button_release_event', self.on_release)
        self.cidmotion = self.point.figure.canvas.mpl_connect(
            'motion_notify_event', self.on_motion)

    def on_press(self, event):
        'on button press we will see if the mouse is over us and store some data'
        if event.inaxes != self.point.axes: return

        contains, attrd = self.point.contains(event)
        if not contains: return
        #print 'event contains', self.point.center
        x0, y0 = self.point.center
        self.press = x0, y0, event.xdata, event.ydata
        #if saverst ==1:
           # bsave.label.set_text('Save')
            #bsave.label.set_fontsize(20)

    def on_motion(self, event):
        'on motion we will move the point if the mouse is over us'
        if self.press is None: return
        if event.inaxes != self.point.axes: return
        x0, y0, xpress, ypress = self.press
        dx = event.xdata - xpress
        dy = event.ydata - ypress
        #print 'x0=%f, xpress=%f, event.xdata=%f, dx=%f, x0+dx=%f'%(x0, xpress, event.xdata, dx, x0+dx)
        #dx=0.
        self.point.center=x0+dx, y0+dy
        drawcurve()
        self.point.figure.canvas.draw()

    def on_release(self, event):
        'on release we reset the press data'
        self.press = None
        drawcurve()
        self.point.figure.canvas.draw()

    def disconnect(self):
        'disconnect all the stored connection ids'
        self.point.figure.canvas.mpl_disconnect(self.cidpress)
        self.point.figure.canvas.mpl_disconnect(self.cidrelease)
        self.point.figure.canvas.mpl_disconnect(self.cidmotion)

### DEFINE HORIZONTALLY DRAGABLE POINT
class DraggableHPoint:
    def __init__(self, point):
        self.point = point
        self.press = None

    def connect(self):
        'connect to all the events we need'
        self.cidpress = self.point.figure.canvas.mpl_connect(
            'button_press_event', self.on_press)
        self.cidrelease = self.point.figure.canvas.mpl_connect(
            'button_release_event', self.on_release)
        self.cidmotion = self.point.figure.canvas.mpl_connect(
            'motion_notify_event', self.on_motion)

    def on_press(self, event):
        'on button press we will see if the mouse is over us and store some data'
        if event.inaxes != self.point.axes: return

        contains, attrd = self.point.contains(event)
        if not contains: return
        #print 'event contains', self.point.center
        x0, y0 = self.point.center
        self.press = x0, y0, event.xdata, event.ydata
        #if saverst ==1:
           ## bsave.label.set_text('Save')
           # bsave.label.set_fontsize(20)

    def on_motion(self, event):
        'on motion we will move the point if the mouse is over us'
        if self.press is None: return
        if event.inaxes != self.point.axes: return
        x0, y0, xpress, ypress = self.press
        dx = event.xdata - xpress
        dy = 0#event.ydata - ypress
        #print 'x0=%f, xpress=%f, event.xdata=%f, dx=%f, x0+dx=%f'%(x0, xpress, event.xdata, dx, x0+dx)
        #dx=0.
        self.point.center=x0+dx, y0+dy
        drawcurve()
        self.point.figure.canvas.draw()

    def on_release(self, event):
        'on release we reset the press data'
        self.press = None
        drawcurve()
        self.point.figure.canvas.draw()

    def disconnect(self):
        'disconnect all the stored connection ids'
        self.point.figure.canvas.mpl_disconnect(self.cidpress)
        self.point.figure.canvas.mpl_disconnect(self.cidrelease)
        self.point.figure.canvas.mpl_disconnect(self.cidmotion)

### DEFINE VERTICALLY DRAGABLE POINT
class DraggableVPoint:
    def __init__(self, point):
        self.point = point
        self.press = None

    def connect(self):
        'connect to all the events we need'
        self.cidpress = self.point.figure.canvas.mpl_connect(
            'button_press_event', self.on_press)
        self.cidrelease = self.point.figure.canvas.mpl_connect(
            'button_release_event', self.on_release)
        self.cidmotion = self.point.figure.canvas.mpl_connect(
            'motion_notify_event', self.on_motion)

    def on_press(self, event):
        'on button press we will see if the mouse is over us and store some data'
        if event.inaxes != self.point.axes: return

        contains, attrd = self.point.contains(event)
        if not contains: return
        #print 'event contains', self.point.center
        x0, y0 = self.point.center
        self.press = x0, y0, event.xdata, event.ydata
        #if saverst ==1:
         #   bsave.label.set_text('Save')
         #   bsave.label.set_fontsize(20)

    def on_motion(self, event):
        'on motion we will move the point if the mouse is over us'
        if self.press is None: return
        if event.inaxes != self.point.axes: return
        x0, y0, xpress, ypress = self.press
        dx = 0#event.xdata - xpress
        dy = event.ydata - ypress
        #print 'x0=%f, xpress=%f, event.xdata=%f, dx=%f, x0+dx=%f'%(x0, xpress, event.xdata, dx, x0+dx)
        #dx=0.
        self.point.center=x0+dx, y0+dy
        drawcurve()
        self.point.figure.canvas.draw()

    def on_release(self, event):
        'on release we reset the press data'
        self.press = None
        drawcurve()
        self.point.figure.canvas.draw()

    def disconnect(self):
        'disconnect all the stored connection ids'
        self.point.figure.canvas.mpl_disconnect(self.cidpress)
        self.point.figure.canvas.mpl_disconnect(self.cidrelease)
        self.point.figure.canvas.mpl_disconnect(self.cidmotion)


### DEFINE LABELED BUTTONS
class LabelledButton:
    def __init__(self):
        pass

    ### RUN BUTTON TO MAKE A CALL TO MISES
    def run(self,event):
	global PARAM,Yp
	#print 'running-mises'
	#time.sleep(5)
	A1 =sl_A1.val
	A2 = sl_A2.val
	M2=sl_mach.val
	FLOW_Target=CLASS_FLOW()
	FLOW_Target.A1 = A1
	FLOW_Target.A2 = A2
	FLOW_Target.M2 = M2	### target exit mach number (used to set inlet mach number for mises)
	FLOW_Target.Re2 = 6.0e5#1e5	### target exit Reynolds number (used to set inlet Reynolds number for mises)
	FLOW_Target.Gamma = 1.4		### gas property gamma (only change if not air)
	FLOW_Target.M1S=M2_to_M1(FLOW_Target.Gamma,FLOW_Target.M2,FLOW_Target.A1,FLOW_Target.A2)

	AVDR = 1.0
	Arat = np.cos(np.radians(A2))/np.cos(np.radians(A1))
	Aratx = np.cos(np.radians(A2))/np.cos(np.radians(0.))

	MR = 1.0
	MRx = 1.0

	for it in range(1000):
		AR = F_AR(MR,M2)
		ARx = F_AR(MRx,M2)
		if AR>Arat:MR=MR*0.99
		else:MR=MR*1.01
		if ARx>Aratx:MRx=MRx*0.99
		else:MRx=MRx*1.01
		if abs((AR-Arat)/Arat)<0.01 and abs((ARx-Aratx)/Aratx)<0.01:break

	Minrat=MR

	### calculate isentropic inlet Mach #
	M1S_GLOBAL = FLOW_Target.M1S									### set as global (possibly not used anymore)
	FLOW_Target.Re1=RE2_to_RE1(FLOW_Target.Re2,FLOW_Target.A1,FLOW_Target.A2)	
    	GEO = PARAM_TO_GEO(PARAM[:42])
        FLOW=RUN_MISES_TO_M2(FLOW_Target,GEO)
	Yp=FLOW.YP

	#Get turning control points
	npts=len(points)
	xcen=np.zeros(npts)
	ycen=np.zeros(npts)
	for i, point in enumerate(points):
		xcen[i]=point.center[0]
		ycen[i]=point.center[1]

   
	Mless = ycen[0]*Minrat*2.
	Mps = ycen[-1]*MRx
	#Mps = ycen[-1]*np.cos(np.radians(A2))
	#Mless= ycen[0]*np.cos(np.radians(A2))/np.cos(np.radians(A1))
	lsss = xcen[1]	### set start of SS ramp
	lsps = xcen[-2]	### set start of PS hold
	leps = xcen[-1]	### set end of PS hold
	gte = 2.0	### set diffusion grad. at TE
	lpeak = xcen[2]	### set Mpeak location
	Mpeak = ycen[2]	### set Mpeak
	(X,Y)=Analytic_CM(Mpeak,Mless,Mps,lpeak,lsss,lsps,leps,gte)
	#Update Plot
	line_TARGET.set_data([X,Y])#shifted for axes purposes
        (FLOW_Target.S,FLOW_Target.CM)=(X,Y)#Analytic_CM(Mpeak,Mless,Mps,lpeak,lsss,lsps,leps,gte)
	#plt.figure()
	#plt.plot(FLOW_Target.S,FLOW_Target.CM,'-xr')
	#plt.show()
	RSME=RMS_COST(FLOW,FLOW_Target)
	#print Yp,RSME
	#print FLOW.S,FLOW.CM
	#line_TARGET.set_data([FLOW_Target.S,FLOW_Target.CM])
	line_TRUE.set_data([FLOW.S,FLOW.CM])

	Yp_calc.set_text('Yp:'+str(round(FLOW.YP*100.,4))+'%')
	alpha_calc.set_text('A2:'+str(round(FLOW.A2,3)))
	P2C_calc.set_text('RMSE:'+str(round(RSME,3)))
	#ax1.set_xlabel('l/L'+str(FLOW.Yp),size=20,rotation='horizontal')
	#drawcurve()
	#print 'running-mises'

    ### REDUNDENT BUTTON PREVIOUSLY USED FOR SWITCHING TO PITCH-TO-CHORD MODE
    def mode(self,event):
	global mode,functions,A1int,A2int,M2int,MPint,LPint,PARAMint
	global A1intN,A2intN,M2intN,MPintN,LPintN
	
	if mode == True:
		print 'switching to controlling P2C'
		mode = False
		bmode.label.set_text('Mode=P2C')

		sl_AVDR.set_val( PARAM[3])
		'''	
		A2s= [64,66,68,70,72,74]#np.linspace(-55,5,6):
		A1s= [-20,-30,-40,-50][::-1]
		M2s=[0.4,0.5,0.6,0.7]
		MPs = [1.2,1.3,1.4]
		LPs = [0.5,0.6]

		A1int = np.zeros((4,6,4,3,2))
		A2int = np.zeros((4,6,4,3,2))
		M2int = np.zeros((4,6,4,3,2))
		MPint = np.zeros((4,6,4,3,2))
		LPint = np.zeros((4,6,4,3,2))
		PARAMint = np.zeros((42,4,6,4,3,2))
		for i in range(4):
			for j in range(6):
			   for k in range(4):
			     for l in range(3):
			      for m in range(2):
				A1 = A1s[i]
				A2 = A2s[j]
				M2 = M2s[k]
				Mp = MPs[l]
				Lp = LPs[m]
				if Lp==0.5:
					savename =  'DB_'+str(A2)+'_'+str(A1)+'_'+str(M2)+'_'+str(Mp)
				else:
					savename = 'DB_'+str(A2)+'_'+str(A1)+'_'+str(M2)+'_'+str(Mp)+'_'+str(Lp)+'_'+str(0.8)+'_'+str(1.2)
				A1int[i,j,k,l,m]=A1
				A2int[i,j,k,l,m]=A2
				M2int[i,j,k,l,m]=M2
				MPint[i,j,k,l,m]=Mp
				LPint[i,j,k,l,m]=Lp
				#PARAMint[:,i,j,k,l,m] = LOAD_PARAM(savename)
		'''
		MPint = PARAMint[3,:,:,:,:,:]

		A1intN = scalarise(A1int.flatten())
		A2intN = scalarise(A2int.flatten())
		M2intN = scalarise(M2int.flatten())
		MPintN = scalarise(MPint.flatten())*2.
		LPintN = scalarise(LPint.flatten())
		functions = []
		for i in range(43):
				f=sciint.Rbf(A1intN,A2intN,M2intN,MPintN,LPintN,PARAMint[i,:,:,::].flatten())#,epsilon=2.)
				functions.append(f)
		return()
	else:
		print 'switching to controlling Mpeak'
		mode = True
		bmode.label.set_text('Mode=Mpeak')

		A2s= [64,66,68,70,72,74]#np.linspace(-55,5,6):
		A1s= [-20,-30,-40,-50][::-1]
		M2s=[0.4,0.5,0.6,0.7]
		MPs = [1.2,1.3,1.4]
		LPs = [0.5,0.6]

		#A1int = np.zeros((4,6,4,3,2))
		#A2int = np.zeros((4,6,4,3,2))
		#M2int = np.zeros((4,6,4,3,2))
		MPint = np.zeros((4,6,4,3,2))
		#LPint = np.zeros((4,6,4,3,2))
		#PARAMint = np.zeros((43,4,6,4,3,2))
		for i in range(4):
			for j in range(6):
			   for k in range(4):
			     for l in range(3):
			      for m in range(2):
				A1 = A1s[i]
				A2 = A2s[j]
				M2 = M2s[k]
				Mp = MPs[l]
				Lp = LPs[m]
				if Lp==0.5:
					savename =  'DB_'+str(A2)+'_'+str(A1)+'_'+str(M2)+'_'+str(Mp)
				else:
					savename = 'DB_'+str(A2)+'_'+str(A1)+'_'+str(M2)+'_'+str(Mp)+'_'+str(Lp)+'_'+str(0.8)+'_'+str(1.2)
				#A1int[i,j,k,l,m]=A1
				#A2int[i,j,k,l,m]=A2
				#M2int[i,j,k,l,m]=M2
				MPint[i,j,k,l,m]=Mp
				#LPint[i,j,k,l,m]=Lp
				#PARAMint[:42,i,j,k,l,m] = LOAD_PARAM(savename)

		#MPint = PARAMint[3,:,:,:,:,:]

		A1intN = scalarise(A1int.flatten())
		A2intN = scalarise(A2int.flatten())
		M2intN = scalarise(M2int.flatten())
		MPintN = scalarise(MPint.flatten())
		LPintN = scalarise(LPint.flatten())
		functions = []
		for i in range(43):
				f=sciint.Rbf(A1intN,A2intN,M2intN,MPintN,LPintN,PARAMint[i,:,:,::].flatten())#,epsilon=2.)
				functions.append(f)
	
### DEFINE SLIDERS
class SliderBox:
    def __init__(self):
        pass
    ### INLET ANGLES
    def A1(self,val):
        global A1
        A1 =sl_A1.val
        drawcurve()
    ### EXIT ANGLES
    def A2(self,val):
        global A2
        A2 = sl_A2.val
        drawcurve()
    ### EXIT MACH NUMBER
    def mach(self,val):
        global mach
        mach = sl_mach.val
        drawcurve()
    ### AVDR - NOT INCLUDED IN VERSION 0
    def AVDR(self,val):
        global AVDR
        AVDR = sl_AVDR.val
        drawcurve()
    ### GTE - TRAILING EDGE PRESSURE GRADIENT - NOT INCLUDED IN VERSION 0
    def GTE(self,val):
        global GTE
        GTE = sl_GTE.GTE
        drawcurve()

### AREA RATIO CALCULATION - NOTE HARD CODED GAMMA IN VERSION 0
def F_AR(MR,M2):
	AR = MR*((1.+0.5*0.4*(MR*M2)**2.0)/(1.+0.5*0.4*(M2)**2.0))**(-0.5*2.4/0.4)
	return(AR)
### NORMALIZATION OF ARRAYS
def scalarise(array):
	array = (array-array.min())/(array.max()-array.min())
	return(array)

### COMMAND TO REDRAW MAIN GUI
def drawcurve():
	#print 'drawing curve'
	global PARAM,Yp,mode,functions,A1int,A2int,M2int,MPint,LPint,PARAMint
	global A1intN,A2intN,M2intN,MPintN,LPintN

	### LOAD VALUES BASED ON CURRENT SLIDER AND POINT POSITIONS
	A1 =sl_A1.val
	A2 = sl_A2.val
	GTE = sl_GTE.val
	M2 = sl_mach.val
	AVDR = sl_AVDR.val
	Arat = np.cos(np.radians(A2))/np.cos(np.radians(A1))
	Aratx = np.cos(np.radians(A2))/np.cos(np.radians(0.))

	MR = 1.0
	MRx = 1.0

	for it in range(1000):
		AR = F_AR(MR,M2)
		ARx = F_AR(MRx,M2)
		if AR>Arat:MR=MR*0.99
		else:MR=MR*1.01
		if ARx>Aratx:MRx=MRx*0.99
		else:MRx=MRx*1.01
		if abs((AR-Arat)/Arat)<0.01 and abs((ARx-Aratx)/Aratx)<0.01:break

	Minrat=MR


	#Get turning control points
	npts=len(points)
	xcen=np.zeros(npts)
	ycen=np.zeros(npts)
	for i, point in enumerate(points):
		xcen[i]=point.center[0]
		ycen[i]=point.center[1]

   
	Mless = ycen[0]*Minrat*2.	### SCALING FOR LE SS BASED ON INLET MACH
	Mps = ycen[-1]*MRx		### SCALING FOR PS BASED ON AXIAL MACH

	lsss = xcen[1]	### set start of SS ramp
	lsps = xcen[-2]	### set start of PS hold
	leps = xcen[-1]	### set end of PS hold
	gte = GTE	### set diffusion grad. at TE
	lpeak = xcen[2]	### set Mpeak location
	Mpeak = ycen[2]	### set Mpeak

	### CALCULATE ANALYTIC LOADING DISTRIBUTION
	(X,Y)=Analytic_CM(Mpeak,Mless,Mps,lpeak,lsss,lsps,leps,gte)
	#Update Plot
	line_TARGET.set_data([X,Y]) #shifted for axes purposes
	
	### SET INLET MACH NUMBER LINE
	XLE,YLE = [0,1],[Minrat,Minrat]
	line_M1.set_data([XLE,YLE])


	
	### PERFORM NORMALISATION OF VARIABLES FOR INTERPOLATION - V0 HAS LIMITS IN HERE TO AVOID EXTRAPOLATION
	A1N = (A1-A1int.min())/(A1int.max()-A1int.min())
	A2N = (A2-A2int.min())/(A2int.max()-A2int.min())
	M2N=(M2-M2int.min())/(M2int.max()-M2int.min())
	MPSN=min(1,max(0,(ycen[-1]-MPSint.min())/(MPSint.max()-MPSint.min())))
	MLESSN=min(1,max(0,(ycen[0]*2.-MLESSint.min())/(MLESSint.max()-MLESSint.min())))
	if mode == True:
		MPN=min(1,max(0,(Mpeak-MPint.min())/(MPint.max()-MPint.min())))
	else:
		MPN=(AVDR-MPint.min())/(MPint.max()-MPint.min())*2.
	LPN=min(1,max(0,(lpeak-LPint.min())/(LPint.max()-LPint.min())))

	### CALCULATE MINIMUM DISTANCE FROM DATABASE ENTRY
	distance = ((((A1intN-A1N)**2.)+((A2intN-A2N)**2.)+((M2intN-M2N)**2.)+((MPintN-MPN)**2.)+((LPintN-LPN)**2.)+((MLESSintN-MLESSN)**2.)+((MPSintN-MPSN)**2.))**0.5).min()

	### APPLY YAW ANGLE CORRECTION - ASK FOR MORE OR LESS BASED ON PREDICTED ERROR
	f2=functions[-1]
	deltaA2 = f2(A1N,A2N,M2N,MPN,LPN,MPSN,MLESSN)
	A2now = A2+deltaA2
	A2N = (A2now-A2int.min())/(A2int.max()-A2int.min())

	### CALL RBF FUNCTIONS TO INTERPOLATE BLADE SHAPE ###
	PARAM=np.zeros((43))
	PARAM2=np.zeros((43))
	for i in range(43):
		f2 = functions[i]
		PARAM[i]=f2(A1N,A2N,M2N,MPN,LPN,MPSN,MLESSN)

	### CONVERT BLADE PARAMETERS INTO A GEOMETRY
    	GEO = PARAM_TO_GEO(PARAM[:42])

	### UPDATE PLOT OF BLADE SECTION
	line_profile1.set_data([GEO.X,GEO.Z-GEO.Z.max()+0.2])
	line_profile2.set_data([GEO.X,GEO.Z+PARAM[3]-GEO.Z.max()+0.2])


	### UPDATE TEXT INDICATORS
	Yp_pred.set_text('Yp:'+str(100.*round(PARAM[42],4))+'%')
	alpha_pred.set_text('A2:'+str(round(A2,3)))
	P2C_pred.set_text('P2C:'+str(round(PARAM[3],3)))
	conf_pred.set_text('Confidence:'+str(100.*(1.-round(distance,3)))+'%')	### CONFIDENCE IS BASED ON DISTANCE FROM NEAREST ENTRY


### CALCULATE RMS ERROR OF LOADING STYLE
def RMS_COST(FLOW,FLOW_Target):	### MODULE: FLOW-> COST FUNCTION - CURRENTLY USED VERSION OF COSTINGS

    ### SAMPLE DISTRIBUTION IN ALL CASES
    Nsamples=200								### NUMBER OF SAMPLES FOR COSTS
    (XS_Target,YS_Target) = F_sampleSURF(FLOW_Target.S,FLOW_Target.CM,Nsamples)	### SAMPLE TARGET DISTRIBUTION
    (XS,YS) = F_sampleSURF(FLOW.S,FLOW.CM,Nsamples)				### SAMPLE ACTUAL DISTRIBUTION
    
    if 1 == 1:								### INVERSE DESIGN MODE		
	    
	    ### USER OPTIONS
	    FSKIP_PS = 0.01					### IGNORE INITIAL FRACTION OF SURFACE -PS
	    FSKIP_SS = 0.01					### IGNORE INITIAL FRACTION OF SURFACE -SS    
	  

	    COSTPS = 0.0	### INITIALISE PS COST = 0
	    COSTSS=0.0		### INITIALISE SS COST = 0

	    ########### PRESSURE SURFACE #############
	    for i in range(len(YS)/2):

		dist = ((YS[i]-YS_Target[i])**2) 			### LEAST SQUARES  COST
		COSTPS=COSTPS+(dist)	
	    ########### SUCTION SURFACE ##############
	    for i in range(len(YS)/2,len(YS)):
						
		dist = ((YS[i]-YS_Target[i])**2)			### LEAST SQUARES COST
		COSTSS=COSTSS+(dist)	

	    COSTPS=COSTPS/Nsamples		### PERFORM AVERAGE SO CHANGING NSAMPLES DOESNT CHANGE TOLERENCES ETC
	    COSTSS=COSTSS/Nsamples
	    COST=(COSTPS+COSTSS)**0.5		### SUM SS AND PS COSTS 

    return(COST)	



#################################################
### Turbine Interpolation and Database system ###
### Written by Chris Clark 2018 #################
#################################################


##########################
####Start Main Routine ###
##########################


### intialise a few global parameters
global PARAM,Yp,mode,functions,A1int,A2int,M2int,MPint,LPint,PARAMint
global A1intN,A2intN,M2intN,MPintN,LPintN,PARAMintN

warnings.filterwarnings("ignore", category=DeprecationWarning) 
PARAM = []
Yp=0.
mode=True
### default loading distribution
Mless = 0.8
Mps = 0.8
lsss = 0.1	### set start of SS ramp
lsps = 0.1	### set start of PS hold
leps = 0.3	### set end of PS hold
gte = 2.0	### set diffusion grad. at TE
lpeak = 0.5	### set Mpeak location
Mpeak = 1.2	### set Mpeak


###########################
### load design database###
###########################

### DATABASE VERSION 0 - AS DESCRIBED IN CLARK 2019
A2s= [65.0,67.5,70.0,72.5]
A1s= [-20,-35,-50][::-1]
M2s=[0.4,0.55,0.7]
MPs = [1.2,1.3,1.4]
LPs = [0.5,0.55,0.6]
MPSs = [0.8,0.9,1.0]
MLESSs = [1.2,1.5,1.8]

A1int = np.zeros((4,3,3,3,3,3,3))
A2int = np.zeros((4,3,3,3,3,3,3))
M2int = np.zeros((4,3,3,3,3,3,3))
MPint = np.zeros((4,3,3,3,3,3,3))
LPint = np.zeros((4,3,3,3,3,3,3))
MPSint = np.zeros((4,3,3,3,3,3,3))
MLESSint = np.zeros((4,3,3,3,3,3,3))
PARAMint = np.zeros((44,4,3,3,3,3,3,3))


### LOAD ADDITIONAL DATA - YP YAW CORRECTION
f3 = open('Database/extra.txt','r')
YPs=[]
YAWTs=[]
linesa = f3.readlines()
for i in range(len(linesa)):
	YPs.append(float(linesa[i].split()[1]))
	YAWTs.append(float(linesa[i].split()[0]))

f3.close()

print 'STARTING CAMBRIDGE UNIVERSITY RAPID TURBINE INVERSE SYSTEM'
print 'Version: 1.0 - July 2019'
print 'written by Chris Clark : cjc95@cam.ac.uk'
print 'booting...'
time.sleep(1.0) ### temporary pause so people see the above title
icount=-1
#print 'loading raw database...'
for j in range(3):
   	for i in range(4):
	   for k in range(3):
	     for l in range(3):
	      for m in range(3):
	      	for n in range(3):
	      	  for o in range(3):
			icount=icount+1
			#icount=icount+1 
			#print icount

			A2 = A2s[i]
			A1 = A1s[j]
		
			M2 = M2s[k]
			Mp = MPs[l]
			Lp = LPs[m]
			Mps_rat = MPSs[n]
			Mless_rat = MLESSs[o]
			savename ='Database/CURTIS_'+str(A2)+'_'+str(A1)+'_'+str(M2)+'_'+str(Mp)+'_'+str(Lp)+'_'+str(Mps_rat)+'_'+str(Mless_rat)
			A1int[i,j,k,l,m,n,o]=A1
			
			A2int[i,j,k,l,m,n,o]=A2#YAWTs[icount]#A2

			M2int[i,j,k,l,m,n,o]=M2
			MPint[i,j,k,l,m,n,o]=Mp
			LPint[i,j,k,l,m,n,o]=Lp
			MPSint[i,j,k,l,m,n,o]=Mps_rat
			MLESSint[i,j,k,l,m,n,o]=Mless_rat

			PARAMint[:-2,i,j,k,l,m,n,o] = LOAD_PARAM(savename)
			PARAMint[-2,i,j,k,l,m,n,o]=YPs[icount]
			PARAMint[-1,i,j,k,l,m,n,o]=A2-YAWTs[icount]

		
A1intN = scalarise(A1int.flatten())
A2intN = scalarise(A2int.flatten())
M2intN = scalarise(M2int.flatten())
MPintN = scalarise(MPint.flatten())
LPintN = scalarise(LPint.flatten())
MLESSintN = scalarise(MLESSint.flatten())
MPSintN = scalarise(MPSint.flatten())
functions = []

############################################################
### USER INPUT - LOAD = 0, AVOIDS REFITTING RBF NETWORKS ###
############################################################
load=1

### set this to 0 after the first run that generated the RBF networks

if load==1:
	### USE NORMALISED DATABASES TO FIT AND PICKLE NEW RBF NETWORKS
	print 'fitting new RBFs...'
	for i in range(44):
			print 'RBF:',i+1,'/43'
			### GENERATE NEW RBF NETWORK AND SAVE
			f=sciint.Rbf(A1intN,A2intN,M2intN,MPintN,LPintN,MPSintN,MLESSintN,PARAMint[i,:,:,:,:,:,:,:].flatten())#,epsilon=2.)
			functions.append(f)
			functionname = 'Pickles/RBF_F2'+str(i)+'.pkl'
	
			
			RBFfile = open(functionname,'wb')
			RBFpickler = cPickle.Pickler(RBFfile,protocol=2)

			# RBF can't be pickled directly, so save everything required for reconstruction
			RBFdict = {}            
			for key in f.__dict__.keys():
			    if key != '_function' and key!= 'norm':
				RBFdict[key] = f.__getattribute__(key)   

			RBFpickler.dump(RBFdict)
			RBFfile.close()
			### STORE RBF AS PICKLED OBJECT

else:
	### LOAD PICKLED RBF NETWORKS - MUCH FASTER
	print 'loading RBFs...'
	for i in range(44):
			#print 'RBF:',i+1,'/43'
			functionname = 'Pickles/RBF_F2'+str(i)+'.pkl'
			
			#This object's data parts are then replaced with the saved data:

			RBFfile = open(functionname,'rb')
			RBFunpickler = cPickle.Unpickler(RBFfile)
			RBFdict = RBFunpickler.load()
			RBFfile.close()
			f = sciint.Rbf(np.array([1,2,3]), np.array([10,20,30]), np.array([1,2,3]), function = RBFdict['function'] )
			## replace rbfi's contents with what was saved ##
			for key,value in RBFdict.iteritems():
			    f.__setattr__(key, value)
			#output = open(functionname, 'rb')
			#f=pickle.load( output)
			functions.append(f)


# RADIAL BASIS FUNCTIONS ALL LOADED NOW

########################
##### initialise GUI ###
########################

#Make Figure --- set title etc
figsize=(24,12)
fig = plt.figure(figsize=figsize,facecolor='lightgrey',edgecolor='grey')
fig.canvas.set_window_title('Cambridge University Rapid Turbine Inverse System')
#Two Subplots
##subplot 1 -design style
ax1 = plt.axes([0.05,0.3,0.4,0.65])
ax1.set_title('Style',size = 20)
ax1.set_xlim([0,1])
ax1.set_ylim([0.0,1.5])
ax1.set_xlabel('Surface Fraction',size=15,rotation='horizontal')
ax1.set_ylabel('Mach/Mach @ TE',size=15,rotation='vertical')

## subplot 2 - blade profiles
ax2 = plt.axes([0.55,0.3,0.4,0.65])
ax2.set_title('Blade Profile',size=20)
ax2.set_xlim([-0.7, 1.7])
ax2.set_ylim([-1.2, 1.8])
ax2.set_xlabel('Meridional distance / Meridional Chord',size=15,rotation='horizontal')
ax2.set_ylabel('Tangential distance / Meridional Chord',size=15,rotation='vertical')


#### ADD LINES TO SUBPLOT 1
line_M1=lines.Line2D([0,1],[0,1],color='r',label = 'M1/M2',linewidth=2.)
line_TARGET=lines.Line2D([0,1],[0,1],color='k',label = 'Target',linewidth=2.)
line_TRUE=lines.Line2D([0,1],[1,1],color='b',label = 'True',linewidth=2.)

ax1.add_line(line_M1)
ax1.add_line(line_TARGET)
ax1.add_line(line_TRUE)
ax1.legend()

### ADD LINES TO SUBPLOT 2
line_profile1=lines.Line2D([0,1],[0,1],color='k',linewidth=2.)
line_profile2=lines.Line2D([0,1],[0,1],color='k',linewidth=2.)
ax2.add_line(line_profile1)
ax2.add_line(line_profile2)


### INITIALISE CONTROL POINTS

ncpts = 5
points=[]
xcen = [0.01,lsss,lpeak,lsps,leps]
ycen = [Mless,1.49,Mpeak,0.01,Mps]

for cpt in range(ncpts):
    circle=Ellipse((xcen[cpt],ycen[cpt]),0.025,0.03,clip_on=False,color='k')
    points.append(ax1.add_patch(circle))

#Make control points draggable - H - HORIZONTAL ONLY, V VERTICAL ONLY, OR BOTH
dps = []
dp = DraggableVPoint(points[0])
dp.connect()
dps.append(dp)
dp = DraggableHPoint(points[1])
dp.connect()
dps.append(dp)
dp = DraggablePoint(points[2])
dp.connect()
dps.append(dp)
dp = DraggableHPoint(points[3])
dp.connect()
dps.append(dp)
dp = DraggablePoint(points[4])
dp.connect()
dps.append(dp)

### Add Buttons ####
callback = LabelledButton()
axrun = plt.axes([0.85,0.15,0.1,0.075])
brun = Button(axrun, 'Validate')
brun.label.set_fontsize(20)
brun.on_clicked(callback.run)


callback = LabelledButton()
axmode= plt.axes([0.85,0.05,0.1,0.075])
bmode = Button(axmode, 'Mode=Mpeak')
bmode.label.set_fontsize(15)
bmode.on_clicked(callback.mode)

### Add DUTY Sliders
fig.text(0.31,0.22,'Duty Definition',size=20,rotation='horizontal',horizontalalignment='right',verticalalignment='center')
axA1 = plt.axes([0.05,0.16,0.4,0.04])
sl_A1 = Slider(axA1, 'Alpha1',float(-50),float(-20),valinit=-30.,color='k')
sl_A1.on_changed(SliderBox().A1)

axA2 = plt.axes([0.05,0.11,0.4,0.04])
sl_A2 = Slider(axA2, 'Alpha2',65,72.5,valinit=70,color='k')
sl_A2.on_changed(SliderBox().A2)

axmach = plt.axes([0.05,0.06,0.4,0.04])
sl_mach = Slider(axmach, 'Mach2',float(0.4),float(0.7),valinit=0.5,color='k')
sl_mach.on_changed(SliderBox().mach)

axAVDR = plt.axes([0.05,0.01,0.4,0.04])
sl_AVDR = Slider(axAVDR, 'P2C',float(1.0),float(1.4),valinit=1.0,color='k')
sl_AVDR.on_changed(SliderBox().AVDR)


### add design style slider
axGTE = plt.axes([0.75,-0.1,0.2,0.04])
sl_GTE = Slider(axGTE, 'GTE',float(0.0),float(3.0),valinit=2.0,color='k')
sl_GTE.on_changed(SliderBox().GTE)

### add permenant text
xleft = 0.55
xright=0.73
fig.text(0.625,0.23,'Performance:',size=20,rotation='horizontal',horizontalalignment='left',verticalalignment='center')
fig.text(xleft,0.2,'Predicted:',size=20,rotation='horizontal',horizontalalignment='left',verticalalignment='center')
Yp_pred =fig.text(xleft,0.15,'Yp:'+str(0),size=20,rotation='horizontal',horizontalalignment='left',verticalalignment='center')
alpha_pred =fig.text(xleft,0.1,'A2:'+str(0),size=20,rotation='horizontal',horizontalalignment='left',verticalalignment='center')
P2C_pred =fig.text(xleft,0.06,'P2C:'+str(0),size=20,rotation='horizontal',horizontalalignment='left',verticalalignment='center')
conf_pred =fig.text(xleft,0.02,'P2C:'+str(0),size=20,rotation='horizontal',horizontalalignment='left',verticalalignment='center')
fig.text(xright,0.2,'Calculated:',size=20,rotation='horizontal',horizontalalignment='left',verticalalignment='center')
Yp_calc =fig.text(xright,0.15,'Yp:'+str(0),size=20,rotation='horizontal',horizontalalignment='left',verticalalignment='center')
alpha_calc =fig.text(xright,0.1,'A2:'+str(0),size=20,rotation='horizontal',horizontalalignment='left',verticalalignment='center')
P2C_calc =fig.text(xright,0.05,'P2C:'+str(0),size=20,rotation='horizontal',horizontalalignment='left',verticalalignment='center')


#Fit splines and draw plot
print 'READY.'
drawcurve()
mng = plt.get_current_fig_manager()
plt.show()

