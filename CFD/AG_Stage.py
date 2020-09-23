import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from ts import ts_tstream_grid, ts_tstream_type,ts_tstream_load_balance
from ts import ts_tstream_reader
import scipy
import ts_tstream_default
import copy
#import PROF_FIT as PF
import math,time
import scipy.optimize as SciOpt
import scipy.interpolate as SciInt
from ts import ts_autogrid_reader 

def set_default(g,):
    # set default application variables
    for name in ts_tstream_default.av:
        val = ts_tstream_default.av[name]
        if type(val) == type(1):
            g.set_av(name, ts_tstream_type.int, val)
        else:
            g.set_av(name, ts_tstream_type.float, val)
            
    # set default block variables
    for bid in g.get_block_ids():
        for name in ts_tstream_default.bv:
	    #print name
            val = ts_tstream_default.bv[name]
            if type(val) == type(1):
                g.set_bv(name, ts_tstream_type.int, bid, val)
            else:
                g.set_bv(name, ts_tstream_type.float, bid, val)


    trans_dyn_vis=(0.84*(2.0+4.6))*1.8E-5
    for bid in g.get_block_ids():
        		b=g.get_block(bid)
        		g.set_bp("trans_dyn_vis",ts_tstream_type.float,bid,np.zeros([b.nk,b.nj,b.ni],np.float32)+1.8e-5)

bcs_path = 'TEST.bcs'
p3d_path = 'TEST.g'###

### gas properties to run
Cpgas=5187
Gamma=1.6625
Rgas = 2067.
### inlet conditions
P0inlet = 14500000.
T0inlet = 950.
ALPHAinlet = 0.0
PITCHinlet = 0.0
#exit condition
pstat = 14300000.
### initial guess 
Vxguess = 50.0
roguess = pstat/(Rgas*T0inlet)

agr= ts_autogrid_reader.AutogridReader()
g = agr.read(bcs_path, p3d_path, sliding=False)

set_default(g)

for bid in g.get_block_ids():
	for pid in g.get_patch_ids(bid):
		p=g.get_patch(bid,pid)
		print p.pid,p.kind
		### set inlet conditions
		if p.kind == 0:
			bid0=bid
			p0=g.get_patch(bid,pid)
			

			pstag = np.zeros((p0.ken-p0.kst, p0.jen-p0.jst, p0.ien-p0.ist), np.float32)
			tstag = np.zeros((p0.ken-p0.kst, p0.jen-p0.jst, p0.ien-p0.ist), np.float32)
			pitch = np.zeros((p0.ken-p0.kst, p0.jen-p0.jst, p0.ien-p0.ist), np.float32)
			yaw = np.zeros((p0.ken-p0.kst, p0.jen-p0.jst, p0.ien-p0.ist), np.float32)
			radius_in=(g.get_bp('r',0)[0,:,0])
			spanfrac_in = (radius_in-radius_in.min())/(radius_in.max()-radius_in.min())

			#print spanfrac_in
			for j in range(len(pstag[:,0])):
			    for k in range(len(pstag[0,:])):
				pstag[j,k,:]=P0inlet#np.interp(spanfrac_in[k],Rinlet,P0)
				tstag[j,k,:]=T0inlet#p.interp(spanfrac_in[k],Rinlet,T0)
				yaw[j,k,:]=ALPHAinlet#-np.interp(spanfrac_in[k],Rinlet,Alpha)
				pitch[j,k,:]=PITCHinlet#np.interp(spanfrac_in[k],Rinlet,Beta)

			g.set_pv("rfin",ts_tstream_type.float,p0.bid,p0.pid,0.5)
			g.set_pv("sfinlet",ts_tstream_type.float,p0.bid,p0.pid,0.0)
		     	g.set_pp("pstag", ts_tstream_type.float, p0.bid, p0.pid, pstag)
		    	g.set_pp("tstag", ts_tstream_type.float, p0.bid, p0.pid, tstag)    	
		    	g.set_pp("pitch", ts_tstream_type.float, p0.bid, p0.pid, pitch)   
		    	g.set_pp("yaw", ts_tstream_type.float, p0.bid, p0.pid, yaw)
			p0=(g.get_pp(bid,pid,'pstag'))
			t0=(g.get_pp(bid,pid,'tstag'))


		if p.kind == 1:
			print 'adding 2d exit plane'
			### add inlet patch

			bid0=bid
			b0 = g.get_block(bid0)
			p0 = g.get_patch(bid,pid)
			p0.kind=1
			g.set_pv("pout",ts_tstream_type.float,bid0,p0.pid,pstat)
			#g.set_pv("pout",ts_tstream_type.float,bid0,p0.pid,147779.2)
			g.set_pv("ipout",ts_tstream_type.int,bid0,p0.pid,+3)
			g.set_pv("throttle_type",ts_tstream_type.int,bid0,p0.pid,1)
			g.set_pv("throttle_target",ts_tstream_type.float,bid0,p0.pid,16.0)
			g.set_pv("throttle_k0",ts_tstream_type.float,bid0,p0.pid,20.)
			g.set_pv("throttle_k1",ts_tstream_type.float,bid0,p0.pid,100.)
			g.set_pv("throttle_k2",ts_tstream_type.float,bid0,p0.pid,200.0)
			

for bid in g.get_block_ids():
	b=g.get_block(bid)
	g.set_bv('pstatin',ts_tstream_type.float,b.bid,pstat)
	g.set_bv('pstatout',ts_tstream_type.float,b.bid,pstat)
	g.set_bv('tstagin',ts_tstream_type.float,b.bid,np.average(T0inlet))
	g.set_bv('tstagout',ts_tstream_type.float,b.bid,np.average(T0inlet))
	g.set_bv('vgridin',ts_tstream_type.float,b.bid,Vxguess)
	g.set_bv('vgridout',ts_tstream_type.float,b.bid,Vxguess)

	#for bp in ["ro","rovx","rovr","rorvt","roe","trans_dyn_vis"]:
	ro = copy.deepcopy(g.get_bp('x',bid))
	ro[:,:,:]=roguess
	g.set_bp("ro", ts_tstream_type.float, bid, ro)
	rovx = copy.deepcopy(g.get_bp('x',bid))
	rovx[:,:,:]=roguess*Vxguess
	g.set_bp("rovx", ts_tstream_type.float, bid, rovx)
	rovr = copy.deepcopy(g.get_bp('x',bid))
	rovr[:,:,:]=0.0
	g.set_bp("rovr", ts_tstream_type.float, bid, rovr)
	rorvt = copy.deepcopy(g.get_bp('x',bid))
	rorvt[:,:,:]=0.0
	g.set_bp("rorvt", ts_tstream_type.float, bid, rorvt)
	roe = copy.deepcopy(g.get_bp('x',bid))
	roe[:,:,:]=roguess*Cpgas*T0inlet
	g.set_bp("roe", ts_tstream_type.float, bid, roe)

RPM = 6782.
for bid in g.get_block_ids()[:7]:
	g.set_bv('rpm',ts_tstream_type.float,bid,0.0)	### BLOCK ALWAYS ROTATES WITH ATTACHED BLADE
	g.set_bv('rpmi1',ts_tstream_type.float,bid,0.0)		### AXIAL INLET FACE DOESNT ROTATE
	g.set_bv('rpmi2',ts_tstream_type.float,bid,0.0)
	g.set_bv('rpmj1',ts_tstream_type.float,bid,RPM)	### set it so it rotates with subsequent row
	g.set_bv('rpmj2',ts_tstream_type.float,bid,0.0) ### ste it so it rotates with previous rotor
	g.set_bv('rpmk1',ts_tstream_type.float,bid,0.0)	### K SURFACES SPIN WITH BLOCK
	g.set_bv('rpmk2',ts_tstream_type.float,bid,0.0)
for bid in g.get_block_ids()[7:23]:
	g.set_bv('rpm',ts_tstream_type.float,bid,0.0)	### BLOCK ALWAYS ROTATES WITH ATTACHED BLADE
	g.set_bv('rpmi1',ts_tstream_type.float,bid,0.0)		### AXIAL INLET FACE DOESNT ROTATE
	g.set_bv('rpmi2',ts_tstream_type.float,bid,0.0)
	g.set_bv('rpmj1',ts_tstream_type.float,bid,0.0)	### set it so it rotates with subsequent row
	g.set_bv('rpmj2',ts_tstream_type.float,bid,0.0) ### ste it so it rotates with previous rotor
	g.set_bv('rpmk1',ts_tstream_type.float,bid,0.0)	### K SURFACES SPIN WITH BLOCK
	g.set_bv('rpmk2',ts_tstream_type.float,bid,0.0)
for bid in g.get_block_ids()[23:30]:
	g.set_bv('rpm',ts_tstream_type.float,bid,RPM)	### BLOCK ALWAYS ROTATES WITH ATTACHED BLADE
	g.set_bv('rpmi1',ts_tstream_type.float,bid,RPM)		### AXIAL INLET FACE DOESNT ROTATE
	g.set_bv('rpmi2',ts_tstream_type.float,bid,RPM)
	g.set_bv('rpmj1',ts_tstream_type.float,bid,0.0)	### set it so it rotates with subsequent row
	g.set_bv('rpmj2',ts_tstream_type.float,bid,RPM) ### ste it so it rotates with previous rotor
	g.set_bv('rpmk1',ts_tstream_type.float,bid,RPM)	### K SURFACES SPIN WITH BLOCK
	g.set_bv('rpmk2',ts_tstream_type.float,bid,RPM)
for bid in g.get_block_ids()[30:46]:
	g.set_bv('rpm',ts_tstream_type.float,bid,0.0)	### BLOCK ALWAYS ROTATES WITH ATTACHED BLADE
	g.set_bv('rpmi1',ts_tstream_type.float,bid,0.0)		### AXIAL INLET FACE DOESNT ROTATE
	g.set_bv('rpmi2',ts_tstream_type.float,bid,0.0)
	g.set_bv('rpmj1',ts_tstream_type.float,bid,0.0)	### set it so it rotates with subsequent row
	g.set_bv('rpmj2',ts_tstream_type.float,bid,0.0) ### ste it so it rotates with previous rotor
	g.set_bv('rpmk1',ts_tstream_type.float,bid,0.0)	### K SURFACES SPIN WITH BLOCK
	g.set_bv('rpmk2',ts_tstream_type.float,bid,0.0)
g.set_av('restart', ts_tstream_type.int,0)
g.set_av('nchange',ts_tstream_type.int,10000)
g.set_av('ga', ts_tstream_type.float,Gamma)
g.set_av('cp', ts_tstream_type.float,Cpgas)
g.set_av('cfl_ko', ts_tstream_type.float,0.4)
g.set_bv('nimixl', ts_tstream_type.int,2, 5)
g.set_av('nstep', ts_tstream_type.int, 500000)
#g.set_bv('dampin_mul', ts_tstream_type.float,1, 0.5)
'''tsr = ts_tstream_reader.TstreamReader()
g2 =tsr.read('OUPUT.hdf5')
for bid in g.get_block_ids():####
	b=g.get_block(bid)
	b2=g2.get_block(bid)
	if b.ni!=b2.ni:print 'Ni different',bid
	if b.nj!=b2.nj:print 'Nj different',bid
	if b.nk!=b2.nk:print 'Nk different',bid

	if b.ni!=b2.ni and b.nj==b2.nj and b.nk==b2.nk:
		for bp in ["ro","rovx","rovr","rorvt","roe","trans_dyn_vis"]:
			data2 = g2.get_bp(bp,bid)
			data = g.get_bp(bp,bid)
			data[:,:,:]=np.mean(data2)
			for i in range(b.ni):
				itemp = int(i*(b2.ni-1)/(b.ni-1))
				data[:,:,i]=data2[:,:,itemp]

			g.set_bp(bp, ts_tstream_type.float, bid, data)	
		print 'interpolating bp in i'
	elif b.ni==b2.ni and b.nj==b2.nj and b.nk!=b2.nk:
		for bp in ["ro","rovx","rovr","rorvt","roe","trans_dyn_vis"]:
			data2 = g2.get_bp(bp,bid)
			data = g.get_bp(bp,bid)
			data[:,:,:]=np.mean(data2)
			for k in range(b.nk):
				ktemp = int(k*(b2.nk-1)/(b.nk-1))
				data[k,:,:]=data2[ktemp,:,:]

			g.set_bp(bp, ts_tstream_type.float, bid, data)	
		print 'interpolating bp in k'
	elif b.ni!=b2.ni or b.nj!=b2.nj or b.nk!=b2.nk:
		for bp in ["ro","rovx","rovr","rorvt","roe","trans_dyn_vis"]:
			data2 = g2.get_bp(bp,bid)
			data = g.get_bp(bp,bid)
			data[:,:,:]=np.mean(data2)
			g.set_bp(bp, ts_tstream_type.float, bid, data)
	else:

		for bp in ["ro","rovx","rovr","rorvt","roe","trans_dyn_vis"]:
			g.set_bp(bp, ts_tstream_type.float, bid, g2.get_bp(bp,bid))'''

ts_tstream_load_balance.load_balance(g, 1)
fname='INPUT'
g.write_hdf5(fname+'.hdf5')
g.write_xdmf(fname+".xdmf",'x','r','rt')














