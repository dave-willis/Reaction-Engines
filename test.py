import numpy as np

#sigma_y = 300*10**6
#rho = 8300
#omega = 6782*2*np.pi/60
#Ro = 0.5
#nu = 0.3
#P = 145*10**5
#
#a = rho*omega**2*(1-nu)
#b = 2*(rho*omega**2*Ro**2*(1+nu)-2*sigma_y)
#c = 4*sigma_y*Ro**2-8*P*Ro**2-rho*omega**2*(3+nu)*Ro**4
#
#dis = 4*np.sqrt(sigma_y**2+rho**2*omega**4*Ro**4+2*rho*omega**2*Ro**2*(P*(1-nu)-sigma_y))
#
#print('a:',a)
#print('b:',b)
#print('c:',c)
#print('discriminant:',sigma_y**2+rho**2*omega**4*Ro**4+2*rho*omega**2*Ro**2*(P*(1-nu)-sigma_y))
#
#print('Dis. check:',dis-np.sqrt(b**2-4*a*c))
#
#Ri = np.sqrt((-b-np.sqrt(b**2-4*a*c))/(2*a))
#
#print('Ri:',Ri)
#
#sigma = rho*omega**2/4*((3+nu)*Ro**2+(1-nu)*Ri**2)+2*P*Ro**2/(Ro**2-Ri**2)
#
#print('Calc. check:',sigma-sigma_y)
#
#Ri=0.21
#
#sig_P = 2*P*Ro**2/(Ro**2-Ri**2)
#
#sig_r = rho*omega**2/4*((3+nu)*Ro**2+(1-nu)*Ri**2)
#
#print('Yield stress:',"{:.2e}".format(sigma_y))
#print('Pressure stress:',"{:.2e}".format(sig_P))
#print('Rotation stress:',"{:.2e}".format(sig_r))
#print('Stress sum:',"{:.2e}".format(abs(sig_r+sig_P)))
#Ri=0.6
#print(np.sqrt(sigma_y*Ri**2/(sigma_y-2*P)))

from Profile import blade_dims

print(blade_dims(0,70,0.01,0.05)[0], blade_dims(0,70,0.01,1)[0]*0.05**2)
print(blade_dims(0,70,0.01,0.05)[1], blade_dims(0,70,0.01,1)[1]*0.05)

#
#def radial_stress(Ri):
#    
#    test_max = Ro**2*Ri**2*(1-8*P/(rho*omega**2*(3+nu)*(Ro**2-Ri**2)))
#  
#    if test_max > 0:
#        r = test_max**0.25
#    else:
#        r = Ro
#    
#    sigma_rot = rho*omega**2/8*(Ro**2+Ri**2-Ro**2*Ri**2/r**2-r**2)
#    
#    sigma_p = P*Ro**2*(r**2-Ri**2)/(r**2*(Ro**2-Ri**2))
#    
#    sigma_rad = sigma_rot-sigma_p
#   
#    return sigma_y-abs(sigma_rad)
#
#
#def tan_stress(Ri):
#
#    test_max = Ro**2*Ri**2*(8*P/(rho*omega**2*(1+3*nu)*(Ro**2-Ri**2))-(3+nu)/(1+3*nu))
#    test_Ri = rho*omega**2/4*((3+nu)*Ro**2+(1-nu)*Ri**2)-2*P*Ro**2/(Ro**2-Ri**2)
#    test_Ro = rho*omega**2/4*((3+nu)*Ri**2+(1-nu)*Ro**2)-2*P*(Ro**2+Ri**2)/(Ro**2-Ri**2)
#  
#    if test_max > 0:
#        r = test_max**0.25
#        
#        sigma_rot = rho*omega**2/8*((3+nu)*(Ro**2+Ri**2+Ro**2*Ri**2/r**2)-(1+3*nu)*r**2)
#        
#        sigma_p = P*Ro**2*(r**2+Ri**2)/(r**2*(Ro**2-Ri**2))
#        
#        sigma_tan = sigma_rot-sigma_p
#        
#    elif abs(test_Ri) > abs(test_Ro):
#        
#        sigma_tan = test_Ri
#        
#    else:
#        
#        sigma_tan = test_Ro
#    
#    return sigma_y-abs(sigma_tan)
#
#def Radius_i(Ri):
#    return -Ri
#
#from scipy.optimize import minimize
##Form the starting point list
#x0 = [0.99*Ro]
##Form the tuple of bounds
#bnds = [(0,Ro)]
##Form the tuple of constraints
#cons = cons = ({'type': 'ineq', 'fun': tan_stress}, {'type': 'ineq', 'fun': radial_stress})
##Find the minimum
#res = minimize(Radius_i, x0, method='SLSQP', bounds=bnds, constraints=cons)
##Extract the optimal variable and return them
#print(res['x'])
#print(tan_stress(res['x'])+sigma_y)
#print(radial_stress(res['x'])+sigma_y)
#
#Ri = res['x']
#
#print(Ro**2*Ri**2*(8*P/(rho*omega**2*(1+3*nu)*(Ro**2-Ri**2))-(3+nu)/(1+3*nu)))
