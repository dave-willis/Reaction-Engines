#import itertools
import numpy as np

#n = len(np.arange(5, 11))
#phi = len(np.arange(0.2, 0.9, 0.3))
#psi = len(np.arange(0.6, 2, 0.3))
#Lambda = len(np.arange(0.2,0.9,0.3))
#AR = len(np.arange(0.9, 2.0, 0.5))
#dho = len(np.arange(1, 1.2, 0.1))
#a1 = len(np.arange(0, 20, 2))
#
#num = AR*phi**2*psi**2*dho**2*Lambda**2
#
#print(num)


import csv
import scipy.interpolate as sciint

table = []
xin = []
xout = []
sslen = []
with open("ss_length.csv") as csvfile:
    reader = csv.reader(csvfile, quoting=csv.QUOTE_NONNUMERIC) # change contents to floats
    for row in reader: # each row is a list
        table.append(row)
        xin.append(row[0])
        xout.append(row[1])
        sslen.append(row[2])
        
xout.insert(1,'')
        
with open('ss_grid.csv', mode='w') as ss_grid:
    table_writer = csv.writer(ss_grid, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    table_writer.writerow(xout[1:81])
    for i in range(1, int(len(xin)/80)):
        row = sslen[80*i+1:80*(i+1)]
        row.insert(0, xin[80*i])
        table_writer.writerow(row)

table = []       
with open("ss_grid.csv") as csvfile:
    reader = csv.reader(csvfile, quoting=csv.QUOTE_NONNUMERIC) # change contents to floats
    for row in reader: # each row is a list
        table.append(row)

for i in range(1, len(table)):
    for j in range(1, len(table[0])):
        if table[i][0] > table[0][j]:
            table[i][j] = table[len(table)-i][len(table[0])-j]
            
with open('ss_grid.csv', mode='w') as ss_grid:
    table_writer = csv.writer(ss_grid, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    for row in table:
        table_writer.writerow(row)
        
#with open("ss_grid.csv") as csvfile:
#    reader = csv.reader(csvfile, quoting=csv.QUOTE_NONNUMERIC) # change contents to floats
#    line = 0
#    for row in reader: # each row is a list
#        table.append(row)
#        if line == 0:
#            xout = row[1:]
#            line += 1
#        else:
#            xin.append(row[0])
#            sslen.append(row[1:])
#            line += 1
#
#xout = np.asarray(xout)
#xin = np.asarray(xin)
#sslen = np.asarray(sslen)
#print(xout.size)
#print(xin.size)
#print(sslen.shape)    
#import time
#start = time.time()       
#f1 = sciint.interp2d(xout, xin, sslen)
#print(time.time()-start)
#print(f1(69.67,31.345))
#print(time.time()-start)
#start = time.time()  
#f2 = sciint.RectBivariateSpline(xin, xout, sslen, kx=1, ky=1)
#print(time.time()-start)
#print(f2(31.345,69.67)[0][0])
#print(time.time()-start)
