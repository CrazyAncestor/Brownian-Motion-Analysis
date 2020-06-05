import numpy as np
import pandas as pd
debug = True
data_source= pd.read_csv("Fe-H2O-1-test.csv")#You may input the wanted file name here
data = np.array(data_source)

#physical constant
dt = 1./6.
dl = 3.26e-6
eta = 8.94e-4
T = 273.15+26.2

def split(df):
    n = 0
    prev = 0
    row = 0
    particles = []
    while row<len(df):
        if df[row][0]!=n:            
            particles.append(df[prev:row][:])
            prev = row
            n+=1
        row+=1
    return np.array(particles)

def calculate(particle):
    x = particle[:,1]
    y = particle[:,2]
    dx = x[1:]-x[:-1]
    dy = y[1:]-y[:-1]
    r2 = dx**2 + dy**2
    r2_ave = np.average(r2)
    r2_var = np.var(r2)

    dim = particle[:,4]
    dim_ave = np.average(dim)
    dim_var = np.var(dim)
    return np.array([r2_ave*(dl**2), r2_var*(dl**2),dim_ave*(dl),dim_var*(dl)])

if debug:
    particles = split(data)
    r2 = []
    a_inv = []
    kb = []
    for i in range(len(particles)):
        [r2_ave, r2_var, dim_ave,dim_var] =calculate(particles[i])
        if r2_ave!=0.0:
            r2.append(r2_ave)
            a_inv.append(1/dim_ave)
            kb.append( r2_ave*dim_ave*(3*np.pi*eta)/(2*T*dt))
            #print(r2_ave*dim_ave*(3*np.pi)*eta/(2*T*dt))
    r2 = np.array(r2)
    a_inv = np.array(a_inv)
    
    result = np.polyfit(a_inv,r2,1)
    print('kb:'+str(result[0]*(3*np.pi*eta)/(2*T*dt)))
    #print(np.average(kb))
