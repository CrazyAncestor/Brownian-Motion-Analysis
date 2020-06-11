import numpy as np
import pandas as pd


#physical constant
dt = 1./30.
dl = 1e-4/305.#3.26e-6/4.
eta_water = 8.94e-4
eta_alcohol = 1.074e-3#3.06e-4


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

def rm0(r):
    u = []
    for i in range(len(r)):
        if r[i]!=0.:
            u.append(r[i])
    return np.array(u)

def calculate(particle):
    x = particle[:,1]
    y = particle[:,2]
    dx = x[1:]-x[:-1]
    dy = y[1:]-y[:-1]
    r2 = dx**2 + dy**2
    
    r2_ave = np.average(r2)
    r2_var = np.var(r2)

    dim = particle[:,4]/2.
    dim_ave = np.average(dim)
    dim_var = np.var(dim)
    return np.array([r2_ave*(dl**2), r2_var*(dl**2),dim_ave*(dl),dim_var*(dl)])

def analysis(solution,metal,time_elapse,T):
    
    filename = metal+str('-')+solution+str('.csv')

    data_source= pd.read_csv(filename)
    data = np.array(data_source)
    particles = split(data)
    kb = []

    n=0
    for i in range(len(particles)):
        [r2_ave, r2_var, dim_ave,dim_var] =calculate(particles[i])
        if r2_ave>0.0 and len(particles[i])>time_elapse:
            n+=1
            temp = 0
            if solution == 'water':
                temp = r2_ave*dim_ave*eta_water*(3*np.pi)/(2*T*dt)
            elif solution == 'alcohol':
                temp = r2_ave*dim_ave*eta_alcohol*(3*np.pi)/(2*T*dt)
            
            kb.append(temp)
            
    print(n)
    print(filename+':'+str(np.average(kb)))

    return np.average(kb)

if __name__ == '__main__':
   T = 273.15+26.2
   fewater= analysis('water','Fe',120,T)
   fealcohol = analysis('alcohol','Fe',25,T)
   eta_inv = np.array([1/eta_water,1/eta_alcohol])
   fig = np.array([fewater,fealcohol])
   print((fewater/eta_water-fealcohol/eta_alcohol)/(1/eta_water-1/eta_alcohol))
   

   T = 273.15+24.7
   mnwater= analysis('water','Mn',90,T)
   mnalcohol = analysis('alcohol','Mn',40,T)
   eta_inv = np.array([1/eta_water,1/eta_alcohol])
   fig = np.array([mnwater,mnalcohol])
   print((mnwater/eta_water-mnalcohol/eta_alcohol)/(1/eta_water-1/eta_alcohol))
   