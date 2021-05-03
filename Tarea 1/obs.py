import numpy as np 
import matplotlib.pyplot as plt
import glob
from scipy.optimize import curve_fit


path = r"C:\Users\Camila\Desktop\Astronomía Experimental\Tarea1AstroExperimental\sec_mierc_sem1_2021"

# Se define funcion gaussiana a fitear
def f_gauss(x,T0,mean,stdv):
    return T0*np.exp(-((x-mean)**2)/(2*(stdv**2)))

fl = sorted(glob.glob(path+'\sdf*'))

#print(fl)
for i in range(len(fl)):
    v,T = np.genfromtxt(fl[i], unpack = True, skip_header=108)
    # First Guess: T0=20, m=10, s=1 obtenidos al inspeccionar espectro ploteado
    fg = [20, 10, 1] 
    coefs,cov = curve_fit(f_gauss,v,T, p0=fg) # Se fitea
    t0,M,S = coefs[0],coefs[1],coefs[2]  # Se extraen los coeficientes fiteados
    # print ('Valores fiteados: (t0,M,S) =',(t0,M,S))

    plt.plot(v,T, label='Data', color = 'blue')
    plt.plot(v,f_gauss(v,t0, M,S), label='Fiteo', color = 'red')
    plt.title('Espectro '+fl[i][-11:], fontsize=18)
    plt.ylabel('Temperatura [K]', fontsize=18)
    plt.xlabel(r'Velocidad [$\frac{km}{s}$]', fontsize=18)
    plt.legend()
    # plt.savefig('Espectro '+fl[i][-11:]+'.pdf')