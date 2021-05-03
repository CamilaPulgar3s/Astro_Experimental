import numpy as np 
import matplotlib.pyplot as plt
import glob
from scipy.optimize import curve_fit
from scipy.integrate import trapz


path = r"C:\Users\Camila\Desktop\Astronomía Experimental\Tarea1AstroExperimental\sec_mierc_sem1_2021"

# Se define funcion gaussiana a fitear
def f_gauss(x,T0,mean,stdv):
    return T0*np.exp(-((x-mean)**2)/(2*(stdv**2)))

fl = sorted(glob.glob(path+'\sdf*'))

Temperaturas = []

lii = []
bii = []

#print(fl)
for i in range(len(fl)):

    file = open(fl[i], "r")
    lines = file.readlines()
    lii.append(float(lines[22].strip()[5:]))
    bii.append(float(lines[23].strip()[5:]))

    v,T = np.genfromtxt(fl[i], unpack = True, skip_header=108)

    Temperaturas.append(T)
    # First Guess: T0=20, m=10, s=1 obtenidos al inspeccionar espectro ploteado
    fg = [20, 10, 1] 
    coefs,cov = curve_fit(f_gauss,v,T, p0=fg) # Se fitea T
    t0,M,S = coefs[0],coefs[1],coefs[2]  # Se extraen los coeficientes fiteados
    # print ('Valores fiteados: (t0,M,S) =',(t0,M,S))

    plt.plot(v,T, label='Data', color = 'blue')
    plt.plot(v,f_gauss(v,t0, M,S), label='Fiteo', color = 'red')
    plt.title('Espectro '+fl[i][-11:], fontsize=18)
    plt.ylabel('Temperatura [K]', fontsize=18)
    plt.xlabel(r'Velocidad [$\frac{km}{s}$]', fontsize=18)
    plt.legend()
    # plt.savefig('Espectro '+fl[i][-11:]+'.pdf')

T_max1 = (Temperaturas[0].max() + Temperaturas[5].max() + Temperaturas[10].max())/3 #arriba
T_max2 = (Temperaturas[1].max() + Temperaturas[6].max() + Temperaturas[11].max())/3 #izquierda
T_max3 = (Temperaturas[2].max() + Temperaturas[7].max() + Temperaturas[12].max())/3 #centro
T_max4 = (Temperaturas[3].max() + Temperaturas[8].max() + Temperaturas[13].max())/3 #derecha
T_max5 = (Temperaturas[4].max() + Temperaturas[9].max() + Temperaturas[14].max())/3 #abajo

T_maxh = np.array([T_max2, T_max3, T_max4])
T_maxv = np.array([T_max1, T_max3, T_max5])

coordlii = np.array([lii[1], lii[2], lii[3]])
coordbii = np.array([bii[0], bii[2], bii[4]])
# 208.9,-19.2   208.8,-19.3   208.9,-19.3  209.1,-19.3 208.9, -19.5

# el centro de la nebulosa tiene mayor temperatura
fg1 = [T_max3, 208.9, 1] 
fg2 = [T_max3, -19.3, 1]
coefs1,cov1 = curve_fit(f_gauss,coordlii,T_maxv, p0=fg1) # Se fitea T
coefs2,cov2 = curve_fit(f_gauss,coordbii,T_maxh, p0=fg2) # Se fitea T
t01,M1,S1 = coefs1[0],coefs1[1],coefs1[2]
t02,M2,S2 = coefs2[0],coefs2[1],coefs2[2]  

clii = np.linspace(208.8, 209.1, 1000)
cbii = np.linspace(-20, -19, 1000)

plt.clf()
plt.plot(coordlii,T_maxv, label='Data', color = 'blue', marker = 'o', fillstyle = 'none', ls = '')
plt.plot(clii,f_gauss(clii,t01, M1,S1), label='Fiteo', color = 'red')
plt.title(r'$T_{max}$ v/s l$^{II}$, b$^{II}$ = -19°', fontsize=18)
plt.ylabel('Temperatura máxima [K]', fontsize=18)
plt.xlabel(r'Longitud [°]', fontsize=18)
plt.legend()
plt.savefig('lii.pdf')

plt.clf()
plt.plot(coordbii,T_maxv, label='Data', color = 'blue', marker = 'o', fillstyle = 'none', ls = '')
plt.plot(cbii,f_gauss(cbii,t02, M2,S2), label='Fiteo', color = 'red')
plt.title(r'$T_{max}$ v/s b$^{II}$, l$^{II}$ = 208°', fontsize=18)
plt.ylabel('Temperatura máxima [K]', fontsize=18)
plt.xlabel(r'Latitud [°]', fontsize=18)
plt.legend()
plt.savefig('bii.pdf')

