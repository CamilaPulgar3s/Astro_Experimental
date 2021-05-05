import numpy as np 
import matplotlib.pyplot as plt
import glob
from scipy.optimize import curve_fit
from scipy.integrate import trapz
from astropy.stats import sigma_clip


path = r"C:\Users\Camila\Desktop\Astronomía Experimental\Tarea1AstroExperimental\sec_mierc_sem1_2021"

# Se define funcion gaussiana a fitear
def f_gauss(x,T0,mean,stdv):
    return T0*np.exp(-((x-mean)**2)/(2*(stdv**2)))

fl = sorted(glob.glob(path+'\sdf*'))

# creamos las listas en las que iremos guardando los valores para 
# graficar en las partes b y c
Temperaturas = []
T_integrada = []
lii = []
bii = []

#print(fl)
for i in range(len(fl)):
    # Para cada archivo dentro de la carpeta del path 
    # sacamos la información que necesitamos 
    # que corresponden a las coordenadas, las velocidades y temperaturas
    file = open(fl[i], "r")
    lines = file.readlines()
    lii.append(float(lines[22].strip()[5:]))
    bii.append(float(lines[23].strip()[5:]))
    v,T = np.genfromtxt(fl[i], unpack = True, skip_header=108)
    
    # agregamos las temperaturas a la lista previamente creada
    Temperaturas.append(T)
    # las temperaturas son decrecientes, entonces hay que poner un signo menos 
    # antes para que las temperaturas integradas sean positivas
    T_integrada.append(-trapz(T,v))
    
    # First Guess: T0=20, m=10, s=1 obtenidos al inspeccionar espectro ploteado
    fg = [20, 10, 1] 
    coefs,cov = curve_fit(f_gauss,v,T, p0=fg) # Se fitea T
    t0,M,S = coefs[0],coefs[1],coefs[2]  # Se extraen los coeficientes fiteados
    # print ('Valores fiteados: (t0,M,S) =',(t0,M,S))
    
    # el siguiente plot graficará los 15 espectros con su gaussiana
    plt.clf()
    plt.plot(v,T, label='Data', color = 'blue')
    plt.plot(v,f_gauss(v,t0, M,S), label='Fiteo', color = 'red')
    plt.title('Espectro '+fl[i][-11:], fontsize=15)
    plt.ylabel('Temperatura [K]', fontsize=12)
    plt.xlabel(r'Velocidad [$\frac{km}{s}$]', fontsize=12)
    plt.legend()
    # plt.savefig('Espectro '+fl[i][-11:]+'.pdf')

# ahora queremos fitear las temperaturas máximas versus coordenadas 
# primero calculamos la temperatura máxima promedio de cada cuadrante 
T_max1 = (Temperaturas[0].max() + Temperaturas[5].max() + Temperaturas[10].max())/3 #arriba
T_max2 = (Temperaturas[1].max() + Temperaturas[6].max() + Temperaturas[11].max())/3 #izquierda
T_max3 = (Temperaturas[2].max() + Temperaturas[7].max() + Temperaturas[12].max())/3 #centro
T_max4 = (Temperaturas[3].max() + Temperaturas[8].max() + Temperaturas[13].max())/3 #derecha
T_max5 = (Temperaturas[4].max() + Temperaturas[9].max() + Temperaturas[14].max())/3 #abajo

# diferenciamos las temperaturas máximas de los cuadrantes horiontales y verticales 
T_maxh = np.array([T_max2, T_max3, T_max4])
T_maxv = np.array([T_max1, T_max3, T_max5])

# creamos arreglos de las coordenadas de cada cuadrante 
coordlii = np.array([lii[1], lii[2], lii[3]])
coordbii = np.array([bii[0], bii[2], bii[4]])
# 208.9,-19.2   208.8,-19.3   208.9,-19.3  209.1,-19.3 208.9, -19.5

# como el centro de la nebulosa tiene mayor temperatura
# la temperatura máxima que recibe como parámetro la función curve_fit es T_max3
fg1 = [T_max3, 208.9, 1] 
fg2 = [T_max3, -19.3, 1]
coefs1,cov1 = curve_fit(f_gauss,coordlii,T_maxv, p0=fg1) # Se fitea T
coefs2,cov2 = curve_fit(f_gauss,coordbii,T_maxh, p0=fg2) # Se fitea T
t01,M1,S1 = coefs1[0],coefs1[1],coefs1[2]
t02,M2,S2 = coefs2[0],coefs2[1],coefs2[2]  

# ejes x para graficar:
clii = np.linspace(207, 210, 1000)
cbii = np.linspace(-20, -19, 1000)

# gráfico tmax vs lii
plt.clf()
plt.plot(coordlii,T_maxv, label='Data', color = 'blue', marker = 'o', fillstyle = 'none', ls = '')
plt.plot(clii,f_gauss(clii,t01, M1,S1), label='Fiteo', color = 'red')
plt.title(r'$T_{max}$ v/s l$^{II}$, b$^{II}$ = -19°', fontsize=15)
plt.ylabel('Temperatura máxima [K]', fontsize=12)
plt.xlabel(r'Longitud [°]', fontsize=12)
plt.legend()
# plt.savefig('lii.pdf')

# gráfico tmax vs bii
plt.clf()
plt.plot(coordbii,T_maxv, label='Data', color = 'blue', marker = 'o', fillstyle = 'none', ls = '')
plt.plot(cbii,f_gauss(cbii,t02, M2,S2), label='Fiteo', color = 'red')
plt.title(r'$T_{max}$ v/s b$^{II}$, l$^{II}$ = 208°', fontsize=15)
plt.ylabel('Temperatura máxima [K]', fontsize=12)
plt.xlabel(r'Latitud [°]', fontsize=12)
plt.legend()
# plt.savefig('bii.pdf')


# integral 
# Análogo al caso anterior, en vez de usar T_max se usa T_integrada
T_int1 = (T_integrada[0] + T_integrada[5] + T_integrada[10])/3 #arriba
T_int2 = (T_integrada[1] + T_integrada[6] + T_integrada[11])/3 #izquierda
T_int3 = (T_integrada[2] + T_integrada[7] + T_integrada[12])/3 #centro
T_int4 = (T_integrada[3] + T_integrada[8] + T_integrada[13])/3 #derecha
T_int5 = (T_integrada[4] + T_integrada[9] + T_integrada[14])/3 #abajo

T_inth = np.array([T_int2, T_int3, T_int4])
T_intv = np.array([T_int1, T_int3, T_int5])

# el centro de la nebulosa tiene mayor temperatura
fg3 = [T_int3, 208.9, 1] 
fg4 = [T_int3, -19.3, 1]
coefs3,cov3 = curve_fit(f_gauss,coordlii,T_intv, p0=fg3) # Se fitea T
coefs4,cov4 = curve_fit(f_gauss,coordbii,T_inth, p0=fg4) # Se fitea T
t03,M3,S3 = coefs3[0],coefs3[1],coefs3[2]
t04,M4,S4 = coefs4[0],coefs4[1],coefs4[2]  

# gráfico T integrada vs lii
plt.clf()
plt.plot(coordlii,T_intv, label='Data', color = 'blue', marker = 'o', fillstyle = 'none', ls = '')
plt.plot(clii,f_gauss(clii,t03, M3,S3), label='Fiteo', color = 'red')
plt.title(r'$T_{integrada}$ v/s l$^{II}$, b$^{II}$ = -19°', fontsize=15)
plt.ylabel(r'Temperatura integrada [$\frac{km}{s}$]', fontsize=12)
plt.xlabel(r'Longitud [°]', fontsize=12)
plt.legend()
# plt.savefig('tintlii.pdf')

# gráfico T integrada vs bii
plt.clf()
plt.plot(coordbii,T_intv, label='Data', color = 'blue', marker = 'o', fillstyle = 'none', ls = '')
plt.plot(cbii,f_gauss(cbii,t04, M4,S4), label='Fiteo', color = 'red')
plt.title(r'$T_{integrada}$ v/s b$^{II}$, l$^{II}$ = 208°', fontsize=15)
plt.ylabel(r'Temperatura integrada [$\frac{km}{s}$]', fontsize=12)
plt.xlabel(r'Latitud [°]', fontsize=12)
plt.legend()
# plt.savefig('tintbii.png')

# Error
# creamos un arreglo de temperaturas 
T = np.array(Temperaturas)
# calculamos el error en las 15 mediciones
rms = np.sqrt(T.mean(axis=1)**2)
# calculamos el error promedio en cada uno de los 5 cuadrantes
mean = (T[0:5]+T[5:10]+T[10:15])/3
rmsmean = np.sqrt(mean.mean(axis=1)**2)
# printeamos la razón existente que debe ser aprox 1/sqrt(3)
print(np.repeat(rmsmean,3) / rms)