from astropy.io import fits
from astropy.stats import sigma_clip
from scipy.optimize import curve_fit
import numpy as np 
import pandas as pd
import scipy as sp
import matplotlib.pyplot as plt 
cubo = fits.open("southgal_fixbadc.fits") #abrir objeto cubo de datos
data = cubo[0].data #extraer matriz de datos
header = cubo[0].header #extraer el header del archivo fits

# funcion para obtener los valores de los ejes l, b, v
def values(h,j):
    N=h['NAXIS'+str(j)];
    val=np.zeros(N)
    for i in range(0,N):
        val[i] = (i+1-float(h['CRPIX'+str(j)]))*float(h['CDELT'+str(j)]) + float(h['CRVAL'+str(j)])
    return val

#Estos seran los tres arreglos con los valores reales de los tres ejes del cubo
velocidad=values(header,1)
longitud=values(header,2)
latitud=values(header,3)

columns = ['longitud l', 'latitud b', 'v_tan']
tabla = pd.DataFrame(columns=columns)

for i_b in range(len(latitud)):
    for i_l in range(len(longitud)):
        T = data[i_b][i_l][:]
        r = sigma_clip(T, sigma_lower=3, sigma_upper=3)
        rms = np.sqrt(np.mean(r**2))
        rmask = r.mask
        if len(velocidad[rmask])==0:
            v_tan = np.nan
        else:
            v_tan = velocidad[rmask][0]
        
        tabla = tabla.append({'longitud l':longitud[i_l], 'latitud b':latitud[i_b], 'v_tan':v_tan}, ignore_index=True) 

for lat in latitud:
    table_b_fix = tabla.loc[tabla['latitud b'] == lat]
    min_vel = table_b_fix['v_tan'].min()
    # print(lat, min_vel)

for lon in longitud:
    table_b_fix = tabla.loc[tabla['longitud l'] == lon]
    min_vel = table_b_fix['v_tan'].min()
    # print(lon, min_vel)

# Se crea una funcion que para una longitud(l) fija, se recorre latitud(b) y se calcula el rms de las
# velocidades
# Esta misma funcion recorre el cubo de las velocidades asociadas a l y b, hasta que se llega a una
# velocidad que es 5 veces mayor que el rms, esta ultima se guarda un arreglo

def fmin(l,latitud,vs):
    #recorre latitud
    for q in range(33):
        T1=data[q][l][:]
        rms=np.sqrt(np.mean(T1**2))   #calcula rms
        #recorre velocidad
        for w in range(306):
            if data[q][l][w]>=5*rms:  #buscamos que no sea ruido
                vs[q]=velocidad[w]    #guardamos la primera v donde T mayor a 5rms
                break

vmin = np.zeros(385)
bvmin = np.zeros(385)
R = np.zeros(385)
Z = np.zeros(385)
R0km = 2.623e+17 # distancia en km
R0 = 8.5 #kPc
vsol = 220
omegasol = vsol/R0km

#maximorum
# Se recorren las longitudes y se busca la velocidad más negativa (mayor en modulo), se guarda esta
# y su latitud asociada
# Se obtiene un arreglo de R con la ecuacion R =| R0 · cos(l π/180 ) |

for i in range(385):
    vs=np.zeros(33)
    fmin(i,latitud,vs)
    v1=vs[0]
    b1=latitud[0]
    for j in range(32):
        if vs[j+1]<v1:
            v1=vs[j+1]
            b1=latitud[j+1]
    vmin[i]=v1
    bvmin[i]=b1
    R[i]=np.abs(R0*np.sin(longitud[i]*np.pi/180.)) #R0 sin(l)  
    Z[i]= b1*np.pi/180*R0*np.cos(longitud[i]*np.pi/180.)
    

    
# Se obtiene la Vtan con Vtan = −Vmin − Vsol · sin(lπ/180 ), donde Vmin es la velocidad mayor en
# modulo para l, y Vsol es la velocidad de rotacion del sol.    
#velocidad de rotacion

vR = np.zeros(385)
for i in range(385):
    vR[i] = vmin[i]*(np.abs(np.sin(longitud[i]*np.pi/180.))/ \
            np.sin(longitud[i]*np.pi/180.)) + \
            np.abs(vsol*np.sin(longitud[i]*np.pi/180.))

omegaR = np.zeros(385)
for i in range(385):
    # 1 kpc = 3.08567758128E+16 km
    # 1 km = 3.2408e-17 kpc
    # 220 km/s = 220
    omegaR[i] = (vR[i]*3.2408e-17)/R[i] + omegasol
'''
CURVA DE ROTACIÓN
'''
plt.plot(R,vR, 'lightcoral')
plt.plot(R,vR, 'k.')
plt.grid
plt.title("Curva de Rotación")
plt.xlabel("R [kpc]")
plt.ylabel(r"V$_{tan}$ [km/s]")
# plt.savefig('curva.png')
plt.show()

plt.plot(R, omegaR, 'lightcoral')
plt.plot(R, omegaR, 'k.')
plt.grid
plt.title("Curva de Rotación")
plt.xlabel("R [kpc]")
plt.ylabel(r"$\omega_{tan}$ [rad/s]")
plt.savefig('curva2.png')
plt.show()



'''
AJUSTE DE MASAS
'''
def masa_puntual(R, M0):
    return np.sqrt(G*M0/R)

def esfera_uniforme(R, rho):
    M = 4/3*np.pi*R**3*rho
    v = np.sqrt(G*M/R)
    return v

def esfera_uniforme_masapuntual(R, rho, M0):
    M = 4/3*np.pi*R**3*rho + M0
    v = np.sqrt(G*M/R)
    return v

def disco_uniforme(R, s):
    M = np.pi*R**2*s
    v = np.sqrt(G*M/R)
    return v

def disco_uniforme_masapuntual(R, s, M0):
    M = np.pi*R**2*s + M0
    v = np.sqrt(G*M/R)
    return v

G = 4.302e-6   # kpc/Msun*(km**2/s**2)


m1,covm1 = curve_fit(masa_puntual, R, vR)
m2,covm2 = curve_fit(esfera_uniforme, R, vR)
m3,covm3 = curve_fit(esfera_uniforme_masapuntual, R, vR)
m4,covm4 = curve_fit(disco_uniforme, R, vR)
m5,covm5 = curve_fit(disco_uniforme_masapuntual, R, vR)

print('M0 masa puntual:', m1[0])
print('rho esfera uniforme:', m2[0])
print('rho, M0 esfera unirme con masa puntual:', m3)
print('s disco uniforme:', m4[0])
print('s, M0 disco uniforme masa puntual:', m5)

fig = plt.figure(figsize= (5, 10))
ax1 = fig.add_subplot(5,1,1)
ax2 = fig.add_subplot(5,1,2, sharex = ax1)
ax3 = fig.add_subplot(5,1,3, sharex = ax1)
ax4 = fig.add_subplot(5,1,4, sharex = ax1)
ax5 = fig.add_subplot(5,1,5, sharex = ax1)


ax1.plot(R,vR,'lightcoral',label='datos')
# ax1.plot(R,vR,'r-',label='datos')
ax1.plot(R,masa_puntual(R,m1[0]),'k',label='masa puntual')

ax2.plot(R,vR,'lightcoral',label='datos')
# ax2.plot(R,vR,'r-',label='datos')
ax2.plot(R,esfera_uniforme(R,m2[0]),'k',label='esfera uniforme')

ax3.plot(R,vR,'lightcoral',label='datos')
# ax3.plot(R,vR,'r-',label='datos')
ax3.plot(R,esfera_uniforme_masapuntual(R,m3[0], m3[1]),'k',label='esfera uniforme mp')

ax4.plot(R,vR,'lightcoral',label='datos')
# ax4.plot(R,vR,'r-',label='datos')
ax4.plot(R,disco_uniforme(R,m4[0]),'k',label='disco uniforme')

ax5.plot(R,vR,'lightcoral',label='datos')
# ax5.plot(R,vR,'r-',label='datos')
ax5.plot(R,disco_uniforme_masapuntual(R,m5[0], m5[1]),'k',label='disco uniforme mp')

ax1.set_title('Ajuste de Masas')
ax5.set_xlabel('R [kpc]')
ax3.set_ylabel('Velocidad [km/s]')
ax1.grid()
ax2.grid()
ax3.grid()
ax4.grid()
ax5.grid()
ax1.legend()
ax2.legend()
ax3.legend()
ax4.legend()
ax5.legend()
plt.subplots_adjust(wspace=0, hspace=0)
# plt.savefig('ajuste.png')
plt.show()

'''
CORRUGACIÓN DEL PLANO
'''
plt.plot(R,Z, 'lightcoral')
plt.plot(R,Z, 'k.')
plt.grid()
plt.title("Corrugación del plano en función de R")
plt.xlabel("R [kpc]", fontsize='12')
plt.ylabel("Z [kpc]", fontsize='12')
plt.tick_params(labelsize='12')
# plt.savefig('corrugacion.png')
plt.show()

# %%
