import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 

df_AntDip = pd.read_excel(r"C:\Users\Camila\Desktop\Astronomía Experimental\Tarea1AstroExperimental\antdip_AE2021A.xlsx", comment = '#')


# Queremos graficar una función lineal y = mx + n
# donde y = ln(dW), x = 1/sen(E_l), m = tau0
# y n = ln(T_hot) = T(300) considerando la temperatura ambiente en Kelvin
# Leemos la columna 3 y 7 del archivo antdip_AE2021A.xlsx que contienen
# los valores de 1/sen(E_l) y ln(dW) respectivamente

ln_dP = df_AntDip.iloc[:-1,7]
sec = -(df_AntDip.iloc[:-1,3]).to_numpy() 

# Queremos hacer un ajuste lineal a los datos
# Para eos definimos una función que reciba x, m y n, 
# y nos entregue la función lineal y = mx + n
def lineal(x, m, n):
    x=np.array(x)
    return m*x + n

# Una vez definidas las variables y la función lineal
# se usa la función polyfit de numpy que entregará los valores 
# de la pendiente m y de n, luego estos se usan para definir 
# y_lineal que corresponde al ajuste lineal.
z = np.polyfit(sec, ln_dP, 1)
x_lineal = [np.min(sec), np.max(sec)]
y_lineal = lineal(x_lineal, z[0], z[1])

# Graficamos el ajuste lineal y los datos obtenidos experimentalmente
plt.plot(x_lineal, y_lineal, 'k', color = 'red', label = 'Ajuste Lineal')
plt.plot(sec, ln_dP, marker ='o', ls = '', label = 'Datos Experimentales')
plt.title(r'Fiteo $\tau_w$')
plt.xlabel(r'$\frac{1}{sen(E_l)}$')
plt.legend()
plt.show()
print(z)