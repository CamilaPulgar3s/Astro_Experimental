# Como solo tenemos que trabajar con 4 valores, se definen como constantes

T_H = 300
T_C = 77
P_H = -44.5 #dBm
P_C = -47.94 #dBm

# Las temperaturas están en kelvin por lo que no requieren 
# un cambio en la unidad de medida
# Luego podemos definir la función que calcule T_rec
def T_rec(y_factor, T_hot, T_cold):
    return (T_hot - y_factor*T_cold)/(y_factor - 1)

# Las potencias están en dBm, por lo que se debe crear una función 
# que las pase a Watts
def dBm_to_watt(P):
    return (10**(P/10))/1000

# Con esto ya podemos calcular W_H y W_C para luego
# obtener el Y factor, con el que podremos calcular el T_rec deseado
W_H = dBm_to_watt(P_H)
W_C = dBm_to_watt(P_C)
Y = W_H/W_C

T_rec = T_rec(Y, T_H, T_C)
print('La temperatura de ruido obtenida experimentalmente es',T_rec)
print('El valor de la potencia W_hot en Watts es', W_H) 
print('El valor de la potencia W_cold en Watts es', W_C) 
print('Por lo tanto, el valor del factor Y es', Y)