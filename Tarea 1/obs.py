import numpy as np 
import matplotlib.pyplot as plt 
import pandas as pd

# Primero leemos la carpeta en la que se encuentran los datos
path = r"C:\Users\Camila\Desktop\Astronomía Experimental\Tarea1AstroExperimental\sec_mierc_sem1_2021"

# definimos la función con la que se graficará
# la gaussiana que queremos graficar
def gauss(x,T0,mean,stdv):
    return T0*np.exp(-((x-mean)**2)/(2*(stdv**2)))
