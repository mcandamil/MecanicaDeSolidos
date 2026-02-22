---
layout: default
title: Librerías para el Estudio de la Elasticidad Plana en Python y Definición de Variables Simbólicas
---

```python
#Librería Simbólica SymPy
import sympy as sp

#Librería Numérica Numpy
import numpy as np

#Mostrar Mejor las Ecuaciones
from IPython.display import display

#Formatear Prints con Símbolos
sp.init_printing()

  

#Definición Simbólica de los Esfuerzos Axiales
sx, sy, sz = sp.symbols('sigma_x, sigma_y, sigma_z')

#Definición Simbólica de los Esfuerzos Tangenciales
txy, txz, tyz = sp.symbols('tau_xy, tau_xz, tau_yz')

#Definición Simbólica de las Deformaciones Longitudinales
ex, ey, ez = sp.symbols('epsilon_x, epsilon_y, epsilon_z')

#Definición Simbólica de las Deformaciones Angulares
gx, gy, gz = sp.symbols('gamma_x, gamma_y, gamma_z')

#Definición Simbólica de Módulos de Elasticidad y Coeficiente de Poisson
E, G, nu = sp.symbols('E, G, nu')
```
