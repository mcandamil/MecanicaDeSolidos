---
title: Deformación Plana en Materiales Isótropos
---

La deformación plana es el concepto de trabajar con solo una deformación longitudinal nula, permitiendo realizar un análisis del elemento en una sola rebanada de toda una sección prismática. Ocurre cuando existe una dimensión mucho más larga que las demás haciendo el elemento alargado.
Esta deformación plana sucede cuando:
- No existen cargas aplicadas al eje z.
- Las cargas aplicadas son perpendiculares al eje z.
- Las cargas están distribuidas uniformemente en la longitud del elemento.

<div align="center">  
<img src="PlaneStrain.png" width="450">  
</div>
<p align="center" style="font-size: 0.85rem; opacity: 0.7;">
  <b>Figura 1.</b> Elemento Sólido en Estado de Deformación Plana.<br>
  Fuente: Álvarez Marín, D. A. (2023). 
  <i>
    <a href="https://bffrepositorio.unal.edu.co/server/api/core/bitstreams/b07bc9a3-7305-4e52-97f3-108813d014a2/content" target="_blank">
      Teoría de la elasticidad usando Matlab y Maxima. Volumen 1: fundamentos.
    </a>
  </i> Universidad Nacional de Colombia.
</p>
Matemáticamente se representa como la nulidad de deformaciones en la dirección del eje z:
$$
\varepsilon_z = \gamma_{xz} = \gamma_{yz} = 0
$$
De esta manera, el tensor de esfuerzos en estado de deformación plana es:
$$
\underline{\underline{\varepsilon}} =
\begin{bmatrix}
\varepsilon_x & \gamma_{xy} & 0 \\
\gamma_{xy} & \varepsilon_y & 0 \\
0 & 0 & 0
\end{bmatrix}
$$
<span style="opacity:0.7; font-size:0.85rem;"><em>Nótese que el rango del tensor disminuye ya que tiene una columna sin aporte dimensional. Esto hace que el espacio donde trabaje sea máximo &#8477&#178 (cuando &epsilon;<sub>x</sub>&epsilon;<sub>y</sub> &minus; &gamma;<sub>xy</sub><sup>2</sup> &ne; 0), &#8477&#185 (cuando &epsilon;<sub>x</sub>&epsilon;<sub>y</sub> &minus; &gamma;<sub>xy</sub><sup>2</sup> = 0) o &#8477&#8304 (cuando todas las deformaciones son nulos). Pero en general, se está llevando el problema tridimensional a uno bidimensional para figuras prismáticas.
</em>
</span>

Cabe recordar las relaciones entre esfuerzos y deformaciones encontradas haciendo uso del módulo de elasticidad lineal $E$ y el coeficiente de Poisson $\nu$:
$$
\begin{alignedat}{2}

\varepsilon_x &= \frac{1}{E}\left(\sigma_x - \nu(\sigma_y + \sigma_z)\right)
\qquad\qquad
&\gamma_{xy} &= \frac{2(1+\nu)}{E}\,\tau_{xy} \\[12pt]

\varepsilon_y &= \frac{1}{E}\left(\sigma_y - \nu(\sigma_x + \sigma_z)\right)
\qquad\qquad
&\gamma_{yz} &= \frac{2(1+\nu)}{E}\,\tau_{yz} \\[12pt]

\varepsilon_z &= \frac{1}{E}\left(\sigma_z - \nu(\sigma_x + \sigma_y)\right)
\qquad\qquad
&\gamma_{xz} &= \frac{2(1+\nu)}{E}\,\tau_{xz}

\end{alignedat}
$$

Al igualar $\varepsilon_z = 0$ se observa que el módulo de elasticidad se anula en la ecuación y queda que $\sigma_z = \nu(\sigma_x+\sigma_y)$. Con esta igualdad, es posible reemplazar $\sigma_z$ en las expresiones de $\varepsilon_x$ y $\varepsilon_y$ :
$$
\varepsilon_x = \dfrac{1}{E}(\sigma_x-\nu(\sigma_y + \nu(\sigma_x + \sigma_y)))
$$
$$
\varepsilon_y = \dfrac{1}{E}(\sigma_y-\nu(\sigma_x + \nu(\sigma_x + \sigma_y)))
$$

Con el siguiente código de Python se hace la factorización de estas expresiones:
	*(Se requiere [[Librerías para el Estudio de la Elasticidad Plana en Python y Definición de Variables Simbólicas|Librerías y Definición de Variables Simbólicas]])*
```python
#Resolver y factorizar para EpsilonX
sol1 = sp.factor(sp.solve(sp.Eq(ex, (1/E)*(sx - nu*(sy + nu*(sx + sy)))), ex))

#Resolver y factorizar para EpsilonY
sol2 = sp.factor(sp.solve(sp.Eq(ey, (1/E)*(sy - nu*(sx + nu*(sx + sy)))), ey))

#Se presentan las soluciones: EpsilonX y EpsilonY respectivamente
sol1, sol2
```

Donde todas las expresiones resultan:
$$
\begin{array}{c c}
\begin{aligned}
\varepsilon_x &= \frac{1+\nu}{E}\left((1-\nu)\sigma_x - \nu\sigma_y\right) \\[8pt]
\varepsilon_y &= \frac{1+\nu}{E}\left((1-\nu)\sigma_y - \nu\sigma_x\right) \\[8pt]
\varepsilon_z &= 0
\end{aligned}
\qquad
&
\begin{aligned}
\gamma_{xy} &= \frac{1}{G}\tau_{xy} \\[8pt]
\gamma_{yz} &= 0 \\[8pt]
\gamma_{xz} &= 0
\end{aligned}
\end{array}
$$

Para despejar los esfuerzos, se hace uso de este código de Python:
```python
#Definir las ecuaciones
eq1 = sp.Eq(ex, (1/E)*(sx - nu*(sy + nu*(sx + sy))))
eq2 = sp.Eq(ey, (1/E)*(sy - nu*(sx + nu*(sx + sy))))

#Resolver el sistema respecto a SigmaX y SigmaY
sol = sp.factor(sp.solve([eq1, eq2], [sx, sy]))

# Mostrar resultados simplificados
sx_sol = sp.simplify(sol[sx])
sy_sol = sp.simplify(sol[sy])

#Reemplazo y simplificación en la ecución de SigmaZ
solz = sp.factor(sp.simplify(nu*(sx_sol + sy_sol)))

#Se presentan las soluciones: SigmaX, SigmaY y SigmaZ respectivamente
sol, solz
```

Resultando en:
$$
\begin{array}{c c}
\begin{aligned}
\sigma_x &= \frac{E}{(1+\nu)(1-2\nu)}
\left[(1-\nu)\varepsilon_x + \nu\varepsilon_y\right] \\[12pt]

\sigma_y &= \frac{E}{(1+\nu)(1-2\nu)}
\left[(1-\nu)\varepsilon_y + \nu\varepsilon_x\right] \\[12pt]

\sigma_z &= \frac{E\nu}{(1+\nu)(1-2\nu)}
\left(\varepsilon_x + \varepsilon_y\right)
\end{aligned}
&
\begin{aligned}
\tau_{xy} &= G\,\gamma_{xy} \\[12pt]

\tau_{yz} &= 0 \\[12pt]

\tau_{xz} &= 0
\end{aligned}
\end{array}
$$

Recordando que los esfuerzos se relacionan proporcionalmente a las deformaciones por medio del tensor de flexibilidad:
$$
\epsilon_{ij} = S_{ijkl}\,\sigma_{kl}
$$

Se puede organizar matricialmente la información obtenida para establecer la relación Deformación - Esfuerzo:
$$
\begin{pmatrix}
\sigma_x \\
\sigma_y \\
\tau_{xy}
\end{pmatrix}
=
\frac{E}{(1+\nu)(1-2\nu)}
\begin{pmatrix}
1-\nu & \nu & 0 \\
\nu & 1-\nu & 0 \\
0 & 0 & \dfrac{1-2\nu}{2}
\end{pmatrix}
\begin{pmatrix}
\varepsilon_x \\
\varepsilon_y \\
\gamma_{xy}
\end{pmatrix}
$$
$$
\underline{\underline{\varepsilon}}=\underline{\underline{D_{DP}}}\;\underline{\underline{\sigma}}
$$

Siendo $\underline{\underline{D_{DP}}}$ la matriz constitutiva de la relación Deformación - Esfuerzo en el estado de deformación plana.
En general las matrices resultantes son:
$$
\begin{aligned}
\boldsymbol{\sigma} &=
\begin{bmatrix}
\sigma_x & \tau_{xy} & 0 \\
\tau_{xy} & \sigma_y & 0 \\
0 & 0 & \sigma_z
\end{bmatrix}
\qquad
&
\boldsymbol{\varepsilon} &=
\begin{bmatrix}
\varepsilon_x & \gamma_{xy} & 0 \\
\gamma_{xy} & \varepsilon_y & 0 \\
0 & 0 & 0
\end{bmatrix}
\end{aligned}
$$

Se observa que, a pesar de no existir deformación en el eje z, si existe un esfuerzo en este eje debido a que la acción de los esfuerzos axiales de los otros ejes generan deformación en el eje z para compensarse volumétricamente. La acción del esfuerzo $\sigma_z$ contrarresta esta deformación y por ello depende de el coeficiente de Poisson y los otros dos esfuerzos axiales.

Ya que, matemáticamente, el tensor de deformaciones se comporta de la misma manera que el tensor de esfuerzos, existe una gran similitud con la [[Tensión Plana En Materiales Isótropos]].
