---
title: Tensión Plana en Materiales Isótropos
---
La tensión plana es un estado de esfuerzos específico cuando se anula (o se desprecia) solo uno de los esfuerzos principales de un sólido. Por lo tanto ocurre en elementos con una dimensión muy pequeña respecto a las demás (placas, muros grandes, etc.).
Esta tensión plana sucede cuando:
- No existen cargas aplicadas en el eje z.
- Las cargas aplicadas son perpendiculares al eje z.
- Las cargas están distribuidas uniformemente en el espesor del elemento.

<div align="center">  
<img src="PlaneStress.png" width="450">  
</div>  
  
<p align="center" style="font-size: 0.85rem; opacity: 0.7;">
  <b>Figura 1.</b> Elemento Sólido en Estado de Tensión Plana.<br>
  Fuente: Álvarez Marín, D. A. (2023). 
  <i>
    <a href="https://bffrepositorio.unal.edu.co/server/api/core/bitstreams/b07bc9a3-7305-4e52-97f3-108813d014a2/content" target="_blank">
      Teoría de la elasticidad usando Matlab y Maxima. Volumen 1: fundamentos.
    </a>
  </i> Universidad Nacional de Colombia.
</p>

Matemáticamente se representa como la nulidad de esfuerzos en la dirección del eje z:
$$
\sigma_z = \tau_{xz} = \tau_{yz} = 0
$$
De esta manera, el tensor de esfuerzos en estado de tensión plana es:
$$
\underline{\underline{\sigma}} =
\begin{bmatrix}
\sigma_x & \tau_{xy} & 0 \\
\tau_{xy} & \sigma_y & 0 \\
0 & 0 & 0
\end{bmatrix}
$$

<span style="opacity:0.7; font-size:0.85rem;"><em>Nótese que el rango del tensor disminuye ya que tiene una columna sin aporte dimensional. Esto hace que el espacio donde trabaje sea máximo &#8477&#178 (cuando &sigma;<sub>x</sub>&sigma;<sub>y</sub> &minus; &tau;<sub>xy</sub><sup>2</sup> &ne; 0), &#8477&#185 (cuando &sigma;<sub>x</sub>&sigma;<sub>y</sub> &minus; &tau;<sub>xy</sub><sup>2</sup> = 0) o &#8477&#8304 (cuando todos los esfuerzos son nulos). Pero en general, se está llevando el problema tridimensional a uno bidimensional para figuras planas.
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

Con las condiciones de tensión plana y recordando que el módulo de elasticidad transversal $G$ se define como $G = \frac{E}{2(1+\nu)}$, se obtiene:
$$
\begin{alignedat}{2}

\varepsilon_x &= \frac{1}{E}\left(\sigma_x - \nu\sigma_y\right)
\qquad\qquad
&\gamma_{xy} &= G\,\tau_{xy} \\[12pt]

\varepsilon_y &= \frac{1}{E}\left(\sigma_y - \nu\sigma_x\right)
\qquad\qquad
&\gamma_{yz} &= 0 \\[12pt]

\varepsilon_z &= -\frac{\nu}{E}\left(\sigma_x + \sigma_y\right)
\qquad\qquad
&\gamma_{xz} &= 0

\end{alignedat}
$$

Con el siguiente código de Python se despejan los esfuerzos resultantes:
	*(Se requiere [[Librerías para el Estudio de la Elasticidad Plana en Python y Definición de Variables Simbólicas|Librerías y Definición de Variables Simbólicas]])*
```python
#Despejar SigmaX de la expresión de EpsilonX
sol1 = sp.solve(sp.Eq(ex, (1/E)*(sx - nu*sy)), sx)

#Despejar SigmaY de la expresión de EpsilonY
sol2 = sp.solve(sp.Eq(ey, (1/E)*(sy + nu*sx)), sy)

#Se conoce que SigmaZ = 0 por estado de esfuerzos en tensión plana
sol3 = 0

#Se presentan las soluciones: SigmaX, SigmaY, SigmaZ respectivamente
sol1, sol2, sol3
```

Resultando en:
$$
\begin{aligned}
\sigma_x &= \frac{E}{1-\nu^2}\left(\varepsilon_x + \nu \varepsilon_y \right)
&\qquad
\tau_{xy} &= G\,\gamma_{xy} \\[10pt]

\sigma_y &= \frac{E}{1-\nu^2}\left(\varepsilon_y + \nu \varepsilon_x \right)
&\qquad
\tau_{yz} &= 0 \\[10pt]

\sigma_z &= 0
&\qquad
\tau_{xz} &= 0
\end{aligned}
$$

Recordando que los esfuerzos se relacionan proporcionalmente a las deformaciones por medio del tensor de rigidez:
$$
\sigma_{ij} = C_{ijkl}\,\varepsilon_{kl}
$$

Se puede organizar matricialmente la información anteriormente obtenida para establecer la relación Esfuerzo - Deformación:
$$
\begin{pmatrix}
\sigma_x \\
\sigma_y \\
\tau_{xy}
\end{pmatrix}
=
\begin{pmatrix}
\dfrac{E}{1-\nu^2} & \dfrac{E\nu}{1-\nu^2} & 0 \\
\dfrac{E\nu}{1-\nu^2} & \dfrac{E}{1-\nu^2} & 0 \\
0 & 0 & \dfrac{E}{2(1+\nu)}
\end{pmatrix}
\begin{pmatrix}
\varepsilon_x \\
\varepsilon_y \\
\gamma_{xy}
\end{pmatrix}
$$
$$
\underline{\underline{\sigma}}=\underline{\underline{D_{TP}}}\;\underline{\underline{\varepsilon}}
$$

Siendo $\underline{\underline{D_{TP}}}$ la matriz constitutiva de la relación Esfuerzo - Deformación en el estado de tensión plana. 
En general, las matrices resultantes son:
$$
\begin{aligned}
\boldsymbol{\sigma} &=
\begin{bmatrix}
\sigma_x & \tau_{xy} & 0 \\
\tau_{xy} & \sigma_y & 0 \\
0 & 0 & 0
\end{bmatrix}
\qquad
&
\boldsymbol{\varepsilon} &=
\begin{bmatrix}
\varepsilon_x & \gamma_{xy} & 0 \\
\gamma_{xy} & \varepsilon_y & 0 \\
0 & 0 & \varepsilon_z
\end{bmatrix}
\end{aligned}
$$

Se observa que, a pesar de no existir esfuerzo en el eje z, si existe una deformación en dicho eje debida a la compensación volumétrica del sólido de manera perpendicular al esfuerzo axial aplicado.

<div align="center">  
<img src="PoissonRatio.png" width="450">  
</div>
<p align="center" style="font-size: 0.85rem; opacity: 0.7;">
  <b>Figura 2.</b> Diagrama de un Elemento de Material Isotrópico Elástico Lineal Deformado en los Tres Ejes por Acción de Esfuerzos en solo uno de los Ejes.<br>
  Fuente: Wikipedia Contributors. (2019). <i>Wikipedia.</i>
  <i>
    <a href="https://es.wikipedia.org/wiki/Coeficiente_de_Poisson" target="_blank">
      Coeficiente de Poisson.
    </a>
  </i>
</p>
Ya que, matemáticamente, el tensor de esfuerzos se comporta de la misma manera que el tensor de deformaciones, existe una gran similitud con la [[Deformación Plana en Materiales Isótropos]].
