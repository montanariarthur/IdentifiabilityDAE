# Identifiability of Differential-Algebraic Equations
Codes for parameter identifiability analysis of models described by differential-algebraic equations (DAE).

# Usage

The following codes are direct implementations of the parameter identifiability analysis conducted for the DAE systems studied in Ref. [1].

- `main_idft_chemicalreactor` : Analysis of the identifiability of parameter `Tc` of a chemical reactor model for different types of measurement signals (sensors). See the tutorial section below for a discussion of this example.

- `main_idft_pendulum` : Analysis of the identifiability of different combinations of parameters of a pendulum equation.

- `main_idft_linearDAE.m`: Analysis of the identifiability of a linear DAE model.

- `DAEobsvmatrix` : Function that computes the observability matrix of a DAE model.


The above codes are dependent on the following functions:

- `DAE models` : Folder contains implementations of the DAEs for the systems investigated in Ref. [1], which are used for numerical integration.
- `sym2str` : Converts symbolic variables to string variables. 

# Tutorial example

This tutorial example demonstrates how identifiability analysis can be applied to a chemical reactor model. The tutorial highlights the analysis can guide sensor placement and also examines certain numerical aspects, including the dependence of identifiability on the local state and the conditioning of the identifiability matrix.

### Chemical reactor model

Consider an exothermic reactor with a single first-order reaction, $$\mathrm{A} \xrightarrow{r(t)} \mathrm{B},$$ where the generated heat is removed through an external cooling circuit. The reactor is modeled by the DAEs

$$
\begin{aligned}
\dot{x}_1
    &= k_1(c_0-x_1)-x_3, \\
\dot{x}_2
    &= k_1(T_0-x_2)+k_2x_3+k_3(T_c-x_2), \\
0
    &= x_3-k_5\exp\left(-\frac{k_4}{x_2}\right)x_1.
\end{aligned}
$$

The state vector is $${x} = [x_1 , x_2 , x_3]^{\mathsf T} \in \mathbb{R}^3,$$ where $x_1$ is the concentration of reactant A, $x_2$ is the reactor temperature, and$x_3$ is the reaction rate per unit volume (algebraic variable). The algebraic equation arises from a conservation law (Arrhenius relation) that imposes the dependence of raction rates on temperature. The model parameters are: the inlet concentration $c_0$ of reactant A, the inlet feed temperature $T_0$, the coolant temperature $T_c$, and other reaction and heat-transfer constants $k_i$. Although the algebraic constraint can be substituted into the differential equations to obtain an ODE model, the MATLAB solver `ode15s` for DAEs, used in the code, directly enforces the constraint during numerical integration. For the parameter values considered below, the reactor dynamics converge to a limit cycle.

### Constructing the identifiability matrix

Let the coolant temperature be the parameter sought to be identified (i.e., $\theta=T_c$). The accuracy of parameter estimation depends on the choice of observable, that is, which state variable $x_i$ is measured by a sensor. To determine the most suitable sensor data for parameter estimation, we now test which observables ensure parameter identifiability.
As a tutorial, we analyze the identifiability of the chemical reactor model step-by-step, considering the output signal $y(t) = x_1(t)$, $t\in[0,10]$ s. First, we augment the model by representing $\theta(t) = T_c(t)$ as an additional state variable with constant dynamics, i.e., $\dot\theta(t) = 0$. Functions $H_{\nu}$ and $F_{\mu}$ are then constructed up to order $\mu=\nu=n+p-1=3$. Thus, $H_{\nu} = [x_1 , \dot x_1 , \ddot x_1 ]^{\mathsf T}$ and $F_{\mu} = [(F)^{\mathsf T} , (dF/dt)^{\mathsf T} , (d^2F/dt^2)^{\mathsf T}]^{\mathsf T}$, where

$$
\begin{aligned}
F &= [k_1(c_0-x_1) - x_3 - \dot x_1 , \quad k_1(T_0-x_2) + k_2x_3 + k_3(\theta - x_2) - \dot x_2 , \quad x_3 - k_5\exp(-\frac{k_4}{x_2})x_1 , \quad \dot{\theta}]^{\mathsf T}, \\
\frac{dF}{dt} &= [-k_1\dot x_1 - \dot x_3 - \ddot x_1, \quad -k_1\dot x_2 + k_2 \dot x_3 + k_3(\dot\theta - \dot x_2) - \ddot x_2, \quad \dot x_3 - k_5{\exp}(-\frac{k_4}{x_2}) (\dot x_1 + k_4\frac{x_1\dot x_2}{x_2^2}), \quad 0]^{\mathsf T},
\end{aligned}
$$

and $d^i F/dt^i$, $i\geq 2$, is omitted for brevity. Note that $F_{\mu}$ and $H_{\nu}$ are functions of $\theta$, $x$, and its time derivatives. By assumption, $\dot\theta = \ddot\theta = \dddot\theta = 0$, $\forall t$. After building the identifiability matrix (8), condition (9) can then be locally evaluated to test if the system is identifiable at any consistent state $x_0\in\mathbb{L}\subset\mathbb{R}^3$ and parameter value $\theta_0\in\mathbb{R}^1$. 

### Local identifiability analysis and sensor placement problem

Fig. 1 shows the identifiable and unidentifiable regions in the state space depending on the observable choice. For the output $y(t)=x_1(t)$ (Fig. 1a), the limit cycle solution ${\Phi}_t(x(0);\theta)$ lies predominantly in an unidentifiable region. Following Definition 3, at an unidentifiable state $x_0\in\mathbb{L}$, there might exist a distinct parameter $\theta'\neq\theta$ within some neighborhood $\mathcal U$ of the true parameter value $\theta$ such that the corresponding measurement signals are indistinguishable: $h\circ{\Phi_t}(x_0;\theta) = h\circ{\Phi_t}(x_0;\theta')$. The prevalence of unidentifiable regions in the state space implies the possible existence of non-unique solutions for $\theta$ in the neighborhood of ${\Phi}_t(x_0;\theta)$, which can directly compromise the parameter estimation. That said, identifiability is still achieved in a small region of the limit cycle around $x = [0.25 , 360 , 0.5]^{\mathsf T}$ (as well as other regions containing the transient response). Theoretically, parameter estimation could be possible by using data collected solely at these identifiable states. However, this may not be possible in practice if the identifiable regions are too small (as in Fig. 1a) due to limited data availability.

After performing the same identifiability analysis for $y(t) = x_2(t)$ and $y(t) = x_3(t)$, Fig. 1b,c shows that the system is locally identifiable everywhere around the limit cycle solution when these other types of measurement signals (sensors) are considered. This result is particularly advantageous from a sensor placement perspective, as temperature sensors ($y = x_2$) are both affordable and practical to install in chemical reactors.

![](fig/fig_chemical.pdf)

**Fig. 1.**  Identifiable regions of the chemical reactor model for measurement signals given by: **(a)** $y = x_1$, **(b)** $y = x_2$, and **(c)** $y = x_3$. The blue (red) colors correspond to states $\bm x$ in which the parameter $T_c$ is identifiable (unidentifiable). The black solid line represents the state trajectory $x(t)$ starting at the initial condition $x(0) = [0.5 , 350 , 0.4995]^{\mathsf T}$ and converging to a limit cycle. The system parameters were set to $(c_0,T_0,T_c) = (1, 350,305)$, and $(k_1,k_2,k_3,k_4,k_5) = (1, 209.205,  2.0921, 8.7503\cdot 10^3,7.2\cdot 10^{10})$. The numerical simulations are shown using the `ode15s` solver in MATLAB.

### Numerical conditioning and sensitivity analysis

We note that the outcome of our identifiability analysis depends on the choice of tolerances used in the rank computations. Regions of identifiable/unidentifiable states may expand or shrink depending on the tolerance value, as illustrated in Fig. 2. If the rank tolerances are not appropriately set, condition (9) may be numerically satisfied for ill-conditioned identifiability matrices, even when the system is *practically* (i.e., almost) unidentifiable.

![](fig/fig_ranktolerance.pdf)

**Fig. 2.**  Identifiable regions of the chemical reactor for different choices of rank tolerances: **(a)** $10^{-1}$, **(b)** $10^{-6}$, **(c)** $10^{-15}$. The blue (red) colors correspond to states $x$ in which the parameter $T_c$ is identifiable (unidentifiable) using sensor $y = x_1$.


# License

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

The full text of the GNU General Public License can be found in the file "LICENSE.txt".

# References
1. A. N. Montanari, F. Lamoline, R. Bereza, J. Goncalves. Identifiability of differential-algebraic equations. *Under review*.
