---
title: 'LorenzCycleToolkit: A Tool for Calculating the Lorenz Energy Cycle'
tags:
  - Python
  - meteorology
  - atmospheric dynamics
  - diagnostic
  - cyclones
authors:
  - name: Danilo Couto de Souza
    orcid: 0000-0003-4121-7583
    affiliation: 1
  - name: Pedro Leite da Silva Dias
    orcid: 0000-0002-4051-2962
    affiliation: 1
  - name: Matheus Bonjour Laviola da Silva
    orcid: 0000-0002-8629-8762
    affiliation: 1
    affiliation: 2
  - name: Ricardo de Camargo
    affiliation: 1
affiliations:
  - name: Institute of Astronomy, Geophysics and Atmospheric Sciences of the São Paulo University, Rua do Matão, 226, Cidade Universitária, 05508-090, São Paulo, Brazil
    index: 1
  - name: OceanPact, Rua da Glória, 122, Rio de Janeiro, 20241-180, Rio de Janeiro, Brazil
    index: 2
date: 2024-07-16
bibliography: paper.bib
---

# Summary

The LorenzCycleToolkit is a Python package for calculating the Lorenz Energy Cycle (LEC) in limited atmospheric regions. Introduced by Edward Lorenz in 1965, the LEC estimates atmospheric energy, including zonal and eddy components of Kinetic and Available Potential Energy, and conversions between these forms. This toolkit allows researchers to analyze energy transformations in various atmospheric phenomena, enhancing our understanding of atmospheric dynamics and aiding in weather and climate prediction.

# Background

The Lorenz Energy Cycle (LEC), introduced by [@lorenz1955], is a framework to describe the processes of energy conversion and transfer within the Earth's atmosphere. It serves as a diagnostic tool for the general circulation of the atmosphere, providing insights into its dynamic and thermodynamic structures, and information about atmospheric flows and the transport of heat, moisture, and angular momentum. The LEC divides atmospheric energy into four main reservoirs: available potential energy (APE), kinetic energy (K), and their respective zonal mean (Z) and eddy (E) components. The LEC can indicate how energy is generated, converted, and dissipated within these reservoirs.

Applications of the LEC are diverse and significant in meteorology and climate science. It aids in diagnosing the performance of numerical weather prediction models by comparing simulated and observed energy cycles [@ulbrich1991global]. Moreover, it is used to study climate change impacts, such as shifts in energy transfer processes due to global warming [@li2007lorenz; @veiga2013global]. Understanding the LEC is also crucial for improving the representation of atmospheric processes in general circulation models (GCMs), which are essential for accurate climate projections [@lorenz1967nature; @boer2008energy].

## Statement of Need

`LorenzCycleToolkit` is a Python package designed for calculating the Lorenz Energy Cycle (LEC) in limited atmospheric regions, supporting both Eulerian and Semi-Lagrangian frameworks. The Eulerian framework computes the LEC for a fixed region in the atmosphere [@dias2011energy], requiring users to provide a "box limits" file containing the minimum and maximum latitudes and longitudes for the desired computational domain. The Semi-Lagrangian framework, on the other hand, uses a domain that changes position with time [@michaelides1999quasi], which is especially useful for studying cyclonic systems. For this option, users must specify the central position of the domain for each timestep in a "track" CSV file, which can also specify the domain width and length. If these specifications are not provided, a default 15°x15° domain is used. The `LorenzCycleToolkit` also offers an option for users to interactively choose the domain position and size for each data timestep.

The program can handle any reanalysis or model data provided in isobaric levels. Users must adjust the "namelist" file, specifying the names for each variable and coordinate in the input file. The requirements for the input file are that it must be in NetCDF format [@rew1990netcdf], the grid must be structured, and it must have isobaric levels as vertical coordinates. These specifications are common in reanalysis, climate, and weather datasets.

The program outputs the results for each term in comma-separated values (CSV) files, generating one CSV file for each term. For energy, conversion, and generation terms, it also provides results for each vertical level, enabling two-dimensional analysis. Additionally, the program can be configured by the user to automatically generate plots. These plots include spatial representations of the computational domain defined for the analysis, time series of the system's central minimum pressure and relative vorticity at 850 hPa, time series of each LEC term, Hovmöller diagrams, and complete energy flux boxes for the LEC.


# External Libraries Used

For handling meteorological constants and variable units, the program uses the MetPy package [@may2022metpy]. Dask is used to optimize the handling of large datasets [@daniel2019data], allowing the slicing of large data chunks more efficiently and helping to manage memory usage and improve performance when dealing with extensive data in the NetCDF file.


# Mathematics

The energy budget equations are as follows:

$$
\begin{aligned}
\frac{\partial A_Z}{\partial t} &= C_K - C_A + BA_Z + \Delta G_Z \\
\frac{\partial A_E}{\partial t} &= C_A - C_E + BA_E + \Delta G_E \\
\frac{\partial K_Z}{\partial t} &= C_K - C_Z + BK_Z - \Delta R_Z \\
\frac{\partial K_E}{\partial t} &= C_E - C_K + BK_E - \Delta R_E 
\end{aligned}
$$


In these equations, APE is divided into zonal ($A_Z$) and eddy ($A_E$) components, as is K ($K_Z$ and $K_E$, respectively). The transformations between these forms of energy are denoted by $C$, with subscripts $Z$ and $E$ for conversions between zonal and eddy forms, and $A$ and $K$ indicating conversions between APE and kinetic energy, respectively. Thus, $C_A$ represents the conversion between $A_Z$ and $A_E$, $C_E$ denotes the conversion from $A_E$ to $K_E$, $C_K$ signifies the transformation from $K_E$ to $K_Z$, and $C_Z$ describes the conversion from $A_Z$ to $K_Z$. The residual terms are defined as:

$$
\begin{align*}
\Delta R_Z &= B\Phi_Z - D_Z + \epsilon_{KZ} \\
\Delta R_E &= B\Phi_E - D_E + \epsilon_{KE} \\
\Delta G_Z &= G_Z + \epsilon_{GZ} \\
\Delta G_E &= G_E + \epsilon_{GE}
\end{align*}
$$

Where APE generation and dissipation of kinetic energy are indicated by $G$ and $D$, with $G_Z$ and $G_E$ marking the generation of $A_Z$ and $A_E$, and $D_Z$ and $D_E$ representing the dissipation of $K_Z$ and $K_E$, respectively.

Firstly, we define the zonal mean of a variable $X$, between longitudes $\lambda_{1}$ and $\lambda_{2}$:

$$
[X]_\lambda = \frac{1}{\lambda_2 - \lambda_1} \int_{\lambda_2}^{\lambda_1} X d\lambda
$$

The eddy component of this variable is its deviation from the zonal mean:

$$
(X)_\lambda =  X - [X]_\lambda
$$

The domain mean of the variable $X$, defined over the computational domain bounded by longitudes $\lambda_1$ and $\lambda_2$, and latitudes $\phi_1$ and $\phi_2$, is given by:

$$
[X]_{\lambda\phi} = \left(\frac{1}{\lambda_2 - \lambda_1}\right)  \left(\frac{1}{\sin\phi_2 - \sin\phi_1}\right) \int_{\lambda_2}^{\lambda_1} X \cos\phi \ d\lambda \ d\phi 
$$

Similarly, we define the deviation of the zonal mean from the domain mean:

$$
([X]_\lambda)_\phi = [X]_\lambda - [X]_{\lambda\phi}
$$

From the definitions above, the four energy components used in the LEC computation are defined as follows:

$$
\begin{align*}
A_Z &= \int_{p_t}^{p_b} \frac{([(T)_\lambda ])_\phi^{2}]_{\lambda \phi}}  {2[\sigma]_{\lambda \phi}} dp \\
A_E &= \int_{p_t}^{p_b} \frac{[(T)_\lambda^{2}]_{\lambda \phi}]}  {2[\sigma]_{\lambda \phi}} dp \\
K_Z &= \int_{p_t}^{p_b} \frac{[[u]_\lambda^2 + [v]_\lambda^2]_{\lambda \phi}}{2g} dp \\
K_E &= \int_{p_t}^{p_b} \frac{[(u)_\lambda^2 + (v)_\lambda^2]_{\lambda \phi}}{2g} dp
\end{align*}
$$

where $p$ is the atmospheric pressure, with subscripts $b$ and $t$ denoting the lower (base) and upper (top) pressure boundaries of the atmosphere, respectively. $T$ represents temperature, $g$ is the acceleration due to gravity, and $u$ and $v$ are the zonal and meridional wind components, respectively. The static stability parameter $\sigma$ is defined as:

$$
\sigma = \left[\frac{gT}{c_p}-\frac{pg}{R}\frac{\partial T}{\partial p}\right]_{\lambda \phi}
$$

where $c_p$ is the specific heat at constant pressure, and $R$ is the ideal gas constant for dry air.

The four conversion terms are defined as follows, integrating over the atmospheric column from the base ($p_b$) to the top ($p_t$) pressures:

$$
\begin{align*}
C_Z &= \int_{p_t}^{p_b} - [\left([T]_\lambda)_\phi ([\omega]_\lambda\right)_\phi]_{\lambda\phi} \ \frac{R}{gp} \ dp \\
C_E &= \int_{p_t}^{p_b} - [(T)_\lambda (\omega)_\lambda]_{\lambda\phi} \ \frac{R}{gp} \ dp \\
C_A &= \int_{p_t}^{p_b} - \left( \frac{1}{2a\sigma}  \left[ (v)_\lambda (T)_\lambda  \frac{\partial  ([T]_\lambda)_\phi}{\partial \phi} \right]_{\lambda\phi} + \frac{1}{\sigma}  \left[ (\omega)_\lambda (T)_\lambda \frac{\partial  ([T]_\lambda)_\phi}{\partial p} \right]_{\lambda\phi} \right) dp \\
C_K &= \int_{p_t}^{p_b} \frac{1}{g} \left( \left[ \frac{\cos\phi}{a} (u)_\lambda (v)_\lambda \frac{\partial}{\partial\phi} \left(\frac{[u]_\lambda}{\cos\phi}\right)\right]_{\lambda\phi} + \left[ \frac{(v)_\lambda^2}{a} \frac{\partial [v]_\lambda}{\partial\phi}  \right]_{\lambda\phi} \right. \\
&\left. + \left[ \frac{\tan\phi}{a} (u)_\lambda^2 [v]_\lambda  \right]_{\lambda\phi} + \left[ (\omega)_\lambda  (u)_\lambda \frac{\partial [u]_\lambda}{\partial p} \right]_{\lambda\phi} + \left[ (\omega)_\lambda  (v)_\lambda \frac{\partial [v]_\lambda}{\partial p} \right]_{\lambda\phi}  \right) dp 
\end{align*}
$$


where $a$ is the Earth's radius and $\omega$ is the vertical velocity in isobaric coordinates. 

The APE generation and K dissipation terms are defined as:

$$
\begin{align*}
G_Z &=  \int_{p_t}^{p_b} \frac{[([q]_\lambda)_\phi ([T]_\lambda)_\phi]_{\lambda \phi}}{c_p[\sigma]_{\lambda \phi}} dp \\
G_E &=  \int_{p_t}^{p_b} \frac{[(q)_\lambda (T)_\lambda]_{\lambda \phi}}{c_p[\sigma]_{\lambda \phi}} dp \\
D_Z &= -  \int_{p_t}^{p_b} \frac{1}{g} [[u]_\lambda [F_\lambda]_\lambda + [v]_\lambda [F_\phi]_\lambda]_{\lambda \phi} dp \\
D_E &= -  \int_{p_t}^{p_b} \frac{1}{g} [(u)_\lambda (F_\lambda)_\lambda + (v)_\lambda (F_\phi)_\lambda]_{\lambda \phi} dp
\end{align*}
$$

Here, $F_{\lambda}$ and $F_{\phi}$ represent the zonal and meridional frictional components, respectively, and $q$ is the diabatic heating term, computed as a residual from the thermodynamic equation:

$$
\frac{q}{c_p} = \frac{\partial T}{\partial t} + \vec{V}_H \cdot \vec{\nabla}_p T - S_p \omega
$$

where $\vec{V}_H \cdot \vec{\nabla}_p T$ represents the horizontal advection of temperature and $S_p$ approximates the static stability, given by:

$$
S_p \equiv -\frac{T}{\theta}\frac{\partial \theta}{\partial p}
$$

where $\theta$ is the potential temperature. 

The boundary terms are given by:

$$
\begin{align}
& \mathrm{BAZ}=c_1 \int_{p_1}^{p_2} \int_{\varphi_1}^{\varphi_2} \frac{1}{2[\sigma]_{\lambda_{\varphi}}}\left(2\left([T]_\lambda\right)_{\varphi}(T)_\lambda u+\left([T]_{\lambda_{\varphi}}\right)_{\varphi}^2 u\right)_{\lambda_1}^{\lambda_2} \nonumber \\
& \times d \varphi d p+c_2 \int_{p_1}^{p_2} \frac{1}{2[\sigma]_{\lambda \varphi}}\left(2\left[(v)_\lambda(T)_\lambda\right]_\lambda\left([T]_\lambda\right)_{\varphi} \cos \varphi \right. \left.+\left([T]_\lambda\right)_{\varphi}^2[v]_\lambda \cos \varphi\right)_{\varphi_1}^{\varphi_2} d p \nonumber \\
& -\frac{1}{2[\sigma]_{\lambda \varphi}}\left(\left[2(\omega)_\lambda(T)_\lambda\right]_\lambda\left([T]_\lambda\right)_{\varphi}+\left[[\omega]_\lambda\left([T]_\lambda\right)_{\varphi}^2\right]_{\lambda_{\varphi}}\right)_{p_1}^{p_2} \\
& \mathrm{BAE}=c_1 \int_{p_1}^{p_2} \int_{\varphi_1}^{\varphi_2} \frac{1}{2[\sigma]_{\lambda \varphi}}\left[u(T)_\lambda^2\right]_{\lambda_1}^{\lambda_2} d \varphi d p \nonumber \\
& +c_2 \int_{p_1}^{p_2} \frac{1}{2[\sigma]_{\lambda \varphi}}\left(\left[(T)_\lambda^2 v\right]_\lambda \cos \varphi\right)_{\varphi_1}^{^{\varphi_2}} d p \\
& -\left(\frac{\left[\omega(T)_\lambda^2\right]_{\lambda \varphi}}{2[\sigma]_{\lambda \varphi}}\right)_{p_1}^{p_2} \nonumber \\
& \mathrm{BKZ}=c_1 \int_{p_1}^{p_2} \int_{\varphi_1}^{\varphi_2} \frac{1}{2 g}\left(u\left[u^2+v^2-(u)_\lambda^2-(v)_\lambda^2\right]\right)_{\lambda_1}^{\lambda_2} \nonumber \\
& \times d \varphi d p+c_2 \int_{p_1}^{p_2} \frac{1}{2 g}\left(\left[v \cos \varphi \left[u^2+v^2\right.\right.\right. \left.\left.-(u)_\lambda^2-(v)_\lambda^2\right]\right]_{\varphi_1}^{\varphi_2} d p  \\
& -\left(\frac{1}{2 g}\left[\omega\left[u^2+v^2-(u)_\lambda^2-(v)_\lambda^2\right]\right]_{\lambda \varphi}\right)_{p_1}^{p_2} \nonumber \\
& \mathrm{BKE}=c_1 \int_{p_1}^{p_2} \int_{\varphi_1}^{\varphi_2} \frac{1}{2 g}\left(u\left[(u)_\lambda^2+(v)_\lambda^2\right]\right)_{\lambda_1}^{\lambda_2} d \varphi d p \nonumber \\
& +c_2 \int_{p_1}^{p_2} \frac{1}{2 g}\left(\left[v \cos \varphi\left[(u)_\lambda^2+(v)_\lambda^2\right]\right]_\lambda\right)_{\varphi_1}^{\varphi_2} d p \\
& -\left(\frac{1}{2 g}\left[\omega\left[(u)_\lambda^2+(v)_\lambda^2\right]\right]_{\lambda \varphi}\right)_{p_1}^{p_2} \nonumber
\end{align}
$$

where $c_1=-\left[a\left(\lambda_2-\lambda_1\right)\left(\sin \varphi_2-\sin \varphi_1\right)\right]^{-1}, c_2=-[a$ $\left.x\left(\sin \varphi_2-\sin \varphi_1\right)\right]^{-1}$.

Lastly, the terms $B\Phi Z$ and $B\Phi E$ are given by:

$$
\begin{align}
\mathrm{B} \Phi \mathrm{Z}= & c_1 \int_{p_1}^{p_2} \int_{\varphi_1}^{\varphi_2} \frac{1}{g}\left([v]_\lambda\left([\Phi]_\lambda\right)_{\varphi}\right)_{\lambda_1}^{\lambda_2} d \varphi d p \nonumber \\
& +c_2 \int_{p_1}^{p_2} \frac{1}{g}\left(\cos \varphi[v]_\lambda\left([\Phi]_\lambda\right)_{\varphi}\right)_{\varphi_1}^{\varphi_2} d p  \\
& -\frac{1} {g}\left(\left[\left([\omega]_\lambda\right)_{\varphi}\left([\Phi]_\lambda\right)_{\varphi}\right]_{\lambda_{\varphi}}\right)_{p_1}^{p_2} \nonumber \\
\mathrm{~B} \Phi \mathrm{E}= & c_1 \int_{p_1}^{p_2} \int_{\varphi_1}^{\varphi_2} \frac{1}{g}\left((u)_\lambda(\Phi)_{\lambda_\lambda}\right)_{\lambda_1}^{\lambda_2} d \varphi d p \nonumber \\
& +c_2 \int_{p_1}^{p_2} \frac{1}{g}\left(\left[(v)_\lambda(\Phi)_{\lambda_\lambda}\right]_\lambda \cos \varphi\right)_{\varphi_1}^{\varphi_2} d p \\
& -\frac{1}{g}\left(\left[(\omega)_\lambda(\Phi)_\lambda\right]_{\lambda_{\varphi}}\right)_{p_1}^{p_2} \nonumber
\end{align}
$$

# References