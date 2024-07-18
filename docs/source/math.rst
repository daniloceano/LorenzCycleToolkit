Mathematics
============

The energy budget equations are as follows:

.. math::
   :nowrap:

   \begin{align*}
   \frac{\partial A_Z}{\partial t} &= C_K - C_A + BA_Z + \Delta G_Z \\
   \frac{\partial A_E}{\partial t} &= C_A - C_E + BA_E + \Delta G_E \\
   \frac{\partial K_Z}{\partial t} &= C_K - C_Z + BK_Z - \Delta R_Z \\
   \frac{\partial K_E}{\partial t} &= C_E - C_K + BK_E - \Delta R_E
   \end{align*}

In these equations, available potential energy (APE) is divided into
zonal (:math:`A_Z`) and eddy (:math:`A_E`) components, as is kinetic
energy (:math:`K_Z` and :math:`K_E`, respectively). The transformations
between these forms of energy are denoted by :math:`C`, with subscripts
:math:`Z` and :math:`E` for conversions between zonal and eddy forms,
and :math:`A` and :math:`K` indicating conversions between APE and
kinetic energy, respectively. Thus, :math:`C_A` represents the
conversion between :math:`A_Z` and :math:`A_E`, :math:`C_E` denotes the
conversion from :math:`A_E` to :math:`K_E`, :math:`C_K` signifies the
transformation from :math:`K_E` to :math:`K_Z`, and :math:`C_Z`
describes the conversion from :math:`A_Z` to :math:`K_Z`. The residual
terms are defined as:

.. math::
   :nowrap:

   \begin{align*}
   \Delta R_Z &= B \Phi_Z - D_Z + \epsilon_{KZ} \\
   \Delta R_E &= B \Phi_E - D_E + \epsilon_{KE} \\
   \Delta G_Z &= G_Z + \epsilon_{GZ} \\
   \Delta G_E &= G_E + \epsilon_{GE}
   \end{align*}

Where APE generation and dissipation of kinetic energy are indicated by
:math:`G` and :math:`D`, with :math:`G_Z` and :math:`G_E` marking the
generation of :math:`A_Z` and :math:`A_E`, and :math:`D_Z` and
:math:`D_E` representing the dissipation of :math:`K_Z` and :math:`K_E`,
respectively.

Firstly, we define the zonal mean of a variable :math:`X`, between
longitudes :math:`\lambda_{1}` and :math:`\lambda_{2}`:

.. math::

   \begin{align*}
   [X]_\lambda &= \frac{1}{\lambda_2 - \lambda_1} \int_{\lambda_2}^{\lambda_1} X d\lambda
   \end{align*}

The eddy component of this variable is its deviation from the zonal mean:

.. math::

   \begin{align*}
   (X)_\lambda &=  X - [X]_\lambda
   \end{align*}

The domain mean of the variable :math:`X`, defined over the computational domain bounded by longitudes :math:`\lambda_1` and :math:`\lambda_2`, and latitudes :math:`\varphi_1` and :math:`\varphi_2`, is given by:

.. math::

   \begin{align*}
   [X]_{\lambda\phi} &= \left(\frac{1}{\lambda_2 - \lambda_1}\right)  \left(\frac{1}{\sin\phi_2 - \sin\phi_1}\right) \int_{\lambda_2}^{\lambda_1} X \cos\phi d\lambda d\phi 
   \end{align*}

Similarly, we define the deviation of the zonal mean from the domain mean:

.. math::

   \begin{align*}
   ([X]_\lambda)_\phi &= [X]_\lambda - [X]_{\lambda\phi}
   \end{align*}

From the definitions above, the four energy components used in the LEC computation are defined as follows:

.. math::
   \begin{align*}
   A_Z &= \int_{p_t}^{p_b} \frac{([(T)_\lambda ])_\phi^{2}]_{\lambda \phi}}  {2[\sigma]_{\lambda \phi}} dp \\
   A_E &= \int_{p_t}^{p_b} \frac{[(T)_\lambda^{2}]_{\lambda \phi}]}  {2[\sigma]_{\lambda \phi}} dp \\
   K_Z &=  \int_{p_t}^{p_b} \frac{[[u]_\lambda^2 + [v]_\lambda^2]_{\lambda \phi}}{2g} dp \\
   K_E &=  \int_{p_t}^{p_b} \frac{[(u)_\lambda^2 + (v)_\lambda^2]_{\lambda \phi}}{2g} dp
   \end{align*}

where :math:`p` is the atmospheric pressure, with subscripts :math:`b` and :math:`t` denoting the lower (base) and upper (top) pressure boundaries of the atmosphere, respectively. :math:`T` represents temperature, :math:`g` is the acceleration due to gravity, and :math:`u` and :math:`v` are the zonal and meridional wind components, respectively. The static stability parameter :math:`\sigma` is defined as:

.. math::
   \begin{align*}
   \sigma &= \left[\frac{gT}{c_p}-\frac{pg}{R}\frac{\partial T}{\partial p}\right]_{\lambda \phi}
   \end{align*}

where :math:`c_p` is the specific heat at constant pressure, and :math:`R` is the ideal gas constant for dry air.

The four conversion terms are defined as follows, integrating over the atmospheric column from the base (:math:`p_b`) to the top (:math:`p_t`) pressures:

.. math::

   \begin{aligned}
       &C_Z = \int_{p_t}^{p_b} - [\left([T]_\lambda)_\phi ([\omega]_\lambda\right)_\phi]_{\lambda\phi} \ \frac{R}{gp} \ dp \label{eq:CZ} \\
       &C_E = \int_{p_t}^{p_b} - [(T)_\lambda (\omega)_\lambda]_{\lambda\phi} \ \frac{R}{gp} \ dp \label{eq:CE} \\
       &C_A = \int_{p_t}^{p_b} - \left( \frac{1}{2a\sigma}  \left[ (v)_\lambda (T)_\lambda  \frac{\partial  ([T]_\lambda)_\phi}{\partial \phi} \right]_{\lambda\phi} + \frac{1}{\sigma}  \left[ (\omega)_\lambda (T)_\lambda \frac{\partial  ([T]_\lambda)_\phi}{\partial p} \right]_{\lambda\phi} \right) dp \label{eq:CA} \\
       &C_K = \int_{p_t}^{p_b} \frac{1}{g} \left(\left[ \frac{\cos\phi}{a} (u)_\lambda (v)_\lambda \frac{\partial}{\partial\phi} \left(\frac{[u]_\lambda}{\cos\phi}\right)\right]_{\lambda\phi} + \left[ \frac{(v)_\lambda^2}{a} \frac{\partial [v]_\lambda}{\partial\phi}  \right]_{\lambda\phi}  \right.
       \left. + \left[ \frac{\tan\phi}{a} (u)_\lambda^2 [v]_\lambda  \right]_{\lambda\phi} + \left[ (\omega)_\lambda  (u)_\lambda \frac{\partial [u]_\lambda}{\partial p} \right]_{\lambda\phi} + \left[ (\omega)_\lambda  (v)_\lambda \frac{\partial [v]_\lambda}{\partial p} \right]_{\lambda\phi}  \right) dp \label{eq:CK}
   \end{aligned}

where :math:`a` is the Earth's radius and :math:`\omega` is the vertical velocity in isobaric coordinates.

The APE generation and K dissipation terms are defined as:

.. math::

   \begin{align*}
   G_Z &=  \int_{p_t}^{p_b} \frac{[([q]_\lambda)_\phi ([T]_\lambda)_\phi]_{\lambda \phi}}{c_p[\sigma]_{\lambda \phi}} dp \\
   G_E &=  \int_{p_t}^{p_b} \frac{[(q)_\lambda (T)_\lambda]_{\lambda \phi}}{c_p[\sigma]_{\lambda \phi}} dp \\
   D_Z &= -  \int_{p_t}^{p_b} \frac{1}{g} [[u]_\lambda [F_\lambda]_\lambda + [v]_\lambda [F_\phi]_\lambda]_{\lambda \phi} dp \\
   D_E &= -  \int_{p_t}^{p_b} \frac{1}{g} [(u)_\lambda (F_\lambda)_\lambda + (v)_\lambda (F_\phi)_\lambda]_{\lambda \phi} dp
   \end{align*}

Here, :math:`F_{\lambda}` and :math:`F_{\varphi}` represent the zonal and meridional frictional components, respectively, and :math:`q` is the diabatic heating term, computed as a residual from the thermodynamic equation:

.. math::
   \begin{align*}
   \frac{q}{c_p} &= \frac{\partial T}{\partial t} + \vec{V}_H \cdot \vec{\nabla}_p T - S_p\omega
   \end{align*}

where :math:`\vec{V}_H \cdot \vec{\nabla}_p T` represents the horizontal advection of temperature and :math:`S_p` approximates the static stability, given by:

.. math::
   \begin{align*}
   S_p &\equiv -\frac{T}{\theta}\frac{\partial \theta}{\partial p}
   \end{align*}

where :math:`\theta` is the potential temperature.

The boundary terms are given by:

.. math::

   \begin{aligned}
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
   \end{aligned}

where :math:`c_1=-\left[a\left(\lambda_2-\lambda_1\right)\left(\sin \varphi_2-\sin \varphi_1\right)\right]^{-1}, c_2=-\left[a\left(\sin \varphi_2-\sin \varphi_1\right)\right]^{-1}`.

Lastly, the terms :math:`B\Phi_Z` and :math:`B\Phi_E` are given by:

.. math::

   \begin{aligned}
   \mathrm{B} \Phi \mathrm{Z}= & c_1 \int_{p_1}^{p_2} \int_{\varphi_1}^{\varphi_2} \frac{1}{g}\left([v]_\lambda\left([\Phi]_\lambda\right)_{\varphi}\right)_{\lambda_1}^{\lambda_2} d \varphi d p \nonumber \\
   & +c_2 \int_{p_1}^{p_2} \frac{1}{g}\left(\cos \varphi[v]_\lambda\left([\Phi]_\lambda\right)_{\varphi}\right)_{\varphi_1}^{\varphi_2} d p  \\
   & -\frac{1} {g}\left(\left[\left([\omega]_\lambda\right)_{\varphi}\left([\Phi]_\lambda\right)_{\varphi}\right]_{\lambda_{\varphi}}\right)_{p_1}^{p_2} \nonumber \\
   \mathrm{~B} \Phi \mathrm{E}= & c_1 \int_{p_1}^{p_2} \int_{\varphi_1}^{\varphi_2} \frac{1}{g}\left((u)_\lambda(\Phi)_{\lambda_\lambda}\right)_{\lambda_1}^{\lambda_2} d \varphi d p \nonumber \\
   & +c_2 \int_{p_1}^{p_2} \frac{1}{g}\left(\left[(v)_\lambda(\Phi)_{\lambda_\lambda}\right]_\lambda \cos \varphi\right)_{\varphi_1}^{\varphi_2} d p \\
   & -\frac{1}{g}\left(\left[(\omega)_\lambda(\Phi)_\lambda\right]_{\lambda_{\varphi}}\right)_{p_1}^{p_2} \nonumber
   \end{aligned}