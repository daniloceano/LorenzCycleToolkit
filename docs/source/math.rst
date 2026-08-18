Mathematics
============

The energy budget equations are as follows:

.. math::
   :nowrap:

   \begin{align*}
   \frac{\partial A_Z}{\partial t} &= -C_A - C_Z + BA_Z + \Delta G_Z \\
   \frac{\partial A_E}{\partial t} &= C_A - C_E + BA_E + \Delta G_E \\
   \frac{\partial K_Z}{\partial t} &= C_K + C_Z + BK_Z + \Delta R_Z \\
   \frac{\partial K_E}{\partial t} &= C_E - C_K + BK_E + \Delta R_E
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
   \Delta R_Z &= B \Phi_Z + D_Z + \epsilon_{KZ} \\
   \Delta R_E &= B \Phi_E + D_E + \epsilon_{KE} \\
   \Delta G_Z &= G_Z + \epsilon_{GZ} \\
   \Delta G_E &= G_E + \epsilon_{GE}
   \end{align*}

Where APE generation and dissipation of kinetic energy are indicated by
:math:`G` and :math:`D`, with :math:`G_Z` and :math:`G_E` marking the
generation of :math:`A_Z` and :math:`A_E`, and :math:`D_Z` and
:math:`D_E` representing the dissipation of :math:`K_Z` and :math:`K_E`,
respectively.

.. note::

   **Sign convention.** These budgets follow Brennan and Vincent (1980,
   Eqs. 5a-6b). All conversion terms are *signed rates*, and a positive value
   means transfer in the direction named by the term: positive :math:`C_A`
   transfers :math:`A_Z \rightarrow A_E`, positive :math:`C_Z` transfers
   :math:`A_Z \rightarrow K_Z`, positive :math:`C_E` transfers
   :math:`A_E \rightarrow K_E`, and positive :math:`C_K` transfers
   :math:`K_E \rightarrow K_Z`. :math:`D_Z` and :math:`D_E` are likewise
   *signed* frictional terms (negative for a loss), not positive-definite loss
   magnitudes; note that Michaelides (1987) writes them with the opposite sign
   convention, as :math:`-D` in the kinetic budgets.

   **Residuals.** :math:`\Delta R_Z` and :math:`\Delta R_E` are what the code
   exports as ``RKz`` and ``RKe``. They are *composite* residuals: each one
   lumps together the boundary pressure work, the frictional dissipation, any
   unresolved-scale transfer, and the accumulated numerical error. They must
   not be interpreted as friction alone. Likewise ``RGz`` and ``RGe``
   (:math:`\Delta G_Z`, :math:`\Delta G_E`) are generation plus numerical
   error. This matches the residual definitions of Dias Pinto and da Rocha
   (2011, their Eqs. 8-11).

   The code implements exactly these budgets in
   ``src/utils/calc_budget_and_residual.py``, which rearranges them to

   .. math::
      \begin{align*}
      RG_Z &= \frac{\partial A_Z}{\partial t} + C_Z + C_A - BA_Z \\
      RK_Z &= \frac{\partial K_Z}{\partial t} - C_Z - C_K - BK_Z \\
      RG_E &= \frac{\partial A_E}{\partial t} - C_A + C_E - BA_E \\
      RK_E &= \frac{\partial K_E}{\partial t} - C_E + C_K - BK_E
      \end{align*}

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
       &C_{\mathrm{overturning}} = -\int_{p_t}^{p_b} [\omega]_{\lambda\phi}\,\alpha\,\frac{dp}{g},
       \qquad \alpha=\frac{R[T]_{\lambda\phi}}{p} \label{eq:Coverlap} \\
       &C_A = \int_{p_t}^{p_b} - \left( \frac{1}{a\sigma}  \left[ (v)_\lambda (T)_\lambda  \frac{\partial  ([T]_\lambda)_\phi}{\partial \phi} \right]_{\lambda\phi} + \frac{1}{\sigma}  \left[ (\omega)_\lambda (T)_\lambda \frac{\partial  ([T]_\lambda)_\phi}{\partial p} \right]_{\lambda\phi} \right) dp \label{eq:CA} \\
       &C_K = \int_{p_t}^{p_b} \frac{1}{g} \left(\left[ \frac{\cos\phi}{a} (u)_\lambda (v)_\lambda \frac{\partial}{\partial\phi} \left(\frac{[u]_\lambda}{\cos\phi}\right)\right]_{\lambda\phi} + \left[ \frac{(v)_\lambda^2}{a} \frac{\partial [v]_\lambda}{\partial\phi}  \right]_{\lambda\phi}  \right.
       \left. + \left[ \frac{\tan\phi}{a} (u)_\lambda^2 [v]_\lambda  \right]_{\lambda\phi} + \left[ (\omega)_\lambda  (u)_\lambda \frac{\partial [u]_\lambda}{\partial p} \right]_{\lambda\phi} + \left[ (\omega)_\lambda  (v)_\lambda \frac{\partial [v]_\lambda}{\partial p} \right]_{\lambda\phi}  \right) dp \label{eq:CK}
   \end{aligned}

where :math:`a` is the Earth's radius and :math:`\omega` is the vertical velocity in isobaric coordinates.
The separately exported :math:`C_{\mathrm{overturning}}` diagnoses the
strength of domain-mean overturning.  It is **not** a missing conversion in
the Lorenz-cycle budget: in the exact pressure-work identity it cancels the
:math:`\overline{\Phi}\,\overline{\omega}` part of the top/bottom
geopotential flux.  It therefore remains outside :math:`RG_Z` and
:math:`RK_Z`.  Under the sign convention above, domain-mean ascent
(:math:`\omega<0`) produces a positive value.  It vanishes for the global
cycle, where mass continuity requires the horizontal mean
:math:`[\omega]_{\lambda\phi}=0`, but need not vanish in a limited area.  It
belongs neither to :math:`C_Z`, which uses the area anomaly
:math:`\omega^*`, nor to :math:`C_E`, which uses the zonal eddy
:math:`\omega'`.

Mass-continuity diagnostic
--------------------------

The toolkit also exports the pressure-coordinate mass-continuity residual

.. math::

   M(p)=\left\langle\frac{u|_{\lambda_e}-u|_{\lambda_w}}
   {a\cos\varphi\,\Delta\lambda}\right\rangle
   +\left\langle\frac{1}{a\cos\varphi}
   \frac{\partial([v]\cos\varphi)}{\partial\varphi}\right\rangle
   +\frac{\partial\overline{\omega}}{\partial p}.

The per-level file ``M_<pressure-coordinate>.csv`` contains :math:`M(p)` in
:math:`\mathrm{s^{-1}}`; the ``M`` column in the main result is
:math:`\int M\,dp/g` in :math:`\mathrm{kg\,m^{-2}\,s^{-1}}`.  The continuum
value is zero.  A nonzero value measures the mismatch between the archived
horizontal winds and archived :math:`\omega` and sets the numerical noise
floor for geopotential-flux diagnostics that use full :math:`\Phi`.

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
   \mathrm{B}\Phi_Z= &\ c_1 \int_{p_t}^{p_b}\!\!\int_{\varphi_s}^{\varphi_n}
   \frac{1}{g}\left([u]_\lambda\,\Delta_{EW}\Phi'
   + \Phi^*\,\Delta_{EW}u'\right) d\varphi\, dp \\
   &+c_2 \int_{p_t}^{p_b} \frac{1}{g}
   \left(\cos\varphi\,[v]_\lambda\,\Phi^*\right)_{\varphi_s}^{\varphi_n} dp
   -\frac{1}{g}\left(\left[\omega^*\Phi^*\right]_{\lambda\varphi}\right)_{p_t}^{p_b} \\
   \mathrm{~B} \Phi \mathrm{E}= & c_1 \int_{p_1}^{p_2} \int_{\varphi_1}^{\varphi_2} \frac{1}{g}\left((u)_\lambda(\Phi)_{\lambda_\lambda}\right)_{\lambda_1}^{\lambda_2} d \varphi d p \nonumber \\
   & +c_2 \int_{p_1}^{p_2} \frac{1}{g}\left(\left[(v)_\lambda(\Phi)_{\lambda_\lambda}\right]_\lambda \cos \varphi\right)_{\varphi_1}^{\varphi_2} d p \\
   & -\frac{1}{g}\left(\left[(\omega)_\lambda(\Phi)_\lambda\right]_{\lambda_{\varphi}}\right)_{p_1}^{p_2} \nonumber
   \end{aligned}

Here :math:`\Phi^*=[\Phi]_\lambda-[\Phi]_{\lambda\varphi}`,
:math:`\omega^*=[\omega]_\lambda-[\omega]_{\lambda\varphi}`,
:math:`\Phi'=\Phi-[\Phi]_\lambda`, :math:`u'=u-[u]_\lambda` and
:math:`\Delta_{EW}X=X|_{\lambda_2}-X|_{\lambda_1}`.

The east and west faces follow the structure of Brennan and Vincent (1980,
pp. 964-965), who write that wall as the flux :math:`u\Phi-u'\Phi'` evaluated
from west to east; expanding the east-minus-west difference gives
:math:`[u]_\lambda\Delta_{EW}\Phi'` together with a second contribution
carrying :math:`\Delta_{EW}u'`.  That second contribution is taken here with
the area departure :math:`\Phi^*` rather than the zonal mean
:math:`[\Phi]_\lambda`: a limited domain carries a net mass flux through its
walls, so with the full geopotential the term would change if an arbitrary
constant were added to :math:`\Phi`, whereas :math:`\Phi^*` leaves it
unchanged.  Both contributions vanish on a periodic (global) domain, where the
two walls coincide.

Michaelides (1987, p. 25) writes this wall as
:math:`c_1\iint([v]_\lambda\Phi^*)|_{\lambda_1}^{\lambda_2}d\varphi\,dp/g`.
Both factors are zonal means and therefore independent of longitude, which
makes that expression identically zero.

The north/south and top/bottom faces follow Muench (1965) and
Michaelides (1987), using the area departures :math:`\Phi^*` and
:math:`\omega^*`.  Brennan and Vincent (1980) write those two faces with the
full geopotential instead; in a limited domain that form additionally carries
the domain-mean overturning conversion proportional to
:math:`\int\overline{\omega}\,\alpha\,dp/g`, which is an interior
conversion rather than a boundary flux and is diagnosed separately as
:math:`C_{\mathrm{overturning}}`.
