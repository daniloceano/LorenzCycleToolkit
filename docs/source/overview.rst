Overview
========

What is the Lorenz Energy Cycle?
---------------------------------

The atmospheric circulation is driven by energy transformations between different forms:

- **Available Potential Energy (APE)**: Energy stored in atmospheric temperature gradients and density differences
- **Kinetic Energy (KE)**: Energy of atmospheric motion (wind)

Both APE and KE can be further divided into:

- **Zonal components** (Az, Kz): Represent large-scale, east-west averaged circulation patterns
- **Eddy components** (Ae, Ke): Represent deviations from the zonal mean, such as storms, cyclones, and weather disturbances

The Lorenz Energy Cycle quantifies how energy flows between these four reservoirs, helping us understand:

- How storms and cyclones intensify and decay
- Energy sources and sinks in weather systems  
- The role of different processes in atmospheric dynamics
- Differences between model simulations and observations

Mathematical Framework
----------------------

The Lorenz Energy Cycle (LEC), introduced by Edward Lorenz in 1955, is an analytical framework used to estimate atmospheric energy. It categorizes energy into zonal and eddy components of Kinetic Energy (Kz and Ke, respectively) and Available Potential Energy (Az and Ae, respectively). The LEC also quantifies conversions between these forms (Ca, Ce, Cz, and Ck), along with generation and dissipation terms (Gz, Ge, Dz, and De). Originally developed for global energetics, the framework has been adapted for regional studies, incorporating calculations for energy transport across boundaries (BAz, BAe, BKz, BKe).

The LEC budget is described by the following equations:

.. math::
   :nowrap:

   \begin{align*}
   \frac{\partial A_Z}{\partial t} &= -C_Z - C_A + G_Z + B A_Z \\
   \frac{\partial K_Z}{\partial t} &= -C_Z + C_K - D_Z + B K_Z + B \Phi_Z \\
   \frac{\partial A_E}{\partial t} &= C_A - C_E + G_E + B A_E \\
   \frac{\partial K_E}{\partial t} &= C_E - C_K - D_E + B K_E + B \Phi_E
   \end{align*}

Due to the difficulty in measuring friction terms for dissipation, both dissipation and generation are often computed as residuals from the budget equations:

.. math::
   :nowrap:

   \begin{align*}
   RG_Z &= G_Z + \varepsilon_{AZ} \\
   RG_E &= G_E + \varepsilon_{AE} \\
   RK_Z &= B \Phi_Z - D_Z + \varepsilon_{KZ} \\
   RK_E &= B \Phi_E - D_E + \varepsilon_{KE}
   \end{align*}

Where ε represents numerical errors. The complete cycle, assuming all terms are positive, is depicted below:

.. image:: https://github.com/daniloceano/lorenz-cycle/assets/56005607/d59eeb31-5cef-46ac-a841-1ba4170fafbd
   :width: 350

Getting Started
---------------

To begin using the LorenzCycleToolkit:

1. Follow the :doc:`installation` instructions to set up the toolkit
2. Configure your :doc:`configuration` files (namelist, box_limits, or track files)
3. See :doc:`usage` for command-line options and running the toolkit
4. Check :doc:`examples_and_tutorials` for practical step-by-step examples
5. Review :doc:`results` to understand the output files

For detailed mathematical formulations, see the :doc:`math` section.
