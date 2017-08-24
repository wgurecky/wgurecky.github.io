DGFE Spatial Discretization of the 1D Transport Equation
========================================================

.. author:: default
.. categories:: none
.. tags:: none
.. comments::



Background
===========

This work is an extension of the
methods introduced in the computational methods in radiation transport
graduate course offered at the University of Texas at Austin.

The theory
and consequences of applying the discrete ordinates method
and the multigroup approximation to transport equation
is left to the bulk of the course material. The text:
Computational Methods of Neutron Transport by E.E. Lewis (E. 2015)
also covers these basics in detail.

DGFE Theory
-----------

We begin with the within group :math:`g`, within
angle :math:`n` 1D transport equation. In the following section we do not write the
group subscript and angle superscript in order to reduce clutter.

.. math::
   \mu \frac{\partial}{\partial x} \psi(x) + \Sigma_t \psi(x) = S(x)
   :label: 1d_tr_dg

First, we multiply both sides of equation :eq:`1d_tr_dg` by a yet to be
defined test function, :math:`v(x)` resulting in:

.. math::
   \mu \frac{\partial}{\partial x} \psi(x) v(x) + \Sigma_t \psi(x)v(x) = S(x)v(x)
   :label: 1d_tr_dg2

Next the equation is integrated over the domain, :math:`x \in [0, L]`:

.. math::
   \int_0^L \mu \frac{\partial \psi(x)}{\partial x} v(x) dx + \int_0^L \Sigma_t \psi(x)v(x) dx =  \int_0^L S(x)v(x) dx
   :label: 1d_tr_weak

This is known as the weak form of equation :eq:`1d_tr_dg`.

Next, we integrate the first term in :eq:`1d_tr_weak` by parts giving:

.. math::
   \mu \psi(x)v(x)|_0^L- \mu \int_0^L  \psi(x) \frac{\partial v(x)}{\partial x} dx + \int_0^L \Sigma_t \psi(x)v(x) dx =  \int_0^L S(x)v(x) dx
   :label: 1d_tr_weak2

We are free to choose a functional form for :math:`v(x)`. In the
Galerkin approach, the test function is taken to be of the same
functional form as the solution approximation :math:`\psi(x)`. The
simple DGFE choice is to take :math:`v(x)` to be a linear combination of
ramp functions. Each ramp function :math:`h_{ei}(x)` is supported only
at one nodal location in the mesh and is defined to be zero at all other
nodes. The figure displays two interior neighboring ramp
functions which are each non-zero over element :math:`e_1`. The ramp
functions are defined to have unit height over their supporting node. An
element is defined to be the region between bounding nodes which are
points at which the solution is supported. For simplicity all elements
are assumed to have width :math:`\Delta x`.

.. figure:: images/single_ele.png
   :alt: Single finite element.
   :width: 5.00000cm

   Single finite element.

.. figure:: images/multi_ele.png
   :alt: Multiple DG finite elements.
   :width: 8.00000cm

   Multiple DG finite elements.

In practice the transport equation is integrated element-by-element,
rather than over the whole domain. The contribution of each element will
be summed together to recover the neutron balance over the whole domain.
For now, we consider an interior element, :math:`e_1` defined on the
sub-region: :math:`[a, b]`. The above figure shows the interior
element :math:`e_1` bounded by two other elements. Note that the
hypothetical DGFE numerical solution :math:`\psi` jumps in value at
element boundaries. As a consequence, at element boundaries the solution
is double valued at mesh edges. This is where the Discontinuous Galerkin
finite element scheme differs from the more commonly known Continuous
Galerkin (CG) FE spatial discretization method.

Now it is useful to formally define the ramp functions and their linear
combination. Over a single element, the solution :math:`\psi_e(x)` is
given by equation :eq:`sol_ele`.

.. math::
   \psi_e(x) = u_{eL}h_{e1}(x) + u_{eR}h_{e2}(x) = \sum_i u_{ei} h_{ei}(x),\ x\in[a,b]
   :label: sol_ele

Where :math:`i` is the edge index, in the one dimension case this
denotes either the left or right face. The ramp functions are given as:

.. math::
   h_{e1}(x) =
     \begin{cases}
                                      \frac{-1}{\Delta x}(x-a) + 1 & \text{, $x\in[a,b]$} \\
                                      0 & \text{, $otherwise$} 
     \end{cases}

and

.. math::
   h_{e2}(x) =
     \begin{cases}
                                      \frac{1}{\Delta x}(x-a) & \text{, $x\in[a,b]$} \\
                                      0 & \text{, $otherwise$} 
     \end{cases}

As previously stated, the Galerkin approach is to enforce :eq:`gal_asm`.

.. math::
   \psi_e(x) = v_e(x)
   :label: gal_asm

on each element. At first glance this appears this is an arbitrary
choice, and indeed, this assumption does not have to be made. One could
use different functional families for :math:`\psi` and :math:`v`,
however we will not investigate this option.

For this case where we have chosen simple ramp functions to represent
our 1D solution approximation, each element has two unknown scalar
values, :math:`\{u_{eL}, u_{eR}\}` that act to scale the ramp functions
over the element.

.. math::
   \mu \psi_e(x)v_e(x)|_a^b- \mu \int_a^b  \psi_e(x) \frac{\partial v_e(x)}{\partial x} dx + \int_a^b \Sigma_t \psi_e(x)v_e(x) dx =  \int_a^b S_e(x)v_e(x) dx
   :label: 1d_tr_weak_ele

Next we apply :eq:`sol_ele` and :eq:`gal_asm` to :eq:`1d_tr_weak_ele`. The
solution over the entire domain is the summation of the piecewise linear
solution approximation over all elements:

.. math:: \psi(x) = \sum_{e=0}^M \psi_e(x)

Where :math:`M` is the number of finite elements used.

The integral terms in equation :eq:`1d_tr_weak_ele` can be expanded to
explicitly show their dependence on the scaling factors. The second term
in :eq:`1d_tr_weak_ele` integrates to :eq:`w_e`.

.. math::
   W_e = -\int_a^b \mu \psi_e \frac{\partial v_e}{\partial x} dx = \frac{-\mu}{2}(u_{eR}^2 - u_{eL}^2) = 
   \frac{-\mu}{2} \mathbf u_e 
   \begin{bmatrix}
       -1      & 1 \\
       -1       & 1 
   \end{bmatrix}
   \mathbf u_e^T
   :label: w_e

With :math:`\mathbf u_e = [u_{eL}, u_{eR}]`. Note that this produces an
asymmetric element matrix. As a consequence, it is required that the
order of the nodes from left to right is preserved.

The third term in :eq:`1d_tr_weak_ele` integrates to :eq:`m_e`.

.. math::
   M_e = \int_a^b \Sigma_t \psi_e(x)v_e(x) dx =
   \frac{\Sigma_t \Delta x}{3} (u_L^2 + u_L u_R + u_R^2) = 
   \frac{\Sigma_t \Delta x}{3} \mathbf u_e 
   \begin{bmatrix}
       1      & 1/2 \\
       1/2      & 1 
   \end{bmatrix}
   \mathbf u_e^T
   :label: m_e

The RHS of equation :eq:`1d_tr_weak_ele` integrates to :eq:`s_e`.

.. math::
   RHS_e = \int_a^b S_e(x)v_e(x) dx =
   \frac{S_e \Delta x}{2} (u_L + u_R) = 
   \frac{S_e \Delta x}{2}
   \begin{bmatrix}
       1     \\
       1 
   \end{bmatrix}
   \mathbf u_e^T
   :label: s_e

Where we take the value :math:`S_e` to be the value of :math:`S_e(x)` at
the element mid-point. This is valid provided that :math:`S_e(x)` is a
linear function since this is equal to the average value of
:math:`S_e(x)` over the element.

Finally, we must deal with the boundary term which arose from
integrating the first term of equation :eq:`1d_tr_weak_ele` by parts.
This term is the only term which will contain information from
neighboring elements in its definition. This is why it is said that the
DGFE technique is “compact”. Let the outward normal at a given element
boundary to be denoted by :math:`\mathbf n`. The left side outward
normal for element :math:`e_1` is depicted in the figure:

.. figure:: images/bound_norm.png
   :alt: Outward normal on left face of element :math:`e_1`.
   :width: 8.00000cm

Outward normal on left face of element :math:`e_1`. As drawn,
:math:`\psi_{1L}^{\uparrow}=u_{e_0, 2}` and
:math:`\psi_{1L}^{\downarrow}=u_{e_1, 1}` in the figure.

It is useful to define a jump and average condition on an element
boundary. The average condition at the junction between two elements is
given by :eq:`avg`.

.. math:: 
    \{\{u\}\}_p = \frac{1}{2} (\lim_{x \to p^+} \psi(x) + \lim_{x \to p^-} \psi(x))
    :label: avg

Where the subscript :math:`p` denotes evaluation at a boundary. Since
:math:`\psi(x)|_p` and therefore :math:`u_p` is double valued at the
element boundaries; the limit approaching from the left is not equal to
the limit approaching from the right.

And the jump is provided by equation :eq:`jmp`.

.. math::
    [[u]]_p =  (\lim_{x \to p^+} \psi(x) - \lim_{x \to p^-} \psi(x))
    :label: jmp

Now it is useful define the “upwind” flux. According to :eq:`upwind`, the
sign of the dot product between the current neutron flow direction,
:math:`\mu` and the boundary normal vector :math:`\mathbf n_{e,p}` can
be used at each edge to determine the upwind flux value.

.. math::
   \psi^{\uparrow} = 
     \begin{cases}
         \psi_k|_p & \text{if $\mu \cdot \mathbf n_e|_p \leq 0$} \\
         \psi_e|_p & \text{if $\mu \cdot \mathbf n_e|_p > 0$} 
     \end{cases}
   :label: upwind

Where :math:`k` represents the neighboring element and :math:`e` is the
current element.

It is unclear what value to choose for the flux at the element
boundaries. This is required to evaluate
:math:`\mu \psi_e(x) v_e(x)|_a^b` in :eq:`1d_tr_weak_ele`.
The numerical flux is introduced
:math:`\mu \cdot \mathbf n \hat{F}` to resolve this issue. The boundary
term becomes :eq:`dg_fe_bound`.

.. math::
   \mu \psi_e(x) v_e(x)|_a^b = \mu \cdot \mathbf n \hat{F}  v_e(x)|_a^b 
   :label: dg_fe_bound

Upwind Formulation
~~~~~~~~~~~~~~~~~~

In this case, when evaluating :math:`\mu \psi_e(x) v(x)|_a^b`,
:math:`\psi(x)` always takes the upwind value at element boundaries.

.. math:: \mu \cdot \mathbf n \hat{F}  = \mu \cdot \mathbf n \psi^{\uparrow}

Equation :eq:`dg_fe_bound` can now be evaluated. If
:math:`\mu \cdot \mathbf n > 0`:

.. math::
   B_{ep_1} = \mu \cdot \hat{\mathbf n} \hat{F}  v(x)|_{p_1} = 
   \mu \cdot \mathbf n (u_e^2)|_p = 
   (\mu \cdot \mathbf n) \mathbf u_p 
   \begin{bmatrix}
       1      & 0 \\
       0      & 0 
   \end{bmatrix}
   \mathbf u_p^T

Where :math:`\mathbf u_p = [u_e, u_k]|_p`. Again, :math:`u_k|_p` is the
value of :math:`\psi` as the boundary from the neighboring element side
(i.e :math:`u_k=\lim_{x \to p^k}\psi(x)`) and likewise for the current
element side: :math:`u_e|_p=\lim_{x \to p^e}\psi(x)`.

If :math:`\mu \cdot \mathbf n \leq 0`:

.. math::
   B_{ep_2} = \mu \cdot \hat{\mathbf n} \hat{F}  v(x)|_{p_2} = 
   \mu \cdot \mathbf n (u_e \cdot u_k)|_p = 
   (\mu \cdot \mathbf n) \mathbf u_p 
   \begin{bmatrix}
       0      & 0 \\
       1      & 0 
   \end{bmatrix}
   \mathbf u_p^T

And the sum over both edges is given by :eq:`ele_sum`.

.. math::
   B_{ep_1} + B_{ep_2} = \mu \cdot \hat{\mathbf n} \hat{F}  v(x)|_a^b
   :label: ele_sum

Average Flux Formulation
~~~~~~~~~~~~~~~~~~~~~~~~

Alternatively, instead of simply taking the upwind flux value at each
element boundary, one may choose to use the average flux,
:math:`\{\{u\}\}_p` at each boundary. This results in the following:

If :math:`\mu \cdot \mathbf n \leq 0`:

.. math::
   B_{ep_1} = \mu \cdot \hat{\mathbf n} \hat{F}  v(x)|_{p_1} = 
   \mu \cdot \mathbf n \frac{u_e}{2} (u_e + u_k)|_p = 
   (\mu \cdot \mathbf n) \mathbf u_p 
   \begin{bmatrix}
       1/2     & 0 \\
       1/2     & 0 
   \end{bmatrix}
   \mathbf u_p^T

and If :math:`\mu \cdot \mathbf n > 0`:

.. math::
   B_{ep_2} = \mu \cdot \hat{\mathbf n} \hat{F}  v(x)|_{p_2} = 
   \mu \cdot \mathbf n  \frac{u_e}{2} (u_e + u_k)|_p = 
   (\mu \cdot \mathbf n) \mathbf u_p 
   \begin{bmatrix}
       1/2     & 0 \\
       1/2     & 0 
   \end{bmatrix}
   \mathbf u_p^T


System Matrix Construction and Boundary Conditions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For each element in the mesh we can write the neutron balance as:

.. math:: B_{ep_1} + B_{ep_2} + W_e + M_e = S_e

Multiplying both sides by :math:`[\mathbf u \mathbf I]^{-1}` we obtain:

.. math:: [b_{ep_1} + b_{ep_2} + w_e + m_e] \mathbf u^T = s_e

Where :math:`w_e \mathbf u^T= [\mathbf u \mathbf I]^{-1}W_e`,
:math:`s_e=[\mathbf u \mathbf I]^{-1} S_e`, and
:math:`m_e \mathbf u^T=[\mathbf u \mathbf I ]^{-1}M_e` and :math:`\mathbf I` is the
identity matrix. Collapsing further:

.. math:: [A_e] \mathbf u^T = s_e

The goal is to find the combination of the scaling factors,
:math:`\mathbf u=\{u_0, u_1, ...\}`, over all elements that best
satisfies the overall weak form of the neutron balance equation
:eq:`1d_tr_weak`. One can think of the finite element method in an
optimization context.  For more information see notes on the weighted residual method.

To assemble the global system matrix :math:`\mathbf A`, the individual element
matrices are “stamped” into :math:`\mathbf A`. Since each node in the
mesh is assigned a *unique ID* the elements of :math:`A_e` can be copied
into the global matrix :math:`\mathbf A`.

After :math:`\mathbf A` is constructed, the discretized, non-multiplying
transport equation can be written as:

.. math:: \mathbf A \mathbf u^T = \mathbf s

:math:`\mathbf A` is a sparse, non symmetric matrix. This linear system
of equations can be solved by GMRES or other iterative techniques.

Up until this point we have disregarded the application of boundary
conditions since we focused on the interior elements. For the first
order transport equation, all boundary conditions (vacuum, reflective,
white) can be described as either fixed or free. A Fixed boundary
condition specifies the value of :math:`\psi` at the boundary. This
arises in the vacuum case where inward facing ordinate fluxes are set
equal to zero at the boundary. This also arises in the reflective and
white cases where the banked outward fluxes from the previous scattering
source iteration are assigned as fixed boundary values for the inward
facing ordinate fluxes. A free boundary arises in cases where the flux
is allowed to escape from the domain. To implement a fixed boundary
condition, the row in the global system matrix, :math:`\mathbf A`, corresponding to the
boundary node is set equal to zero at all elements except at the
diagonal where the diagonal entry is set equal to 1. On the right hand
side the specified value for the flux at that node is set. Free boundary
conditions require no action.

Results
=======

In all cases the group structure boundaries of:

.. math:: [1.{E^-3} ,1.{E^-2}, 1.{E^-1}, 1.E0, 1.E1, 1.E2, 1.E3, 1.E4, 1.E5, 1.E6, 1.E7](eV)

were used to generate a 10-group cross section library. The infinite
dilution multigroup cross sections were generated with NJOY for this
work (al. 2017). For plotting, the group scalar fluxes are recovered
from the angle-dependant flux by the quadrature rule:

.. math:: \phi_g = \frac{1}{2}\sum_{n=1}^N w_n \psi_g^n(x)

For all presented results, :math:`S_8` Gauss-Legendre quadrature was
used for the angular flux decomposition by the discrete ordinates
method. Accordingly, the scattering cross section was approximated with
the first :math:`8` Legendre moments (thus retaining the first 8 terms
in the Legendre expansion of the scattering kernel). For consistency,
all cases were executed with 160 scattering source iterations to
converge the angle and energy neutron distribution.

The first result shown in figure below demonstrates the discontinuous
nature of the solution approximation. Neutrons are introduced on the left face traveling
to the right with an initial energy of :math:`1e7eV` and with a
source flux of :math:`1.27E6 [n/cm^2s]` into a 50cm thick graphite
block. The upwind formulation was used for the numerical flux. The
graphite was pure :math:`^{12}C` with a density of 2.23\ :math:`[g/cc]`.
The case was executed using a relatively coarse 20 element spatial mesh
for visual clarity of the discontinuities.

.. figure:: results/scflux_graphite_beam_1.png
   :alt: Group scalar fluxes for a high energy beam incident on a graphite block. The :math:`y` axis units are in :math:`n/cm^2s`. Upwind numerical flux used.
   :width: 12.00000cm

The same problem was re-run this time with the average numerical flux
formulation. This resulted in the following figure.

.. figure:: results/scflux_graphite_beam_2.png
   :alt: Group scalar fluxes for a high energy beam incident on a graphite block. The :math:`y` axis units are in :math:`n/cm^2s`. Average numerical flux used at element boundaries.
   :width: 12.00000cm

Interestingly, for this problem it appears the upwind strategy provides
a more accurate result. Qualitatively, the expected far-field
exponential decay of the highest energy group flux is more accurately
captured by the upwind flux formulation.

In the next case, a thin (0.5:math:`[mm]`) sheet of highly absorptive
pure :math:`^{10}B` with a density of 5\ :math:`[g/cc]` was inserted
into the graphite block at 15\ :math:`[cm]`. Shown in figures below,
this effectively eliminated the majority of the thermal neutron current
passing through the region resulting in a sharp dip in thermal flux near
the sheet, followed by thermal neutron recovery further away since there
are still neutrons down scattering into the lower energy groups over the
whole domain. Expectedly, the boron had little influence on the higher
energy groups. The first case was executed with 20 elements followed by
a fine mesh run with 60 elements.

.. figure:: results/scflux_graphite_beam_3.png
   :alt: Coarse mesh solution.
   :width: 12.00000cm

.. figure:: results/scflux_graphite_beam_4.png
   :alt: Fine mesh solution.
   :width: 12.00000cm

Conclusion
==========

The DGFE method was introduced and implemented in 1D. The current
implementation could serve as a starting point to more detailed
investigations.

It was shown that DGFE allows for a flexible definition of the numerical
flux and that this choice has a significant impact on the resulting
numerical approximation.

Improving the order of accuracy of the finite element discretization is
a potential avenue for future work. This would involve increasing the
polynomial order of the ramp basis functions over each element from 1 to
2.

Others have shown that the DGFE method “locks” in the optically thick
diffusion limit, meaning, the flux is artificially depressed in regions
that are highly opaque and highly diffusive to neutrons. For most
practical problems this is not a concern, however, it could be
interesting to investigate the work performed by J. Guermond et. al
(2014) (Guermond, Kanschat, and Ragusa 2014) on this subject. Guermond
et. al. present a method to adaptively choose between the unwinding and
averaging formulation in each element independently based on the local
scattering cross section and cell width. This has been shown to
effectively eliminate this issue with the DGFE method without
significant additional computational overhead.

The code is available online at https://github.com/wgurecky/spyTran.

References
===========

.. raw:: html

   <div id="refs" class="references">

.. raw:: html

   <div id="ref-mac17">

al., R. Macfarlane et. 2017. “The Njoy Nuclear Data Processing System.”
*Los Alamos National Laboratory (LANL)* LA-UR-17-20093.

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Lewis">

E., Lewis. 2015. *Numerical Methods for Radiation Transport*. CRC Press.

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Guermond2014">

Guermond, Jean-Luc, Guido Kanschat, and Jean C. Ragusa. 2014.
“Discontinuous Galerkin for the Radiative Transport Equation.” In
*Recent Developments in Discontinuous Galerkin Finite Element Methods
for Partial Differential Equations: 2012 John H Barrett Memorial
Lectures*, edited by Xiaobing Feng, Ohannes Karakashian, and Yulong
Xing, 181–93. Cham: Springer International Publishing.
doi:\ `10.1007/978-3-319-01818-8\_7 <https://doi.org/10.1007/978-3-319-01818-8_7>`__.

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-lesaint">

Lesaint, P., and P. Raviart. 1974. “On a Finite Element Method for
Solving the Neutron Transport Equation.” *Mathmatical Aspects of Finite
Elements in Partial Differential Equations* 33.

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-reed">

Reed, W., and T. Hill. 1973. “Triangular Mesh Methods for the Neutron
Transport Equation.” *Los Alamos National Lab* LA-UR-73-479.

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-riviere">

Riviere, B. 2008. “Discontinuous Galerkin Method for Solving Elliptic
and Parabolic Equations: Theory and Implementation.” *SIAM Frontiers in
Applied Mathematics*.

.. raw:: html

   </div>

.. raw:: html

   </div>
