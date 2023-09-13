##################
Governing equation
##################

I consider a domain which is filled by an incompressible liquid, which is governed by the following equations

.. math::

   \pder{u_i}{x_i}
   =
   0,

.. math::

   \rho^f
   \left(
      \pder{u_i}{t}
      +
      u_j \pder{u_i}{x_j}
   \right)
   =
   \pder{\sigma_{ij}}{x_j}
   +
   \rho^f
   g_i
   +
   \rho^f
   a_i^{IBM}.

Here :math:`a_i^{IBM}` is a vector field, whose concrete formula is not necessary here.
The only one important thing is that its role is to **constraint the behaviour of the surrounding fluid as if there were an object inside the domain**.
This is to be determined later, such that a desired boundary condition is fulfilled on the surface of the object.

On the other hand, I consider another domain, in which an object sits and surrounded by liquid.
The governing equation describing the behaviour of this one is the Newton's law:

.. math::

   \rho^p V^p \tder{U_i^p}{t}
   =
   \sum_{\forall} \left( \text{force} \right)_i.

I consider two contributions to the right-hand-side forcing term: the hydrodynamic force acting on the particle:

.. math::

   \int_{\partial V^p} \sigma_{ij} n_j dS,

and the gravitational force:

.. math::

   \int_{\partial V^p} \rho^p g_i dV.

Now I consider to couple the fluid behaviour and the particulate motion.
To do so, I integrate the momentum balance of the fluid inside the same volume containing this object:

.. math::

   \int_{V^p}
   \left\{
      \rho^f
      \left(
         \pder{u_i}{t}
         +
         u_j \pder{u_i}{x_j}
      \right)
      -
      \pder{\sigma_{ij}}{x_j}
      -
      \rho^f
      g_i
      -
      \rho^f
      a_i^{IBM}
   \right\}
   dV
   =
   0,

giving

.. math::

   \int_{\partial V^p} \sigma_{ij} n_j dS
   =
   \int_{V^p}
   \left\{
      \rho^f
      \left(
         \pder{u_i}{t}
         +
         u_j \pder{u_i}{x_j}
      \right)
      -
      \rho^f
      g_i
      -
      \rho^f
      a_i^{IBM}
   \right\}
   dV,

where the Gauss theorem is used.
Thus I have

.. math::

   \rho^p V^p \tder{U_i^p}{t}
   =
   \int_{V^p}
   \left\{
      \rho^f
      \left(
         \pder{u_i}{t}
         +
         u_j \pder{u_i}{x_j}
      \right)
      -
      \rho^f
      g_i
      -
      \rho^f
      a_i^{IBM}
   \right\}
   dV
   +
   \int_{\partial V} \rho^p g_i dV^p,

which is simplified as

.. math::

   \rho^p V^p \tder{U_i^p}{t}
   =
   \int_{V^p}
   \left\{
      \rho^f
      \left(
         \pder{u_i}{t}
         +
         u_j \pder{u_i}{x_j}
      \right)
      -
      \rho^f
      a_i^{IBM}
   \right\}
   dV
   +
   \left(
      \rho^p
      -
      \rho^f
   \right)
   g_i V^p,

In summary, a set of the governing equations describing the dynamics of the whole system results in

.. math::

   \pder{u_i}{x_i}
   =
   0,

.. math::

   \pder{u_i}{t}
   +
   u_j \pder{u_i}{x_j}
   =
   \pder{\sigma_{ij}}{x_j}
   +
   g_i
   +
   a_i^{IBM},

and

.. math::

   \rho^p V^p \tder{U_i^p}{t}
   & =
   \int_{V^p}
   \left(
      \pder{u_i}{t}
      +
      u_j \pder{u_i}{x_j}
   \right)
   dV \\
   & -
   \int_{V^p}
   a_i^{IBM}
   dV \\
   & +
   \left(
      \rho^p
      -
      1
   \right)
   g_i V^p,

where I fix :math:`\rho^f \equiv 1` for notational simplicity, and as a result :math:`\rho^p` is now used to tell the density ratio.

