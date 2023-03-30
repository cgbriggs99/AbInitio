Basis Sets
==========

Contains a base class for basis sets.

.. cpp:namespace:: compchem

.. cpp:class:: BasisOrbital

    A class to represent orbitals.

    .. cpp:function:: virtual ~BasisOrbital() = default

        The default destructor.

    .. cpp:function:: virtual double eval(double x, double y, double z) const = 0

        :param x,y,z: The coordinates for the evaluation.
        :return: The value of the orbital at the given point.

    .. cpp:function:: virtual BasisOrbital *copy() const = 0

        :return: A copy of the orbital.

    .. cpp:function:: virtual double laplacian(double x, double y, double z) const = 0

        :param x,y,z: The coordinates for the evaluation.
        :return: The Laplacian of the orbital at the point.

        This method should return :math:`\nabla_{\mathbf{r}}^2\psi\left(\mathbf{r}\right)`, not the kinetic energy.
