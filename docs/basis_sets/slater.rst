Slater-type Orbitals
====================

Represents a Slater-type orbital. These are orbitals that look like this.

.. math::

   S_{n,l,m}(\mathbf{r}) = (2\zeta)^{n} \sqrt{\frac{2\zeta}{(2n)!}} r^{n-1} Y_l^m(\theta, \phi) e^{-\zeta r}

.. cpp:namespace:: compchem

.. cpp:class:: SlaterOrbital : public BasisOrbital

    .. cpp:member:: double Zeff

        The effective atomic number, usually found using Slater's rules.

    .. cpp:member:: int n

        The energy quantum number of the orbital.

    .. cpp:member:: int l

        The angular momentum quantum number.

    .. cpp:member:: int ml

        The magnetic quantum number.

    .. cpp:member:: Polynomial<3> *harms

        The spherical harmonic part of the orbital.

    .. cpp:function:: SlaterOrbital(double zeff, int n, int l, int ml)

        :param zeff: The effective atomic number of the orbital.
        :param n: The orbital energy quantum number.
        :param l: The angular momentum quantum number
        :param ml: The magnetic quantum number.

    .. cpp:function:: SlaterOrbital(const SlaterOrbital &copy)

        :param copy: The orbital to copy.

        Creates a copy of a Slater-type orbital.

    .. cpp:function:: virtual ~SlaterOrbital()

        The destructor.

    .. cpp:function:: double eval(double x, double y, doubel z) const override

        :param x, y, z: The position to evaluate the orbital at.
        :return: The value of the orbital at the given position.

    .. cpp:function:: double getZeff() const

        :return: The effective atomic number of the orbital.

    .. cpp:function:: int getn() const

        :return: The orbital energy quantum number of the orbital.

    .. cpp:function:: int getl() const

        :return: The angular momentum quantum number.

    .. cpp:function:: int getml() const

        :return: The magnetic quantum number.

    .. cpp:function:: const Polynomial<3> &getharms() const

        :return: The spherical harmonic part of the orbital.

    .. cpp:function:: BasisOrbital *copy() const override

        :return: A copy of this orbital.

Non Member Functions
--------------------

.. cpp:function:: double slater_rule(int n, int l, const GSConfig &conf)

    :param n: The orbital energy quantum number of the orbital.
    :param l: The angular momentum quantum number.
    :param conf: The ground-state configuration of the atom.
    :return: The effective atomic number for a shell, computed using Slater's rules.



    
