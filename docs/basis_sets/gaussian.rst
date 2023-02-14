Gaussian Orbitals
=================

A class to represent orbitals built off of linear combinations of Gaussians. These have the following form.

.. math::

   G_{l,m}(\mathbf{r}) = \sum_{k=1}^{N} \frac{\left(2\alpha_k\right)^{3/4}}{2\left(\pi^{1/4}\right)}\sqrt{\frac{\left(8\alpha_k\right)^l l! (2l + 1)}{(2l)!}} r^l Y_l^m(\theta, \phi) e^{-\alpha_k r^2}

.. cpp:namespace:: compchem

.. cpp:class:: GaussianOrbital : public BasisOrbital

    Represents a linear combination of Gaussian orbital functions.

    .. cpp:member:: std::vector<double> coefs

        Holds the coefficients for each term.

    .. cpp:member:: std::vector<double> alphas

        Holds the stretching factors for each term.

    .. cpp:member:: int l

        The angular momentum number for the orbital.

    .. cpp:member:: int ml

        The magnetic number for the orbital.

    .. cpp:member:: Polynomial<3> *harms

        The spherical harmonic for the orbital.

    .. cpp:function:: GaussianOrbital(int l, int ml, const std::vector<double> &coefs, const std::vector<double> &alphas)

        :param l: The angular momentum quantum number.
        :param ml: The magnetic quantum number.
        :param coefs: The coefficients for the linear combination.
        :param alphas: The stretch factor for each term.

        Constructs a new Gaussian orbital.

    .. cpp:function:: GaussianOrbital(const GaussianOrbital &copy)

        :param copy: The orbital to copy.

        Creates a copy of the given orbital.

    .. cpp:function:: virtual ~GaussianOrbital()

        The destructor.

    .. cpp:function:: const std::vector<double> &getcoefs() const

        :return: The coefficients for the linear combination.

    .. cpp:function:: const std::vector<double> &getalphas() const

        :return: The stretch factors for the terms.

    .. cpp:function:: double getcoef(int index) const

        :param index: The index of the term.
        :return: The coefficient for the term at the index.

    .. cpp:function:: double getalpha(int index) const

        :param index: The index of the term.
        :return: The stretch factor for the term.

    .. cpp:function:: int getnterms() const

        :return: The number of terms in the linear combination.

    .. cpp:function:: int getl() const

        :return: The angular momentum quantum number.

    .. cpp:function:: int getml() const

        :return: The magnetic quantum number.

    .. cpp:function:: const Polynomial<3> &getharms() const

        :return: The spherical harmonic part of the orbital.

    .. cpp:function:: double eval(double x, double y, double z) const override

        :param x, y, z: The coordinates to evaluate at.
        :return: The value of the orbital at the given point.

    .. cpp:function:: BasisOrbital *copy() const override

        :return: A copy of this orbital.
