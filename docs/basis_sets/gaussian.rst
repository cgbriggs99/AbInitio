Gaussian Orbitals
=================

A class to represent orbitals built off of linear combinations of Gaussians. These have the following form.

.. math::

   G_{l,m}(\mathbf{r}) = \sum_{k=1}^{N} C_k \frac{\left(2\alpha_k\right)^{3/4}}{2\left(\pi^{1/4}\right)}\sqrt{\frac{\left(8\alpha_k\right)^l l! (2l + 1)}{(2l)!}} r^l Y_{l,m}(\theta, \phi) e^{-\alpha_k r^2}

In this definition, :math:`N` is the number of terms in the linear combination represented by :cpp:member:`compchem::GaussianOrbital::size`; :math:`C_k` is the coefficient for the :math:`k` th term, which is stored in :cpp:member:`compchem::GaussianOrbital::coefs`; :math:`\alpha_k` is the stretching factor for the :math:`k` th term, which is stored in :cpp:member:`compchem::GaussianOrbital::alphas`; :math:`\frac{\left(2\alpha_k\right)^{3/4}}{2\left(\pi^{1/4}\right)}\sqrt{\frac{\left(8\alpha_k\right)^l l! (2l + 1)}{(2l)!}}` is the normalization factor for the :math:`k` th term, which is stored in :cpp:member:`compchem::GaussianOrbital::norms`; :math:`r^l Y_{l,m}(\theta, \phi)` is the spherical harmonic for the orbital, and can be represented as a Cartesian polynomial, which is stored in :cpp:member:`compchem::GaussianOrbital::harms`; :math:`l` is the angular momentum quantum number which is stored as :cpp:member:`compchem::GaussianOrbital::l`; and :math:`m` is the magnetic quantum number, which is stored as :cpp:member:`compchem::GaussianOrbital::ml`.

.. cpp:namespace:: compchem

.. cpp:class:: GaussianOrbital : public BasisOrbital

    Represents a linear combination of Gaussian orbital functions.

    .. cpp:member:: double *coefs

        Holds the coefficients for each term.

    .. cpp:member:: double *alphas

        Holds the stretching factors for each term.

    .. cpp:member:: double *norms

        Holds the normalization factors for each term. This is calculated by the constructor, and is useful to have separate for some calculations.

    .. cpp:member:: int l

        The angular momentum quantum number for the orbital.

    .. cpp:member:: int ml

        The magnetic quantum number for the orbital.

    .. cpp:member:: Polynomial<3> *harms

        The spherical harmonic for the orbital.

    .. cpp:member:: int size;

        The number of terms in the linear combination.

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

        The destructor. It is virtual to allow for extensions.

    .. cpp:function:: const double *getcoefs() const

        :return: The coefficients for the linear combination.

    .. cpp:function:: const double *getalphas() const

        :return: The stretch factors for the terms.

    .. cpp:function:: double getcoef(int index) const

        :param index: The index of the term.
        :return: The coefficient for the term at the index.
        :raises out_of_range: Raises :cpp:class:`out_of_range` when passed an index less than zero or greater than or equal to the size of the linear combination.

    .. cpp:function:: double getalpha(int index) const

        :param index: The index of the term.
        :return: The stretch factor for the term.
        :raises out_of_range: Raises :cpp:class:`out_of_range` when passed an index less than zero or greater than or equal to the size of the linear combination.


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

    .. cpp:function:: double laplacian(double x, double y, double z) const override

        :param x, y, z: The coordinates to find the Laplacian at.
        :return: The value of the Laplacian of the orbital at the given point.

        Overrides :cpp:func:`compchem::BasisOrbital::laplacian`. Returns the Laplacian, and not the kinetic energy.

    .. cpp:function:: BasisOrbital *copy() const override

        :return: A copy of this orbital.

    .. cpp:function:: GaussianOrbital &operator=(const GaussianOrbital &other)

        :param other: The right-hand side of the assignment.
        :return: A reference to :cpp:expr:`this`.

    .. cpp:function:: private void sort()

        Sorts the terms of the linear combination so that the terms with the largest :math:`\alpha` are first.


Non-member Functions
====================

.. cpp:function:: static void sort_internal(double *alphas, double *coefs, int size)

    :param alphas: The array of alphas to sort by.
    :param coefs: The array of coefs, which get swapped around to keep with their alphas.
    :param size: The size of the arrays.

    Does the actual sorting of the arrays for :cpp:func:`compchem::GaussianOrbital::sort`.

