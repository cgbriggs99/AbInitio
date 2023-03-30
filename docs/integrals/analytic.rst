Analytic Integration
====================

Methods to compute analytic integrals of Gaussian orbitals.

.. cpp:namespace:: compchem

.. cpp:class:: AnalyticIntegral

    Computes integrals of Gaussian orbitals analytically. Does not extend :cpp:class:`compchem::IntegralMethod`, since this class needs to ensure the inputs are Gaussian orbitals.

    .. cpp:member:: protected compchem::OptionList &opts

        Options to pass to the integration routines.

    .. cpp:function:: AnalyticIntegral(const OptionList &opts)

        Initializes the options with a copy of the passed list.

    .. cpp:function:: AnalyticIntegral()

        Initializes the options with a copy of the global options.

    .. cpp:function:: virtual double overlap(const compchem::GaussianOrbital *o1, const compchem::GaussianOrbital *o2, std::array<double, 3> center1, std::array<double, 3> center2) const

        :params o1, o2: The orbitals to integrate.
        :params center1, center2: The centers for each orbital.
        :return: The overlap between the two orbitals.

        Computes the overlap, or :math:`S`, integral between two orbitals.

    .. cpp:function:: virtual double kinetic(const compchem::GaussianOrbital *o1, const compchem::GaussianOrbital *o2, std::array<double, 3> center1, std::array<double, 3> center2) const

        :params o1, o2: The orbitals to integrate.
        :params center1, center2: The centers for each orbital.
        :return: The kinetic energy between the two orbitals.

        Computes the kinetic energy, or :math:`T`, integral between two orbitals. The Laplacian is taken on the second orbital, though due to the symmetry of the Laplacian operator, it does not matter which is used.

    .. cpp:function:: virtual double attraction(const compchem::GaussianOrbital *o1, const compchem::GaussianOrbital *o2, std::array<double, 3> center1, std::array<double, 3> center2, const compchem::Atom &atom) const

        :params o1, o2: The orbitals to integrate.
        :params center1, center2: The centers for each orbital.
        :param atom: The atom for the charge and center for the integral.
        :return: The attraction of an electron to the specified atom.

        Computes the Coulombic attraction between an electron and a nucleus.

    .. cpp:function:: virtual double repulsion(const compchem::Gaussianorbital *o1, const compchem::GaussianOrbital *o2, const compchem::GaussianOrbital *o3, const compchem::GaussianOrbital *o4, std::array<double, 3> center1, std::array<double, 3> center2, std::array<double, 3> center3, std::array<double, 3> center4) const

        :params o1, o2, o3, o4: The orbitals to integrate.
        :params center1, center2, center3, center4: The centers for each orbital.
        :return: The repulsion integral between two electrons.

        Computes the two-electron repulsion integral, represented as :math:`(ab|cd)`.

    .. cpp:function:: double boys_square(int j, double T) const

        :param j: The order of the Boys integral.
        :param T: The scaling of the Boys integral.
        :return: The Boys integral.

        Computes the Boys integral.

        .. math::

            F_j(T) = \int_0^1 u^{2j} e^{-Tu^2} du = \frac{1}{2T^{j + \frac{1}{2}}}\gamma\left(j + \frac{1}{2}, T\right)

        When the ``ANALYTIC BOYS`` option is cleared, the integral will be computed numerically, rather than with the analytic solution. Otherwise, if :math:`T` is large, then the solution can be represented as

        .. math::

            F_j(T) = \frac{1}{2T^{j + \frac{1}{2}}} \left(\Gamma\left(j + \frac{1}{2}\right)\mathrm{erf}\left(\sqrt{T}\right) - (-1)^{j - 1}e^{-T}\sqrt{T}\sum_{k = 0}^{j - 1} \left(\frac{1}{2} - j\right)_{j - k - 1}(-T)^k\right)

        This representation has stability problems for small :math:`T`, so a secondary sum is used if :math:`T < 0.001`, which is

        .. math::

            F_j(T) = \sum_{k = 0}^{\infty} \frac{\left(-T\right)^k}{\left(2k + 2j + 1\right)k!}

        **Options**

        ``ANALYTIC BOYS``: *default:* :cpp:expr:`true`.

        If true, compute the Boys integral using a series. If false, use numeric integration.

        ``BOYS POINTS``: *default:* :cpp:expr:`32`

        The number of points to use when calculating the numeric integral for the Boys function. Only referenced if ``ANALYTIC BOYS`` is false.


Protected Methods
^^^^^^^^^^^^^^^^^

    .. cpp:function:: protected double os_attr(const std::array<int, 7> &index, std::map<std::array<int, 7>, double> &ints, const std::array<double, 3> &c1, const std::array<double, 3> &c2, const std::array<double, 3> &ca, double Rx, double Ry, double Rz, double Px, double Py, double Pz, double zeta) const

        :param index: The index of the integral.
        :param ints: A collection of integrals already computed.
        :params c1, c2: The centers for the first and second orbitals.
        :param ca: The position of the atom.
        :params Rx, Ry, Rz: Components of the one-center offset of the integral between the orbitals and the atom.
        :params Px, Py, Pz: Components of the offset of the integral between only the orbitals.
        :param zeta: The sum of the stretching constants for the two Gaussian functions.
        :return: The value of the Obara-Saika intermediate with the specified index.

        The Obara-Saika intermediates are calculated using the recursion formula

        .. math::

            \left[\mathbf{a}\middle|\mathbf{b}\right]^N = \left(P_i - A_i\right) \left[\mathbf{a} - \mathbf{1}_i\middle|\mathbf{b}\right]^N + \frac{1}{2\zeta}\left(\left(a_i - 1\right)\left[\mathbf{a} - \mathbf{2}_i\middle|\mathbf{b}\right]^N + b_i\left[\mathbf{a} - \mathbf{1}_i\middle|\mathbf{b} - \mathbf{1}_i\right]^N\right) \\
            - \left(P_i - C_i\right)\left[\mathbf{a} - \mathbf{1}_i\middle|\mathbf{b}\right]^{N + 1} - \frac{1}{2\zeta}\left(\left(a_i - 1\right)\left[\mathbf{a} - \mathbf{2}_i\middle|\mathbf{b}\right]^{N + 1} + b_i\left[\mathbf{a} - \mathbf{1}_i\middle|\mathbf{b} - \mathbf{1}_i\right]^{N + 1}\right)

        .. math::

            \left[\mathbf{0}\middle|\mathbf{b}\right]^N = \left(P_i - B_i\right)\left[\mathbf{0}\middle|\mathbf{b} - \mathbf{1}_i\right]^N + \frac{1}{2\zeta}\left(b_i - 1\right) \left[\mathbf{0}\middle|\mathbf{b} - \mathbf{2}_i\right]^N \\
            - \left(P_i - C_i\right)\left[\mathbf{0}\middle|\mathbf{b} - \mathbf{1}_i\right]^{N+1} - \frac{1}{2\zeta} \left(b_i - 1\right) \left[\mathbf{0}\middle|\mathbf{b} - \mathbf{2}_i\right]^{N+1}

        .. math::

            \left[\mathbf{0}\middle|\mathbf{0}\right]^N = F_N\left(\zeta R^2\right)

        Where :math:`i = x, y, z`; :math:`\mathbf{1}_i = \left(\delta_{ix}, \delta_{iy}, \delta_{iz}\right)` and :math:`\mathbf{2}_i = 2\mathbf{1}_i`; and where :math:`A_x, A_y, A_z` are stored in :cpp:expr:`c1`; :math:`B_x, B_y, B_z` are stored in :cpp:expr:`c2`; :math:`C_x, C_y, C_z` are stored in :cpp:expr:`ca`; :math:`R_x, R_y, R_z` are stored in :cpp:expr:`Rx, Ry, Rz`; :math:`R^2 = R_x^2 + R_y^2 + R_z^2`; :math:`P_x, P_y, P_z` are stored in :cpp:expr:`Px, Py, Pz`; :math:`\zeta` is stored in :cpp:expr:`zeta`; and the layout for :cpp:expr:`index` is :math:`\left\{a_x, a_y, a_z, b_x, b_y, b_z, N\right\}`. The Boys function, :math:`F_j(T)`, is computed using :cpp:func:`compchem::AnalyticIntegral::boys_square`.

    .. cpp:function:: double attr_integral(const int *pows1, const int *pows2, const std::array<double, 3> &c1, const std::array<double, 3> &c2, const compchem::GaussianOrbital *o1, const compchem::GaussianOrbital *o2, const compchem::Atom &atom) const

        :params pows1, pows2: The exponents on the given Cartesian orbital's radial part.
        :params c1, c2: The centers of the orbitals.
        :params o1, o2: The orbitals for the integrals.
        :param atom: The atom for the attraction integral.
        :return: The attraction integral for the given Cartesian Gaussian orbital.

        Worker function for the attraction integral.

