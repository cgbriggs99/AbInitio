Integral Method Base
====================

This represents the base class for integral methods. Not all integral methods need to implement this class, but any integral method must have methods with names that match.

.. cpp:namespace:: compchem

.. cpp:class:: IntegralMethod

    Base class for integral computers.

    .. cpp:member:: protected compchem::OptionList &opts

        Options to pass to the method.

    .. cpp:function:: IntegralMethod(const OptionList &opts)

        Create an instance of a calculator with the given options.

    .. cpp:function:: IntegralMethod()

        Create an instance of a calculator with a copy of the :cpp:class:`compchem::GlobalOptions` as the options.

    .. cpp:function:: virtual ~IntegralMethod()

        Destructor.

    .. cpp:function:: virtual double overlap(const compchem::BasisOrbital *o1, const compchem::BasisOrbital *o2, std::array<double, 3> center1, std::array<double, 3> center2) const = 0

        :params o1, o2: The orbitals to integrate.
        :params center1, center2: The centers for each orbital.
        :return: The overlap between the two orbitals.

        Computes the overlap, or :math:`S`, integral between two orbitals.

    .. cpp:function:: virtual double kinetic(const compchem::BasisOrbital *o1, const compchem::BasisOrbital *o2, std::array<double, 3> center1, std::array<double, 3> center2) const = 0

        :params o1, o2: The orbitals to integrate.
        :params center1, center2: The centers for each orbital.
        :return: The kinetic energy between the two orbitals.

        Computes the kinetic energy, or :math:`T`, integral between two orbitals. The Laplacian is taken on the second orbital, though due to the symmetry of the Laplacian operator, it does not matter which is used.

    .. cpp:function:: virtual double attraction(const compchem::BasisOrbital *o1, const compchem::BasisOrbital *o2, std::array<double, 3> center1, std::array<double, 3> center2, const compchem::Atom &atom) const = 0

        :params o1, o2: The orbitals to integrate.
        :params center1, center2: The centers for each orbital.
        :param atom: The atom for the charge and center for the integral.
        :return: The attraction of an electron to the specified atom.

        Computes the Coulombic attraction between an electron and a nucleus.

    .. cpp:function:: virtual double repulsion(const compchem::BasisOrbital *o1, const compchem::BasisOrbital *o2, const compchem::BasisOrbital *o3, const compchem::BasisOrbital *o4, std::array<double, 3> center1, std::array<double, 3> center2, std::array<double, 3> center3, std::array<double, 3> center4) const = 0

        :params o1, o2, o3, o4: The orbitals to integrate.
        :params center1, center2, center3, center4: The centers for each orbital.
        :return: The repulsion integral between two electrons.

        Computes the two-electron repulsion integral, represented as :math:`(ab|cd)`.
