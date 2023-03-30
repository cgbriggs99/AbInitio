Base Wavefunction
=================

This is the base class representation of a wavefunction.

.. cpp:namespace:: compchem

.. cpp:class:: Wavefunction

    Represents a wavefunction.

    .. cpp:member:: protected int electrons

        The number of electrons represented by the wavefunction.

    .. cpp:member:: protected int multiplicity

        The multiplicity of the system represented by the wavefunction.

    .. cpp:function:: Wavefunction(int electrons)

        Initializes the wavefunction with the number of electrons and a singlet multiplicity.

    .. cpp:function:: Wavefunction(int electrons, int multiplicity)

        Initializes the wavefunction with the number of electrons and a given multiplicity.

    .. cpp:function:: virtual ~Wavefunction() = default

        Destructor.

    .. cpp:function:: virtual const double *getoverlap(int *dim = nullptr) const = 0

        :param dim: The output for the number of rows in the matrix.
        :return: The overlap matrix for this wavefunction.

    .. cpp:function:: virtual const double *getkinetic(int *dim = nullptr) const = 0

        :param dim: The output for the number of rows in the matrix.
        :return: The kinetic energy matrix for this wavefunction.

    .. cpp:function:: virtual const double *getpotential(int *dim = nullptr) const = 0

        :param dim: The output for the number of rows in the matrix.
        :return: The potential energy matrix for this wavefunction.

    .. cpp:function:: virtual const double *getfock(int *dim = nullptr) const = 0

        :param dim: The output for the number of rows in the matrix.
        :return: The Fock matrix.

    .. cpp:function:: virtual int getnorbs() const = 0

        :return: The number of orbitals in this wavefunction. This is usually what is placed in the :cpp:expr:`dim` variables of the other functions.

    .. cpp:function:: virtual int getelectrons() const

        :return: The number of electrons in the orbital.

    .. cpp:function:: virtual int getmultiplicity() const

        :return: The multiplicity of the wavefunction.
