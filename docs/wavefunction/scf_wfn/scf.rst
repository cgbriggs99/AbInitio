SCF Base Wavefunction
=====================

The base class for SCF wavefunctions.

.. cpp:namespace:: compchem

.. cpp:class:: SCFWfn : public compchem::Wavefunction

    Base class for SCF wavefunctions.

    .. cpp:member:: protected int dim

        The number of rows in the arrays.

    .. cpp:member:: protected const double *S

        A pointer to the overlap matrix.

    .. cpp:member:: protected const double *T

        A pointer to the kinetic energy matrix.

    .. cpp:member:: protected const double *V

        A pointer to the potential energy matrix.

    .. cpp:member:: protected const TEIArray *tei

        A pointer to the electron repulsion integrals.

    .. cpp:member:: protected double *Ca
    .. cpp:member:: protected double *Cb

        Contains the coeficient matrices for alpha and beta spins.


    .. cpp:member:: protected double *Da
    .. cpp:member:: protected double *Db

        Contains the alpha and beta spin density matrices.

    .. cpp:member:: protected double *Fa
    .. cpp:member:: protected double *Fb

        Contains the alpha and beta spin Fock matrices.

    .. cpp:member:: protected double energy

        Represents the total energy.

    .. cpp:member:: protected double *es

        Contains the energies of each orbital.

    .. cpp:member:: protected int orbs

        Contains the number of orbitals.

    
