Rotor Calculation
=================

Computes the rotational constants and rotor type of a molecule.

.. cpp:namespace:: compchem

.. cpp:enum:: rotor_t

    Represents the kind of rotor.

    .. cpp:enumerator:: SPERICAL

        A spherical rotor.

    .. cpp:enumerator:: LINEAR

        A linear rotor.

    .. cpp:enumerator:: OBLATE

        An oblate rotor.

    .. cpp:enumerator:: PROLATE

        A prolate rotor.

    .. cpp:enumerator:: ASSYMETRIC

        An assymetric rotor.


.. cpp:function:: rotor_t comprotor(const Molecule &mol, std::array<double, 3> *out)

    :param mol: The molecule to compute on.
    :param out: After exiting, this parameter will contain the rotational constants in :math:`\textbf{J}`.
    :return: The type of rotor.
