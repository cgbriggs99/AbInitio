molecule.hpp
============

Contains a class to represent a molecule.

.. cpp:namespace:: compchem

.. cpp:class:: Molecule : Copyable

    Represents a molecule.

    .. cpp:member:: std::vector<Atom *> atoms

        A list of atoms in the molecule.

    .. cpp:function:: Molecule()

        Empty constructor for forming arrays.

    .. cpp:function:: Molecule(const std::vector< Atom *> &atoms)

        Construct a molecule from a list of atoms.

    .. cpp:function:: int getsize() const

        Returns the number of atoms.

    .. cpp:function:: const std::vector< Atom *> &getatoms() const

        Returns the list of atoms in the molecule.

    .. cpp:function:: const Atom &getatom(int index) const

        Returns the atom at the given index.

    .. cpp:function:: Molecule &copy() override

        See :cpp:func:`Copyable::copy`.
