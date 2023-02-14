molecule.hpp
============

Contains a class to represent a molecule.

.. cpp:namespace:: compchem

.. cpp:class:: Molecule

    Represents a molecule.

    .. cpp:member:: std::vector<Atom *> atoms

        A list of atoms in the molecule.

    .. cpp:member:: double comx
    .. cpp:member:: double comy
    .. cpp:member:: double comz

        The components of the center of mass.

    .. cpp:function:: Molecule()

        Empty constructor for forming arrays.

    .. cpp:function:: Molecule(const std::vector< Atom *> &atoms)
        
        Construct a molecule from a list of atoms.

    .. cpp:function:: Molecule(const Molecule &copy)
        
        Copy constructor.

    .. cpp:function:: ~Molecule()

        Destructor for the molecule.

    .. cpp:function:: int getsize() const
        
        Returns the number of atoms.

    .. cpp:function:: const std::vector< Atom *> &getatoms() const
        
        Returns the list of atoms in the molecule.

    .. cpp:function:: const Atom &getatom(int index) const
        
        Returns the atom at the given index.

    .. cpp:function:: void addatom(Atom &atom)
        
        Adds the specified atom to the molecule.

    .. cpp:function:: Molecule *copy() const

        Returns a copy of this molecule.
        
    .. cpp:function:: private void calc_com()

        Calculates the center of mass and sets the appropriate fields.

    .. cpp:function:: getcomx() const
    .. cpp:function:: getcomy() const
    .. cpp:function:: getcomz() const

        Returns the requested component of the center of mass.

    .. cpp:function:: void translate_to_com()

        Translates the molecule to its center of mass.
