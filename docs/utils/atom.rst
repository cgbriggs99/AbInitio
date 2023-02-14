 
atom.hpp
========

This header file contains atom definitions.

.. cpp:namespace:: compchem

.. cpp:function:: int getZFromSymb(const std::string &symb)

    Returns the atomic number of a given symbol. If passed "D" or "T", it will return `1`, as these refer to deuterium and tritium. If the string is not an official symbol, it will return `0`.

.. cpp:function:: double getAbundantMass(int Z)

    Returns the mass of the most abundant isotope of the element with the given atomic number.

.. cpp:class:: Atom

    Represents an atom.

    .. cpp:member:: int Z

        The atomic number of the atom

    .. cpp:member:: int charge

        The total charge on the atom.

    .. cpp:member:: double x
    .. cpp:member:: double y
    .. cpp:member:: double z

        The position of the atom in space in Bohr radii.

    .. cpp:member:: double mass

        The mass of the nucleus in atomic units.

    .. cpp:member:: int norbitals

        The number of orbitals attached to the atom.

    .. cpp:member:: std::vector<BasisOrbital *> orbitals

        A list of orbitals.

    .. cpp:function:: Atom()

        An empty constructor to use when making arrays.
	
    .. cpp:function:: Atom(int Z, double x, double y, double z)
        
        Creates an atom with the specified atomic number and position. The mass will be set to the mass listed on the periodic table, and the orbitals will be initialized to the null pointer. The charge is assumed to be zero.

    .. cpp:function:: Atom(int Z, double mass, double x, double y, double z)
        
        Creates an atom with the specified atomic number, mass, and position. The orbitals will be initialized to the null pointer. The charge is assumed to be zero.

    .. cpp:function:: Atom(int Z, double x, double y, double z, int norbs, const BasisOrbital *orbitals)
        
        Creates an atom with the specified atomic number, position, and list of orbitals. The orbitals will be copied. The mass will be initialized to the mass listed on the periodic table. The charge is assumed to be 0.

    .. cpp:function:: Atom(int Z, double mass, double x, double y, double z, int norbs, const BasisOrbital *orbitals)
        
        Creates an atom with the specified atomic number, mass, position, and orbitals. The orbitals will be copied, and the charge assumed to be 0.

    .. cpp:function:: Atom(int Z, int charge, double x, double y, double z)
        
        Creates an atom with the specified atomic number, charge, and position. The mass will be set to the mass listed on the periodic table, and the orbitals will be initialized to the null pointer.

    .. cpp:function:: Atom(int Z, int charge, double mass, double x, double y, double z)
        
        Creates an atom with the specified atomic number, charge, mass, and position. The orbitals will be initialized to the null pointer.

    .. cpp:function:: Atom(int Z, int charge, double x, double y, double z, int norbs, const BasisOrbital *orbitals)
        
        Creates an atom with the specified atomic number, charge, position, and list of orbitals. The orbitals will be copied. The mass will be initialized to the mass listed on the periodic table.

    .. cpp:function:: Atom(int Z, int charge, double mass, double x, double y, double z, int norbs, const BasisOrbital *orbitals)
        
        Creates a fully specified atom. The orbitals will be copied.

    .. cpp:function:: ~Atom()

        Deconstructor. It is virtual to allow for extending of this class.


    .. cpp:function:: double getx() const
    .. cpp:function:: double gety() const
    .. cpp:function:: double getz() const
		      
        Returns the x, y, or z coordinate of the atom.

    .. cpp:function:: void setx(double x)
    .. cpp:function:: void sety(double y)
    .. cpp:function:: void setz(double z)
        
        Sets the x, y, or z coordinate of the atom.


    .. cpp:function:: double getmass() const
        
        Returns the atomic mass of the nucleus.

    .. cpp:function:: int getZ() const
        
        Returns the atomic number of the atom.

    .. cpp:function:: int getcharge() const
        
        Returns the charge of an atom.

    .. cpp:function:: int getnorbitals() const
        
        Returns the number of orbitals on an atom.

    .. cpp:function:: const std::vector<BasisOrbital *> &getOrbitals() const
        
        Returns the orbitals.

    .. cpp:function:: const BasisOrbital &getorbital(int index) const
        
        Returns the basis function at the given index.

    .. cpp:function:: Atom *copy() const
		      
        Returns a copy of this atom.
        
