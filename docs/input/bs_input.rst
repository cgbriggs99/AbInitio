Basis Set Input
===============

Contains a reader for Psi4 .gbs basis set files.

.. cpp:namespace:: compchem

.. cpp:function:: std::vector<compchem::BasisOrbital *> *readPsi4file(std::FILE *fp, int Z, int charge = 0)

    :param fp: A file pointer to the opened .gbs file.
    :param Z: The atomic number to look for.
    :param charge: The charge of the atom, which may change the basis set.
    :return: A vector containing pointers to basis orbitals.
    :raise runtime_error: Raises runtime errors when a line is malformed (may be changed to a syntax error at some point in the future), or when a basis set is not found for the given atom.

    This function reads a Psi4 .gbs file and returns the appropriate basis functinos. This function rewinds the file pointer before execution, but not after.

GBS File Format
---------------

The format of a .gbs file is discussed here. A file must start with a kind specification, followed by data. Anything after a :code:`!` is considered a comment, and is ignored.

Kind Specification
^^^^^^^^^^^^^^^^^^

At the beginning of a .gbs file, there should be a line containing only the word `spherical` or `cartesian`. This determines whether a spherical basis set or a cartesian basis set is used. At the moment, this is ignored, and a spherical basis set is assumed.

Basis Set Specification
^^^^^^^^^^^^^^^^^^^^^^^

After specifying the kind of basis set this file represents, the parameters for the basis set are specified. This starts with a line of asterisks, followed by the element's symbol and (what is presumed to be) the charge on one line. After that, the different orbitals are specified. First, a line containing the type of orbital, the number of Gaussian terms, and the normalization factor appears. After this, the parameters show up as an alpha-coefficient pair, or for SP orbitals, an alpha, s-coefficient, p-coefficient triplet. An example taken from Psi4's STO-3G basis set is shown below.::


    ****
    C 0
    S 3 1.00
         71.6168370 0.15432897
         13.0450960 0.53532814
          3.5305122 0.44463454
    SP 3 1.00
          2.9412494 -0.09996723 0.15591627
          0.6834831 0.39951283 0.60768372
          0.2222899 0.70011547 0.39195739

This represents a neutral carbon atom with three terms for each of the orbitals. The middle column for the SP orbital contains the S coefficients, while the third column contains the P coefficients. The alphas are the same for both.


