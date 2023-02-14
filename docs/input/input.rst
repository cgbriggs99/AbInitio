Molecular Input
===============

Contains various ways to input molecules.

.. cpp:namespace:: compchem

.. cpp:function:: Molecule &parseXYZ(std::FILE *fp)

    :param fp: The .xyz file to parse. Only a certain subset can be handled at this time.
    :return: The molecule represented by the .xyz file.
