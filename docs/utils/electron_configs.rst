Electronic Configuration
========================

Represents the ground-state electronic configuration of an element.


.. cpp:namespace:: compchem

.. cpp:class:: GSConfig

    A general ground-state electronic configuration.

    .. cpp:member:: int *confs

        The occupations of the orbitals. These are in order of filling, so 1s, 2s, 2p, 3s, 3p, 4s, 3d, 4p, 5s, 4d, 5p, 6s, 4f, 5d, 6p, 7s, 5f, 6d, 7p.

    .. cpp:member:: int size_confs

        The number of shells.


    .. cpp:member:: int Z

        The atomic number of the represented element.

    .. cpp:function:: GSConfig(int z, std::initializer_list<int> conf)

        :param z: The atomic number of the element.
        :param conf: The occupations of the shells.

        Creates a ground-state configuration.

    .. cpp:function:: int getShells() const

        :return: The number of shells.

    .. cpp:function:: int getZ() const

        :return: The atomic number being represented.

    .. cpp:function:: int operator[](int index) const

        :param index: The index to look at.
        :return: The occupation of the requested shell.

Non Member Functions
--------------------

.. cpp:function:: const GSConfig &getconfig(int Z)

    :param Z: The atomic number to query.
    :return: The ground-state configuration of the requested element.
