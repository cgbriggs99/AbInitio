Integral Factory
================

This represents a class that produces matrices that can be used in energy calculations.

.. cpp:namespace:: compchem

.. cpp:class:: template<typename Ints> IntegralFactory

    A factory that can produce integral matrices.

    .. cpp:member:: protected compchem::OptionList &opts

        Options to use for computations.

    .. cpp:function:: IntegralFactory()

        Constructor that initializes the options with a copy of the global options.

    .. cpp:function:: IntegralFactory(const OptionList &opts)

        Sets up the factory with the given options.

    .. cpp:function:: virtual ~IntegralFactory()

        Deconstructor.

    .. cpp:function:: void Smatrix(const compchem::Molecule *mol, double *out, int *dim)

        :param mol: The molecule to use when calculating the integrals.
        :param out: The output array for the integral values. Must have at least as many elements as the square of the number of orbitals on the atom.
        :param dim: Output for the number of rows or columns in the array.

        Computes the overlap matrix. It can be multithreaded if the ``THREADS`` option is greater than 1.

        **Options**

        ``THREADS``: *default:* :cpp:expr:`4`

        The number of threads the program can use.

    .. cpp:function:: void Tmatrix(const compchem::Molecule *mol, double *out, int *dim)

        :param mol: The molecule to use when calculating the integrals.
        :param out: The output array for the integral values. Must have at least as many elements as the square of the number of orbitals on the atom.
        :param dim: Output for the number of rows or columns in the array.

        Computes the kinetic energy matrix. It can be multithreaded if the ``THREADS`` option is greater than 1.

        **Options**

        ``THREADS``: *default:* :cpp:expr:`4`

        The number of threads the program can use.

    .. cpp:function:: void Vmatrix(const compchem::Molecule *mol, double *out, int *dim)

        :param mol: The molecule to use when calculating the integrals.
        :param out: The output array for the integral values. Must have at least as many elements as the square of the number of orbitals on the atom.
        :param dim: Output for the number of rows or columns in the array.

        Computes the potential energy matrix. It can be multithreaded if the ``THREADS`` option is greater than 1.

        **Options**

        ``THREADS``: *default:* :cpp:expr:`4`

        The number of threads the program can use.

    .. cpp:function:: compchem::TEIArray *TEIints(const compchem::Molecule *mol)

        :param mol: The molecule for the computation of the two-electron repulsion integrals.
        :return: The repulsion integrals.

        Computes the two-electron repulsion integrals. It can be multithreaded if the ``THREADS`` option is greater than 1.

        **Options**

        ``THREADS``: *default:* :cpp:expr:`4`

        The number of threads the program can use.

Protected Methods
-----------------

    .. cpp:function:: protected static void s_routine(const std::vector<const compchem::GaussianOrbital *> *orbs, const std::vector<std::array<double, 3> *> *centers, double *out, int dim, int thread_num, int threads, compchem::OptionList &opts)

        :param orbs: List of orbitals.
        :param centers: List of centers.
        :param out: The output array.
        :param dim: The number of rows in the output array.
        :param thread_num: The number of the current thread in its pool.
        :param threads: The number of threads in the pool.
        :param opts: Options for the integrals.

        Worker function for :cpp:func:`compchem::IntegralFactory::Smatrix`.

    .. cpp:function:: protected static void t_routine(const std::vector<const compchem::GaussianOrbital *> *orbs, const std::vector<std::array<double, 3> *> *centers, double *out, int dim, int thread_num, int threads, compchem::OptionList &opts)

        :param orbs: List of orbitals.
        :param centers: List of centers.
        :param out: The output array.
        :param dim: The number of rows in the output array.
        :param thread_num: The number of the current thread in its pool.
        :param threads: The number of threads in the pool.
        :param opts: Options for the integrals.

        Worker function for :cpp:func:`compchem::IntegralFactory::Tmatrix`.

    .. cpp:function:: protected static void v_routine(const std::vector<const compchem::GaussianOrbital *> *orbs, const std::vector<std::array<double, 3> *> *centers, double *out, int dim, int thread_num, int threads, compchem::OptionList &opts, const compchem::Molecule *mol)

        :param orbs: List of orbitals.
        :param centers: List of centers.
        :param out: The output array.
        :param dim: The number of rows in the output array.
        :param thread_num: The number of the current thread in its pool.
        :param threads: The number of threads in the pool.
        :param opts: Options for the integrals.
        :param mol: The molecule for the computation.

        Worker function for :cpp:func:`compchem::IntegralFactory::Vmatrix`.

    .. cpp:function:: protected static void tei_routine(const std::vector<const compchem::GaussianOrbital *> *orbs, const std::vector<std::array<double, 3> *> *centers, double *out, int dim, int thread_num, int threads, compchem::OptionList &opts)

        :param orbs: List of orbitals.
        :param centers: List of centers.
        :param out: The output array. It is indexed like the internal indices for :cpp:class:`compchem::TEIArray`.
        :param dim: The number of rows in the output array.
        :param thread_num: The number of the current thread in its pool.
        :param threads: The number of threads in the pool.
        :param opts: Options for the integrals.

        Worker function for :cpp:func:`compchem::IntegralFactory::TEIints`.
