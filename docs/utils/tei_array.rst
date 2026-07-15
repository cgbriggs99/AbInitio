Two-Electron Integral Arrays
============================

This file contains a storage class for the two-electron integrals that respects the symmetry in the arguments.

.. cpp:namespace:: compchem

.. cpp:class:: TEIArray

    Represents a collection of two electron integrals. This collection obeys the following symmetry.

    .. math::

        \left(\mu\nu\middle|\lambda\sigma\right) = \left(\nu\mu\middle|\lambda\sigma\right) = \left(\mu\nu\middle|\sigma\lambda\right) = \left(\nu\mu\middle|\sigma\lambda\right) = \left(\lambda\sigma\middle|\mu\nu\right) = \left(\sigma\lambda\middle|\mu\nu\right) = \left(\lambda\sigma\middle|\nu\mu\right) = \left(\sigma\lambda\middle|\nu\mu\right)

    .. cpp:member:: double *data

        The internal data buffer.

    .. cpp:member:: int dim

        The number of orbitals represented.

    .. cpp:function:: TEIArray(double *data, int dim)

        :param data: A pointer to a data array. Ownership is taken over by the class, and all memory management, including freeing is handled. It must be restricted.
        :param dim: The number of orbitals.

        Constructor that takes over an already set up data buffer.

    .. cpp:function:: TEIArray(int dim)

        :param dim: The number of orbitals.

        Creates a new empty array.

    .. cpp:function:: TEIArray(const compchem::TEIArray &copy)

        :param copy: The instance to copy.

        Copies the data from one array to another.

    .. cpp:function:: ~TEIArray()

        Destructor.

    .. cpp:function:: double at(int mu, int nu, int lam, int sig) const
    .. cpp:function:: double &at(int mu, int nu, int lam, int sig)
    .. cpp:function:: double operator()(int mu, int nu, int lam, int sig) const
    .. cpp:function:: double &operator()(int mu, int nu, int lam, int sig)

        :params mu, nu, lam, sig: The indices of the integral to query or set.
        :return: The value of the integral at the given point, or a reference to this value for manipulation.

        Returns a value stored in the array, or a reference to a place in the array.

    .. cpp:function:: double at_direct(int index) const
    .. cpp:function:: double &at_direct(int index)

        :param index: The index to the data buffer.
        :return: Thhe value of the integral, or a reference to this value.

        Uses a direct index to the data buffer, as opposed to the four-part index in the :cpp:func:`at` functions.

    .. cpp:function::  const double *getdata() const

        :return: A pointer to the data buffer.

        Returns the data buffer.

    .. cpp:function:: int getdim() const

        :return: The number of orbitals represented in the array.

    .. cpp:function:: int getindex(int mu, int nu, int lam, int sig) const

        :params mu, nu, lam, sig: The four-part indices.
        :return: A direct index to the array.

        To respect the symmetry of the array, the following procedure is used. If :math:`\nu > \mu`, swap them. If :math:`\sigma > \lambda`, swap them.  This should ensure that :math:`\mu \ge \nu` and :math:`\lambda \ge \sigma`. Then, compute the secondary values, :math:`\mu\nu = \frac{\mu\left(\mu + 1\right)}{2} + \nu` and :math:`\lambda\sigma = \frac{\lambda\left(\lambda + 1\right)}{2} + \sigma`. Then, if :math:`\lambda\sigma > \mu\nu`, swap them. This should ensure that :math:`\mu\nu \ge \lambda\sigma`. Finally, the index can be calculated as :math:`\mu\nu\lambda\sigma = \frac{\mu\nu\left(\mu\nu + 1\right)}{2} + \lambda\sigma`.

    .. cpp:function:: int getsize() const

        :return: The number of slots in the array.

        To find this, let :math:`n` be the result of :cpp:func:`getdim`. Then, let :math:`s = \frac{n(n - 1)}{2} + n - 1` be an intermediate. Then the total size will be :math:`N = \frac{s (s + 1)}{2} + s + 1`.

    .. cpp:function:: void indextoquad(int index, int *mu, int *nu, int *lam, int *sig) const

        :param index: The direct index to separate.
        :param mu, nu, lam, sig: The output values for the four-part index.

        Takes a direct index ans splits it into a four-part index. This is done by using :cpp:func:`biggest_trinum_leq` and :cpp:func:`biggest_triind_leq`.

Non-member Functions
--------------------

.. cpp:function:: int biggest_trinum_leq(int val)

    :param val: The bound for the triangular number.
    :return: The largest triangular number less than or equal to `val`.

    This function works by first using :cpp:func:`biggest_triind_leq`, then passing the value to :cpp:func:`triangular_num`.

.. cpp:function:: int biggest_triind_leq(int val)

    :param val: The bound for the triangular number.
    :return: The index of the triangular number.

    Finds the index of the largest triangular number less than or equal to the value. That is, it finds the largest :math:`n` that satisfies :math:`T_n \le x`. This can be found using the formula :math:`n = \left\lfloor\frac{\sqrt{1 + 8x} - 1}{2}\right\rfloor`.

.. cpp:function:: int triangular_num(int index)

    :param index: The index of the triangular number.
    :return: The triangular number with the given index.

    To find the :math:`n\mathrm{th}` triangular number, the following formulae may be used.

    .. math::

        T_n = \sum_{k = 0}^n k = \frac{n (n + 1)}{2}
