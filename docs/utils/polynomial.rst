Polynomials
===========

Definitions for polynomials.

.. cpp:namespace:: compchem

.. cpp:class:: template<int n> Polynomial

    :tparam n: The number of variables in the polynomial.

    Represents a polynomial with :code:`n` independent variables.

    .. cpp:member:: std::vector<std::array<int, n> > pows

        Holds the exponents for each term.

    .. cpp:member:: std::vector<double> coefs

        Holds the coefficients for each term.

    .. cpp:function:: explicit Polynomial(double scalar)

        :param scalar: The scalar to create from.

        Creates a polynomial that represents a scalar.

    .. cpp:function:: Polynomial(const std::vector<std::array<int, n> > &pows, const std::vector<double> &coefs)

        :param pows: The array of powers for each term.
        :param coefs: The array of coefficients for each term.

        Creates a polynomial with the given specification.

    .. cpp:function:: Polynomial(const Polynomial<n> &copy)

        :param copy: The polynomial to copy.

        Creats a polynomial that is a copy of the given polynomial.

    .. cpp:function:: Polynomial()

        Creates an empty polynomial.

    .. cpp:function:: ~Polynomial()

        Deconstructor.

    .. cpp:function:: const std::array<int, n> &gettermorder(int index) const

        :param index: The index of the term.
        :return:  An array containing the variable exponents of the requested term.

    .. cpp:function:: double getcoef(int index) const

        :param index: The index of the term.
        :return: The coefficient of the term.

    .. cpp:function:: int getsize() const

        :return: The number of terms in the polynomial.

    .. cpp:function:: protected void reduce()

        Gets rid of empty terms.

    .. cpp:function:: Polynomial<n> &operator+=(double rh)

        :param rh: The right-hand side of the operator.
        :return: A reference to the sum.

        Adds a constant to the polynomial.

    .. cpp:function:: Polynomial<n> &operator-=(double rh)

        :param rh: The right-hand side of the operator.
        :return: A reference to the difference.

        Subtracts a constant from the polynomial.

    .. cpp:function:: Polynomial<n> &operator*=(double rh)

        :param rh: The right-hand side of the operator.
        :return: A reference to the product.

        Multiplies the polynomial by a constant.

    .. cpp:function:: Polynomial<n> &operator/=(double rh)

        :param rh: The right-hand side of the operator.
        :return: A reference to the quotient.

        Divides the polynomial by a constant.

    .. cpp:function:: Polynomial<n> &operator+=(const Polynomial<n> &rh)

        :param rh: The right-hand side of the operator.
        :return: A reference to the sum.

        Adds a polynomial to this polynomial.

    .. cpp:function:: Polynomial<n> &operator-=(const Polynomial<n> &rh)

        :param rh: The right-hand side of the operator.
        :return: A reference to the difference.

        Subtracts a polynomial from this polynomial.

    .. cpp:function:: Polynomial<n> &operator*=(const Polynomial<n> &rh)

        :param rh: The right-hand side of the operator.
        :return: A reference to the product.

        Multiplies this polynomial by another polynomial.

    .. cpp:function:: Polynomial<n> &operator-() const

        :return: The negative of the polynomial.

        Takes the additive inverse of the polynomial.

    .. cpp:function:: Polynomial<n> &operator=(double d)

        :param d: The constant to set the polynomial to.
        :return: A reference to the polynomial.

        Sets the polynomial to the constant provided. The polynomial-polynomial version is not needed, as default behavior is assumed to be sufficient.

    .. cpp:function:: bool operator==(double d) const

        :param d: The constant to compare.
        :return: Whether the polynomial equals the constant.

        This checks whether the polynomial has no terms, or only one term corresponding to a constant term. Then, it checks this constant term against the constant.

    .. cpp:function:: bool operator!=(double d) const

        :param d: The constant to compare.
        :return: Whether the polynomial is different from the constant.

        This is the opposite of :cpp:func:`bool Polynomial::operator==(double)`.

    .. cpp:function:: bool operator==(const Polynomial<n> &rh) const

        :param rh: The other polynomial to compare.
        :return: Whether the two polynomials are equal.

    .. cpp:function:: bool operator!=(const Polynomial<n> &rh) const

        :param rh: The other polynomial to compare.
        :return: Whether the two polynomials are not equal.

    .. cpp:function:: double eval(double x, ...) const

        :param x...: The coordinates to evaluate the polynomial at.
        :return: The value of the polynomial at the given position.

    .. cpp:function:: Polynomial<n> &translate(double x, ...)

        :param x...: The new center coordinates for the polynomial.
        :return: The shifted polynomial.

        This function shifts a polynomial by the given amount in the directions provided. It also returns a reference to the polynomial so that it can be chained with other calculations.


Non Member Functions
--------------------

.. cpp:function:: template<int m> Polynomial<m> &operator+(double lh, const Polynomial<m> &rh)

    :param lh: The left-hand side of the operation.
    :param rh: The right-hand side of the operation.
    :return: The sum of the arguments.

.. cpp:function:: template<int m> Polynomial<m> &operator+(Polynomial<m> lh, double rh)

    :param lh: The left-hand side of the operation.
    :param rh: The right-hand side of the operation.
    :return: The sum of the arguments.

.. cpp:function:: template<int m> Polynomial<m> &operator+(Polynomial<m> lh, const Polynomial<m> &rh)

    :param lh: The left-hand side of the operation.
    :param rh: The right-hand side of the operation.
    :return: The sum of the arguments.

.. cpp:function:: template<int m> Polynomial<m> &operator-(double lh, const Polynomial<m> &rh)

    :param lh: The left-hand side of the operation.
    :param rh: The right-hand side of the operation.
    :return: The difference of the arguments.

.. cpp:function:: template<int m> Polynomial<m> &operator-(Polynomial<m> lh, double rh)

    :param lh: The left-hand side of the operation.
    :param rh: The right-hand side of the operation.
    :return: The difference of the arguments.

.. cpp:function:: template<int m> Polynomial<m> &operator-(Polynomial<m> lh, const Polynomial<m> &rh)

    :param lh: The left-hand side of the operation.
    :param rh: The right-hand side of the operation.
    :return: The difference of the arguments.

.. cpp:function:: template<int m> Polynomial<m> &operator*(double lh, const Polynomial<m> &rh)

    :param lh: The left-hand side of the operation.
    :param rh: The right-hand side of the operation.
    :return: The product of the arguments.

.. cpp:function:: template<int m> Polynomial<m> &operator*(Polynomial<m> lh, double rh)

    :param lh: The left-hand side of the operation.
    :param rh: The right-hand side of the operation.
    :return: The product of the arguments.

.. cpp:function:: template<int m> Polynomial<m> &operator*(Polynomial<m> lh, const Polynomial<m> &rh)

    :param lh: The left-hand side of the operation.
    :param rh: The right-hand side of the operation.
    :return: The product of the arguments.

.. cpp:function:: template<int m> Polynomial<m> &operator/(Polynomial<m> lh, double rh)

    :param lh: The left-hand side of the operation. Only a polynomial.
    :param rh: The right-hand side of the operation. Only a constant.
    :return: The quotient of the arguments.

.. cpp:function:: template<int m> bool operator==(double lh, const Polynomial<m> &rh)

    :param lh: The left-hand side of the operation.
    :param rh: The right-hand side of the operation.
    :return: Whether the two are equal. See :cpp:func:`bool Polynomial::operator==(double)`.

.. cpp:function:: template<int m> Polynomial<m> &pow(const Polynomial<m> &mant, int expo)

    :param mant: The mantissa, a polynomial.
    :param expo: The exponent. Can only be a non-negative integer.
    :return: The polynomial raised to a power.


.. cpp:function:: Polynomial<3> &sphereharm(int l, int ml)

    :param l: The angular momentum quantum number.
    :param ml: The magnetic quantum number.
    :return: A polynomial representing the real spherical harmonics in cartesian coordinates.
