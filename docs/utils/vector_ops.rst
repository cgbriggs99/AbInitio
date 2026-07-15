Vector Operations
=================

This file contains a few vector operations.

.. cpp:namespace:: compchem

.. cpp:function:: double rmsnorm(const std::vector<double> &v1, const std::vector<double> &v2)

    :params v1, v2: The vectors to norm.
    :return: The root-mean-squared norm of the difference of two vectors.

    Computes :math:`\sqrt{\frac{||v_1 - v_2||^2}{n}}` for two :math:`n` dimensional vectors.

.. cpp:function:: double dotprod(const std::vector<double> &v1, const std::vector<double> &v2)

    :params v1, v2: The two vectors to take the dot product.
    :return: The dot product of the two vectors.

    Computes :math:`v_1 \cdot v_2`.

.. cpp:function:: double vecmag(const std::vector<double> &vec)

    :param vec: The vector to find the magnitude of.
    :return: The magnitude of the vector.

    Computes :math:`||v||`.
