copyable.hpp
============

This contains a virtual class which guarantees a copy method for its children.

.. cpp:namespace:: compchem

.. cpp:class:: Copyable

    Represents a copyable object.

    .. cpp:function:: virtual Copyable &copy()

        Return a copy of the object.
