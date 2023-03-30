Default Options
===============

These are the default options set by the program.

Boolean Options
---------------

.. csv-table::
    :header: "Name", "Value", "Reference"

    "``ANALYTIC BOYS``", "``true``", ":cpp:func:`compchem::AnalyticIntegral::boys_square`"


Integer Options
---------------

.. csv-table::
    :header: "Name", "Value", "Reference"

    "``THREADS``", "``4``", ":cpp:func:`compchem::IntegralFactory::Smatrix`, :cpp:func:`compchem::IntegralFactory::Tmatrix`, :cpp:func:`compchem::IntegralFactory::Vmatrix`, :cpp:func:`compchem::IntegralFactory::TEIints`"
    "``BOYS POINTS``", "``32``", ":cpp:func:`compchem::AnalyticIntegral::boys_square`"
    "``MAX SCF CYCLES``", "``100``", ""

Floating-point Options
----------------------

.. csv-table::
    :header: "Name", "Value", "Reference"

    "``SCF RMS CONVERGENCE``", "``1e-6``", ""
    "``SCF ENERGY CONVERGENCE``", "``1e-7``", ""

String Options
--------------


Initializer
-----------

This is the class that initializes the default options.


.. cpp:namespace:: compchem

.. cpp:class:: DefaultOptionsFactory

    .. cpp:function:: static void initializeoptions(OptionList &opts)

        :param opts: The option list to initialize with the defaults.

        Initializes the option list with the default options specified above.
