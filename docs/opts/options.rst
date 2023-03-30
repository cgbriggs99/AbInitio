Options
=======

The way options are passed to routines is through option list classes. Options are case-insensitive, but are stored in upper-case strings.

Option Lists
------------

.. cpp:namespace:: compchem

.. cpp:class:: OptionList

    Represents a list of options. 

    .. cpp:member:: protected std::map<std::string, bool> bool_opts

        Contains options with Boolean values.

    .. cpp:member:: protected std::map<std::string, int> int_opts

        Contains options with integer values.

    .. cpp:member:: protected std::map<std::string, double> double_opts

        Contains options with floating-point values.

    .. cpp:member:: protected std::map<std::string, std::string> string_opts

        Contains options with string values.

    .. cpp:function:: protected static std::string to_upper(const std::string &str)
    .. cpp:function:: protected static std::string to_upper(const char *str)

        :param str: A string to convert to upper-case.
        :return: The input string converted to upper-case.

        Converts a string to upper-case using :cpp:func:`std::isupper` and :cpp:func:`std::toupper` on each character. Thus, only characters affected by this are changed. Notable characters not changed include umlauts (ä, ö, ü), and presumably others. Polymorphic to accept C-type strings and C++ type strings.

    .. cpp:function:: OptionList()

        Constructs an empty list of options.

    .. cpp:function:: OptionList(const OptionList &copy)

        Creates a copy of an option list.

    .. cpp:function:: virtual ~OptionList()

        Clears and deletes a list of options.

    .. cpp:function:: int isoption(const std::string &str) const
    .. cpp:function:: int isoption(const char *str) const

        :param str: The option to check.
        :return: A value that determines whether the string represents an option, and if so, which kind of option.

        +--------------+--------------------+
        | Return Value | Type               |
        +==============+====================+
        | 0            | Not an option      |
        +--------------+--------------------+
        | 1            | :cpp:expr:`bool`   |
        +--------------+--------------------+
        | 2            | :cpp:expr:`int`    |
        +--------------+--------------------+
        | 3            | :cpp:expr:`double` |
        +--------------+--------------------+
        | 4            | :cpp:expr:`string` |
        +--------------+--------------------+

    .. cpp:function:: bool isoptionbool(const std::string &str) const
    .. cpp:function:: bool isoptionbool(const char *str) const

        :param str: The option to check
        :return: Whether the option is a Boolean option. Returns :cpp:expr:`true` if the option is a Boolean option. Returns :cpp:expr:`false` if the option is not a Boolean option, or the option string does not represent an available option.

    .. cpp:function:: bool isoptionint(const std::string &str) const
    .. cpp:function:: bool isoptionint(const char *str) const

        :param str: The option to check
        :return: Whether the option is an integer option. Returns :cpp:expr:`true` if the option is an integer option. Returns :cpp:expr:`false` if the option is not an integer option, or the option string does not represent an available option.

    .. cpp:function:: bool isoptiondouble(const std::string &str) const
    .. cpp:function:: bool isoptiondouble(const char *str) const

        :param str: The option to check
        :return: Whether the option is a floating-point option. Returns :cpp:expr:`true` if the option is a floating-point option. Returns :cpp:expr:`false` if the option is not a floating-point option, or the option string does not represent an available option.

    .. cpp:function:: bool isoptionstring(const std::string &str) const
    .. cpp:function:: bool isoptionstring(const char *str) const

        :param str: The option to check
        :return: Whether the option is a string option. Returns :cpp:expr:`true` if the option is a string option. Returns :cpp:expr:`false` if the option is not a string option, or the option string does not represent an available option.

    .. cpp:function:: bool getbooloption(const std::string &str) const
    .. cpp:function:: bool getbooloption(const char *str) const

        :param str: The option to get.
        :return: Returns the value of the option.
        :raises out_of_range: Throws an :cpp:expr:`std::out_of_range` exception if the string does not represent a Boolean option.

    .. cpp:function:: int getintoption(const std::string &str) const
    .. cpp:function:: int getintoption(const char *str) const

        :param str: The option to get.
        :return: Returns the value of the option.
        :raises out_of_range: Throws an :cpp:expr:`std::out_of_range` exception if the string does not represent an int option.

    .. cpp:function:: double getdoubleoption(const std::string &str) const
    .. cpp:function:: double getdoubleoption(const char *str) const

        :param str: The option to get.
        :return: Returns the value of the option.
        :raises out_of_range: Throws an :cpp:expr:`std::out_of_range` exception if the string does not represent a floating-point option.

    .. cpp:function:: const std::string &getstringoption(const std::string &str) const
    .. cpp:function:: const std::string &getstringoption(const char *str) const

        :param str: The option to get.
        :return: Returns the value of the option.
        :raises out_of_range: Throws an :cpp:expr:`std::out_of_range` exception if the string does not represent a string option.

    .. cpp:function:: void setbooloption(const std::string &str, bool value)
    .. cpp:function:: void setbooloption(const char *str, bool value)

        :param str: The option to add or set.
        :param value: The new value for the option.

        Sets the value of a Boolean option. If the option does not exist, it is added.

    .. cpp:function:: void setintoption(const std::string &str, int value)
    .. cpp:function:: void setintoption(const char *str, int value)

        :param str: The option to add or set.
        :param value: The new value for the option.

        Sets the value of an integer option. If the option does not exist, it is added.

    .. cpp:function:: void setdoubleoption(const std::string &str, double value)
    .. cpp:function:: void setdoubleoption(const char *str, double value)

        :param str: The option to add or set.
        :param value: The new value for the option.

        Sets the value of a floating-point option. If the option does not exist, it is added.

    .. cpp:function:: void setstringoption(const std::string &str, const std::string &value)
    .. cpp:function:: void setstringoption(const char *str, const std::string &value)
    .. cpp:function:: void setstringoption(const std::string &str, const char *value)
    .. cpp:function:: void setstringoption(const char *str, const char *value)

        :param str: The option to add or set.
        :param value: The new value for the option.

        Sets the value of a Boolean option. If the option does not exist, it is added.

    .. cpp:function:: virtual bool isglobal() const

        :return: Whether the options are the global options. For the base class, returns false.

Global Options
--------------

Some options can be made available to all parts of the program. These are stored in the :cpp:class:`compchem::GlobalOptions` class. They can be accessed by calling :cpp:func:`compchem::GlobalOptions::getsingleton`.

.. cpp:class:: GlobalOptions :public OptionList

    Stores options available to all parts of the program. Changes from one part should change all parts.

    .. cpp:function:: private GlobalOptions()

        Constructs the global options dictionary and sets it to be empty.

    .. cpp:function:: private GlobalOptions(const GlobalOptions &copy) = delete
    .. cpp:function:: private GlobalOptions(const GlobalOptions &&move) = delete

        Remove move and copy constructors, as there should only be one version of these options. Copies should be made using the copy constructor from :cpp:class:`compchem::OptionList`.

    .. cpp:function:: private ~GlobalOptions() = default

        Default destructor for the global options. It is made private so that the options don't get deleted by accident elsewhere.

    .. cpp:function:: static GlobalOptions &getsingleton()

        :return: A reference to the unique instance of this class. It is initialized with the default values by :cpp:func:`compchem::DefaultOptionsFactory::initializeoptions`.

    .. cpp:function:: bool isglobal() const override

        :return: This returns true, since the reference is the global option list.
