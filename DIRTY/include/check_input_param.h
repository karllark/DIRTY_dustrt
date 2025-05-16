#ifndef _DIRTY_CHECK_INPUT_PARAM_
#define _DIRTY_CHECK_INPUT_PARAM_
#include "compat.h"

//**********************************************************************
// template function definitions

// check for required numerical parameter with min/max bounds
template <class FI, class FI2, class FI3> void check_input_param(string const &param, FI x, FI2 min_x, FI3 max_x)
{
    // first check if the value is set
    if ((x == -99) || (!finite(x)))
    {
        cout << "Required parameter " << param << " not set." << endl;
        exit(8);
    }
    else if ((x < min_x) || (x > max_x))
    { // then check it is between the limits
        cout << "Parameter '" << param << "' (= " << x << ") is outside of allowed range = [";
        cout << min_x << "," << max_x << "]" << endl;
        exit(8);
    }
}

// set required string parameter to default value
template <class FI> void check_input_param(string const &param, FI x, string const &def_x)
{
    // if not set, then set to default and print warning
    if (x == "NULL")
    {
        cout << "Setting required parameter '" << param << "' to default value of '";
        cout << def_x << "'" << endl;
    }
}

//**********************************************************************
// external function definitions

#endif
