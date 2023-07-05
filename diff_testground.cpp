// C++ includes
#include <iostream>
#include <ostream>
#include <vector>
#include <cmath>
#include <limits>
#include <numbers>
#include <iomanip>

// TNL includes
#include "TNL/Containers/Vector.h"

// autodiff includes
#include <autodiff/forward/dual.hpp>
using namespace autodiff;

using Vector = TNL::Containers::Vector< dual >;

const double EPSILON = 1e-10;

double poly( double x )
{
    return x*x*x-2*x*x+3*x-4;
}

double d_poly( double x )
{
    return 2*x-2;
}

int main()
{
    std::vector< double > val = {1.0, 2.0, 3.0, 10.0, 100.0 };
    // printing in format of a LaTeX table
    for( double x: val )
    {
        std::cout << x << " & ";
        double delta = ( poly( x + EPSILON ) - poly( x ) ) / EPSILON;
        std::cout << std::setprecision(6) << delta 
                  << " & " << d_poly( x ) << " \\\\";
    }
}
