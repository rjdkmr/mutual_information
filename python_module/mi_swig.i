%module mi

%{
    #include "../src/Mutual_information.hpp"
%}

%include "../src/Mutual_information.hpp"
%include "carrays.i"

%array_functions(double, doubleArray);
%{
    double mi(double *x, double *y, int n, int k);
%}
