#ifndef UTILS_H
#define UTILS_H

struct Point                    // structure to hold a vertex
{
    float  x;
    float  y;
    float  z;
};

struct Derivative
{
    float     x;                  //  the first partial deriv with respect to x  (x-gradient)
    float     y;                  //  the first partial deriv with respect to y  (y-gradient)
    float     z;                  //  the first partial deriv with respect to z  (z-gradient)
    float     xy;                 //  the second partial deriv with respect to xy
    float     xz;                 //  the second partial deriv with respect to xz
    float     yz;                 //  the second partial deriv with respect to yz
    float     xx;                 //  the second partial deriv with respect to xx
    float     yy;                 //  the second partial deriv with respect to yy
    float     zz;                 //  the second partial deriv with respect to zz
};


struct Curvature
{
    float     gaussian;           //  the gaussian curvature for a point  (Kg)
    float     mean;               //  the mean curvature for a point      (Ka)
    float     min;                //  the minimun curvature for a point   (Km)
    float     max;                //  the maximum curvature for a point   (KM)
};


#endif  // UTILS_H


