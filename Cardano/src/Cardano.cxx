// cubic equation solver example using Cardano's method
#include "itkImage.h"
#include <iostream>
#include <utils.h>
#include <cmath>
#define THIRD 0.333333333333333
#define ROOTTHREE 1.73205080756888     

using namespace std;

// this function returns the cube root if x were a negative number aswell
double cubeRoot(double x)
{
    if (x < 0)
        return -pow(-x, THIRD);
    else
        return pow(x, THIRD);
}

void solveCubic(double a, double b, double c, double d)
{
    // find the discriminant
    double f, g, h;
    f = (3 * c / a - pow(b, 2) / pow(a, 2)) / 3;
    g = (2 * pow(b, 3) / pow(a, 3) - 9 * b * c / pow(a, 2) + 27 * d / a) / 27;
    h = pow(g, 2) / 4 + pow(f, 3) / 27;
    // evaluate discriminant
    if (f == 0 && g == 0 && h == 0)
    {
        // 3 equal roots
        double x;
        // when f, g, and h all equal 0 the roots can be found by the following line
        x = -cubeRoot(d / a);
        // print solutions
        cout
            << "x = " << endl
            << " " << x << endl
            << " " << x << endl
            << " " << x << endl << endl;
    }
    else if (h <= 0)
    {
        // 3 real roots
        double q, i, j, k, l, m, n, p;
        // complicated maths making use of the method
        i = pow(pow(g, 2) / 4 - h, 0.5);
        j = cubeRoot(i);
        k = acos(-(g / (2 * i)));
        m = cos(k / 3);
        n = ROOTTHREE * sin(k / 3);
        p = -(b / (3 * a));
        // print solutions
        cout
            << "x = " << endl
            << " " << 2 * j * m + p << endl
            << " " << -j * (m + n) + p << endl
            << " " << -j * (m - n) + p << endl << endl;
    }
    else if (h > 0)
    {
        // 1 real root and 2 complex roots
        double r, s, t, u, p;
        // complicated maths making use of the method
        r = -(g / 2) + pow(h, 0.5);
        s = cubeRoot(r);
        t = -(g / 2) - pow(h, 0.5);
        u = cubeRoot(t);
        p = -(b / (3 * a));
        // print solutions
        cout
            << "x = " << endl
            << " " << (s + u) + p << endl
            << " " << -(s + u) / 2 + p << " +" << (s - u) * ROOTTHREE / 2 << "i" << endl
            << " " << -(s + u) / 2 + p << " " << -(s - u) * ROOTTHREE / 2 << "i" << endl << endl;
    }
}


int main()
{
    double a, b, c, d;
    // introduction
    cout << "Cubic Equation Solver" << endl;
    cout << "ax^3 + bx^2 + cx + d = 0" << endl;
    cout << "with a, b, c, d are real, and a is not zero" << endl << endl;
    // infinite loop
    while (1)
    {
        // get the co-efficients of x
        cout << "a = ";
        cin >> a;
        cout << "b = ";
        cin >> b;
        cout << "c = ";
        cin >> c;
        cout << "d = ";
        cin >> d;
        solveCubic(a, b, c, d);
    }
    return 0;
}
