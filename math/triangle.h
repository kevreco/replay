#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "vector3.h"

struct Triangle
{

    Vector3 a,b,c;

    Triangle() {}
    Triangle(const Triangle& t) : a(t.a), b(t.b), c(t.c) {}
    Triangle(const Vector3& a, const Vector3& b, const Vector3& c) : a(a), b(b), c(c) {}

}; // Triangle


struct Tetrahedron
{
    Triangle abc, acd, adb , dbc;

    Tetrahedron(const Vector3& a, const Vector3& b, const Vector3& c, const Vector3& d)
    {
        abc = Triangle(a, b, c);
        acd = Triangle(a, c, d);
        adb = Triangle(a, d, b);
        dbc = Triangle(d, b, c);
    }

}; // Tetrahedron

#endif // TRIANGLE_H
