#ifndef SPHERE_H
#define SPHERE_H

#include "vector3.h"

struct Sphere
{

    Vector3 center;
    float radius;


    inline Sphere() : center(), radius(1.0) {}
    inline Sphere(const Sphere& sphere) : center(sphere.center), radius(sphere.radius) {}
    inline Sphere(const Vector3& center, float radius = 1.0f) : center(center), radius(radius) {}

    void print() const {
        printf("[c:%g,%g,%g;r:%g]", center.x, center.y, center.z, radius);
    }
};

#endif // SPHERE_H
