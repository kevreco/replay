#ifndef PLANE_H
#define PLANE_H

#include "vector3.h"
#include "maths.h"

struct Plane2
{
    Vector2 origin;
    Vector2 normal;
    float dist;

    Plane2(const Vector2& pos, const Vector2& normal) :
        origin(pos),
        normal(normal),
        dist(0)
    {
        this->normal.normalize();
        dist = (this->normal.dotProduct(pos));
    }

    static Vector2 ProjectedPoint(const Plane2& p, const Vector2& point)
    {
        double dist = point.x * p.normal.x + point.y * p.normal.y - p.dist;

        return (point - (p.normal * dist));
    }
}; // struct Plane2

struct Plane
{
    Vector3 origin; //a b c
    double d;       //d
    Vector3 normal; //x y z

    Plane();
    Plane(const Vector3& pos, const Vector3& normal);
    Plane(const Vector3& v1, const Vector3& v2, const Vector3& v3);

    void setWithNormal(const Vector3& pos, const Vector3& normal);

    Plane& operator =(const Plane& p);

    inline bool operator==(const Plane& plane) const
    {
        return (origin == plane.origin && d == plane.d);
    }

    void getEquation(double eq[4]) const;
    void getEquation(float eq[4]) const;

    /*inline friend std::ostream& operator<<(std::ostream& os, const Plane& p)
    {
        os << "Plane[o:"<<p.origin<<" n:"<<p.normal<<" d:"<<p.d<<"]";
        return os;
    }
    */

    void print()
    {
        printf("Plane(%.2f, %.2f, %.2f; %.2f, %.2f, %.2f; %.2f)",
               origin.x, origin.y, origin.z,
               normal.x, normal.y, normal.z,
               d);
    }

    Vector3 project(const Vector3& point) const
    {
        double dist = point.x*normal.x + point.y*normal.y + point.z*normal.z - d;

        Vector3  projected = point - (normal * dist);
        return projected;
    }

    static Vector3 ProjectedPoint(const Plane& p, const Vector3& point)
    {
        double dist = point.x * p.normal.x + point.y * p.normal.y + point.z * p.normal.z - p.d;

        return (point - (p.normal * dist));
    }


}; // struct Plane

#endif // PLANE_H
