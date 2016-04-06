#include "plane.h"

#include "matrix4.h"

Plane::Plane() :
origin(),
d(0.0),
normal(0.0,1.0,0.0)
{
    d = (this->normal.dot(origin));
}

Plane::Plane(const Vector3& pos, const Vector3& normal) :
origin(pos),
d(0),
normal(normal)
{
    this->normal.normalize();
    d = (this->normal.dot(pos));
}

Plane::Plane(const Vector3& v1, const Vector3& v2, const Vector3& v3)
{
    normal = (v1 - v2).normalized().cross( (v1 - v3).normalized()) ;
    normal.normalize();
    //origin = v2;
    //set()
    //d = (normal.dot(origin));
    setWithNormal(v2, normal);
}
 void Plane::setWithNormal(const Vector3& pos, const Vector3& normal)
 {
    origin = pos;
    this->normal.normalize();
    d = (this->normal.dot(pos));
 }

Plane& Plane::operator= (const Plane& p)
{
    if (this != &p)
    {
        this->normal = p.normal;
        this->d = p.d;
        this->origin = p.origin;
    }
    return *this;
}

void Plane::getEquation(double eq[4]) const
{
    eq[0] = normal.x;
    eq[1] = normal.y;
    eq[2] = normal.z;
    eq[3] = (origin.dot(normal));
}


void Plane::getEquation(float eq[4]) const
{
    eq[0] = normal.x;
    eq[1] = normal.y;
    eq[2] = normal.z;
    eq[3] = (origin.dot(normal));
}
