#ifndef RAY_H
#define RAY_H

#include "vector2.h"
#include "vector3.h"

// 2D
struct Ray2 {

    Vector2 origin;
    Vector2 direction;

    // inline
    static Ray2 FromDirection(const Vector2& p, const Vector2& dir) {
        Ray2 r;
        r.origin = p;
        r.direction = dir;
        return r;
    }

    // inline
    static Ray2 FromPoint(const Vector2& p1, const Vector2& p2) {
        return FromDirection(p1, p2 - p1);
    }

    inline Vector2 pointAt(float t) const { return origin + direction.normalized() * t; }

    inline bool isZero() const { return (origin == Vector2::Zero && direction == Vector2::Zero);}

    // inline
    void print() const
    {
        printf("Ray(%.2f, %.2f; %.2f, %.2f)",
               origin.x, origin.y,
               direction.x, direction.y);
    }
}; // struct Ray2

// 3D

struct Ray
{

    Vector3 origin;
    Vector3 direction;

    inline Ray() : origin(0.0f, 0.0f, 0.0f),direction(0.0f, 0.0f, 1.0f){ }
    explicit inline Ray(const Vector3 & origin, const Vector3& direction) : origin(origin), direction(direction) { }

    inline Vector3 pointAt(float t) const { return origin + direction.normalized() * t; }

    inline bool isZero() const { return (origin == Vector3::Zero && direction == Vector3::Zero);}


    // inline
    void print() const
    {
        printf("Ray(%.2f, %.2f, %.2f; %.2f, %.2f, %.2f)",
               origin.x, origin.y, origin.z,
               direction.x, direction.y, direction.z);
    }

}; // struct Ray

#endif // RAY_H
