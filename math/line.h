#ifndef LINE_H
#define LINE_H

#include "vector2.h"
#include "vector3.h"

// Line is considered as infinite line

// 2D
struct Line2 {

    Vector2 point;
    Vector2 direction;

    static Line2 FromDirection(const Vector2& p, const Vector2& dir) {
        Line2 l;
        l.point = p;
        l.direction = dir;
        return l;
    }

    static Line2 FromPoint(const Vector2& p1, const Vector2& p2) {
        return FromDirection(p1, p2 - p1);
    }
};

// 3D
struct Line3
{

    Vector3 point;
    Vector3 direction;

    Line3(const Vector3& p, const Vector3& dir) :
    point(p),
    direction(dir) {

    }

};

#endif // LINE_H
