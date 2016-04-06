#ifndef SEGMENT_H
#define SEGMENT_H

#include "vector2.h"

// 2D
struct Segment2
{
    Vector2 a;
    Vector2 b;

    Segment2(float a, float b, float c, float d) :
        a(a, b),
        b(c, d)
    {
    }

    Segment2(const Vector2& p1, const Vector2& p2) :
        a(p1),
        b(p2)
    {
    }

    void print() const
    {
        printf("(%.2f, %.2f; %.2f, %.2f)", a.x, a.y, b.x, b.y);
    }

}; // Segment2

#endif // SEGMENT_H
