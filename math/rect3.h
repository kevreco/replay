#ifndef RECT3_H
#define RECT3_H

#include "vector3.h"
#include "vector2.h"

struct Rect3
{
    Vector3 center; // center point of rectangle
    Vector3 axis[2]; // unit vectors determining local x and y axes for the rectangle
    Vector2 size; // the halfwidth extents of the rectangle along the axes

    inline Rect3() : center() {
        axis[0] = Vector3::Xpos; //local x axis
        axis[1] = Vector3::Ypos; //local y axis
        size[0] = size[1] = 0.5f;
    }

    inline Rect3(const Vector3 & center, float width, float height) : center(center)
    {
        size[0] = width;
        size[1] = height;
    }

    inline Rect3(const Vector3 & center, Vector3 axis[2]) : center(center)
    {
        this->axis[0] = axis[0];
        this->axis[1] = axis[1];
        size[0] = size[1] = 0.5;
    }

    inline Rect3(const Vector3 & center, const Vector3& axeX, const Vector3& axeY) : center(center)
    {
        this->axis[0] = axeX;
        this->axis[1] = axeY;
        size[0] = size[1] = 0.5;
    }

    inline Rect3(const Vector3 & center, const Vector3& axeX, const Vector3& axeY, float width, float height) : center(center)
    {
        this->axis[0] = axeX;
        this->axis[1] = axeY;
        this->size[0] = width;
        this->size[1] = height;

    }

    inline Rect3(const Vector3 & center, const Vector3& axeX, const Vector3& axeY, float size[2] ) : center(center)
    {
        this->axis[0] = axeX;
        this->axis[1] = axeY;
        this->size[0] = size[0];
        this->size[1] = size[1];
    }


}; //end class Rect3

#endif // RECT3_H
