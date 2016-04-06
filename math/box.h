#ifndef BOX_H
#define BOX_H

#include "vector2.h"
#include "vector3.h"
#include "matrix4.h" // for Box3

// Axis-aligned bounding box 2D

struct Box2
{

    Vector2 min, max;

    Box2() :
        min(),
        max() {
    }

    Box2(float left, float bottom, float right, float top) :
        min(left, bottom), max(right, top)
    {
    }

    explicit Box2(const Vector2& min, const Vector2& max) :
        min(min),
        max(max)
    {
    }

    Box2(const Vector2& position, float width, float height) : min(position)
    {
        set_width(width);
        set_height(height);
    }

    inline void set_width(float w) { max.x = min.x+w;}
    inline float width() const { return max.x - min.x;}

    inline void set_height(float h) { max.y = min.y+h;}
    inline float height() const { return max.y-min.y;}

    inline Vector2 size() const { return Vector2(width(),height());}

    // center of x axis
    inline float center_x() const { return min.x + (0.5f) * width();}

    // center of y  axis
    inline float center_y() const { return min.y + (0.5f) * height();}

    inline Vector2 center() const {return Vector2(center_x(), center_y());}

    inline void set_min(const Vector2& min) { this->min = min;}
    inline void set_max(const Vector2& max) { this->max = max;}

    inline bool intersects(const Box2& b) const
    {
        return (min <= b.max && max >= b.min);
    }

    inline bool intersects_exclusive(const Box2& b) const
    {

        if (max.x <= b.min.x) return false; // a is left of b
         if (min.x >= b.max.x) return false; // a is right of b
         if (max.y <= b.min.y) return false; // a is above b
         if (min.y >= b.max.y) return false; // a is below b
         return true; // boxes overlap
    }


    // Normalized Box to get positive values
    Box2 normalized() const
    {
        Box2 box;
        if (max.x < min.x) {
            box.min.x = max.x;
            box.max.x = min.x;
        } else {
            box.min.x = min.x;
            box.max.x = max.x;
        }

        if (max.y < min.y) {
            box.min.y = max.y;
            box.max.y = min.y;
        } else {
            box.min.y = min.y;
            box.max.y = max.y;
        }
        return box;
    }

    void translate(const Vector2& v) {
        min += v;
        max += v;
    }

    // If the point if contains inside the box or on the side
    bool contains_inclusive(const Vector2& p) const
    {
        return (p.x >= min.x && p.x <= max.x &&
                p.y >= min.y && p.y <= max.y);
    }

    // If the point is strictly inside the box
    bool contains(const Vector2& p) const
    {
        return (p.x > min.x && p.x < max.x &&
                p.y > min.y && p.y < max.y);
    }

    void normalize() { *this = normalized();}

    // inline
    void print() const
    {
        printf("box(%g, %g; %g, %g)", min.x, min.y, max.x, max.y);
    }

}; // struct Box2


// Axis-aligned bounding box 3D
// Also called Box3, square prism or rectangular parallelepiped class

struct Box3
{

    Vector3 min, max;

    Box3()
    {
    }

    Box3(float left, float bottom, float front, float right, float top, float back) :
        min(left, bottom, front), max(right, top, back)
    {
    }

    explicit Box3(const Vector3& min, const Vector3& max) : min(min), max(max)
    {
    }

    Box3(const Vector3& position, float width, float height, float length) : min(position)
    {
        set_width(width);
        set_height(height);
        set_length(length);
    }

    inline void set_width(float w) { max.x = min.x+w;}
    inline float width() const { return max.x - min.x;}

    inline void set_height(float h) { max.y = min.y+h;}
    inline float height() const { return max.y-min.y;}

    inline void set_length(float l) { max.z = min.z+l;}
    inline float length() const { return max.z-min.z;}

    inline Vector3 size() const { return Vector3(width(),height(),length());}

    // center of X axis
    inline float center_x() const { return min.x + (0.5f)*width();}

    // center of Y axis
    inline float center_y() const { return min.y + (0.5f) * height();}

    // center of Z axis
    inline float center_z() const { return min.z + (0.5f) * length();}

    inline Vector3 center() const {return Vector3(center_x(),center_y(),center_z());}

    inline void set_min(const Vector3& min) { this->min = min;}
    inline void set_max(const Vector3& max) { this->max = max;}

    inline bool intersects(const Box3& b) const
    {
        return (min <= b.max && max >= b.min);
    }


    // Normalized Box to get positive values
    Box3 normalized() const
    {
        Box3 box;
        if (max.x < min.x) {
            box.min.x = max.x;
            box.max.x = min.x;
        } else {
            box.min.x = min.x;
            box.max.x = max.x;
        }

        if (max.y < min.y) {
            box.min.y = max.y;
            box.max.y = min.y;
        } else {
            box.min.y = min.y;
            box.max.y = max.y;
        }

        if (max.z < min.z) {
            box.min.z = max.z ;
            box.max.z = min.z;
        } else {
            box.min.z = min.z;
            box.max.z = max.z ;
        }
        return box;
    }

    void expand_bounds_with(const Vector3& p)
    {
        if (p.x < min.x)
            min.x = p.x;
        else if (p.x > max.x)
            max.x = p.x;
        if (p.y < min.y)
            min.y = p.y;
        else if (p.y > max.y)
            max.y = p.y;
        if (p.z < min.z)
            min.z = p.z;
        else if (p.z > max.z)
            max.z = p.z;
    }

    void expand_bounds_with(const Box3& b)
    {
        expand_bounds_with(b.min);
        expand_bounds_with(b.max);
    }
    // Expanded by transform all points of Box3 and create a fit Box3 of the new "obb"
    Box3 expand_by_transform(const Matrix4& m) const
    {
        Box3 box;
        box.expand_bounds_with(m * min);
        box.expand_bounds_with(m * Vector3(min.x, min.y, max.z));
        box.expand_bounds_with(m * Vector3(min.x, max.y, max.z));
        box.expand_bounds_with(m * Vector3(min.x, max.y, min.z));
        box.expand_bounds_with(m * Vector3(max.x, min.y, min.z));
        box.expand_bounds_with(m * Vector3(max.x, max.y, min.z));
        box.expand_bounds_with(m * Vector3(max.x, min.y, max.z));
        box.expand_bounds_with(m * max);
        return box;

    }
    void merge(const Vector3& point)
    {

        max.max(point);
        min.min(point);

    }

    void translate(const Vector3& v) {
        min += v;
        max += v;
    }

    void transform_box_ex(const Matrix4& m) {

        if (m.isIdentity())
            return;


        const Vector3 Amin = min;
        const Vector3 Amax = max;

        min.x = m[12];
        min.y = m[13];
        min.z = m[14];
        max = min;


        for (size_t i = 0; i < 3; ++i) {
            for (size_t j = 0; j < 3; ++j) {
                const float value = m[j * 4 + i];
                const float a = value * Amin[j];
                const float b = value * Amax[j];

                if (a < b) {
                    min[i] += a;
                    max[i] += b;
                } else{
                    min[i] += b;
                    max[i] += a;
                }
            }
        }
    }


    // If the point if contains inside the box or on the side
    bool contains_inclusive(const Vector3& p) const
    {
        return (p.x >= min.x && p.x <= max.x &&
                p.y >= min.y && p.y <= max.y &&
                p.z >= min.z && p.z <= max.z);
    }

    // If the point is strictly inside the box
    bool contains(const Vector3& p) const
    {
        return (p.x > min.x && p.x < max.x &&
                p.y > min.y && p.y < max.y &&
                p.z > min.z && p.z < max.z);
    }


    void normalize() { *this = normalized();}

    void print() const
    {
        printf("(%.2f, %.2f, %.2f; %.2f, %.2f, %.2f)", min.x, min.y, min.z, max.x, max.y, min.z);
    }

}; // struct Box3


#endif // BOX_H
