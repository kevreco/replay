#ifndef VECTOR2_H
#define VECTOR2_H

#include "maths.h"

#include "stdio.h"
#include "assert.h"

struct Vector2
{
    typedef float float2[2];

    union {
        struct {
            float x, y;
        };
        float2 data;
    };

    static const Vector2 Xpos; // Positive X
    static const Vector2 Xneg; // Negative X
    static const Vector2 Ypos;
    static const Vector2 Yneg;
    
    Vector2():
        x(),
        y()
    {}
    
    Vector2(float x, float y)
    {
        this->x = x;
        this->y = y;
    }

    Vector2(const Vector2 & v):
        x(v.x),
        y(v.y)
    {}

    //public:


    /* Vector2 operator-() const
    {
        return Vector2(-x, -y);
    }*/

    Vector2& operator=(const Vector2& vec)
    {
        x = vec.x; y = vec.y; return *this;
    }

    Vector2  operator +(const Vector2& vec) const
    {
        return Vector2(x + vec.x, y + vec.y);
    }

    Vector2& operator +=(const Vector2& vec)
    {
        x+=vec.x; y+=vec.y; return *this;
    }

    Vector2  operator +(const float val) const
    {
        return Vector2(x+val, y+val);
    }

    Vector2& operator +=(const float val)
    {
        x+=val; y+=val;  return *this;
    }

    Vector2  operator -(const Vector2& vec) const
    {
        return Vector2(x - vec.x, y - vec.y);
    }

    Vector2& operator -=(const Vector2& vec)
    {
        x-=vec.x; y-=vec.y; return *this;
    }

    Vector2  operator -(const float val) const
    {
        return Vector2(x - val, y - val);
    }

    Vector2& operator -=(const float val)
    {
        x-=val; y-=val; return *this;
    }

    Vector2  operator *(const Vector2& vec) const
    {
        return Vector2(x*vec.x, y*vec.y);
    }

    Vector2& operator *=(const Vector2& vec)
    {
        x*=vec.x; y*=vec.y; return *this;
    }

    Vector2  operator *(const float v) const
    {
        return Vector2(x * v, y * v);
    }

    Vector2& operator *=(const float v)
    {
        x*=v; y*=v; return *this;
    }

    Vector2  operator /(const Vector2& vec) const
    {
        return Vector2(x/vec.x, y/vec.y);
    }

    Vector2& operator /=(const Vector2& vec)
    {
        x/=vec.x; y/=vec.y;  return *this;
    }

    Vector2  operator /(const float val) const
    {
        assert(fabs(val) > 1.0E-10 && "vector2::operator / : dividing by a null value");
        float iVal = 1.0f / val;
        return Vector2(x*iVal, y*iVal);
    }
    Vector2& operator /=(const float val)
    {
        assert(fabs(val) > 1.0E-10 && "vector2::operator / : dividing by a null value");
        float iVal = 1.0f / val;
        x *= iVal; y *= iVal;
        return *this;
    }

    inline bool operator == (const Vector2 & v) const
    {
        return (x == v.x && y == v.y);
    }

    inline bool operator != (const Vector2 & v) const
    {
        return !((*this) == v);
    }

    bool operator <(const Vector2& vec) const
    {
        return (x < vec.x && y < vec.y);
    }
    bool operator >(const Vector2& vec) const
    {
        return (x > vec.x && y > vec.y);
    }
    bool operator <=(const Vector2& vec) const
    {
        return (x <= vec.x && y <= vec.y);
    }
    bool operator >=(const Vector2& vec) const
    {
        return (x >= vec.x && y >= vec.y);
    }


    inline Vector2 operator -() const
    {
        return Vector2(-x, -y);
    }

    inline Vector2 operator + () const
    {
        return *this;
    }

    //inline operator float* () const
    //	{
    //return (float*) this;
    //}

    //inline operator const float* () const
    //	{
    //    return (const float*) this;
    //	}

    float operator[](int i) const
    {
        return data[i];
    }

    float& operator[](int i)
    {
        return data[i];
    }

    inline void set(float x, float y)
    {
        this->x = x;
        this->y = y;
    }

    Vector2 ortho() const {
        return Vector2(y, -x);
    }

    static float Dot(const Vector2& v1, const Vector2& v2) {
        return v1.x * v2.x + v1.y * v2.y;
    }

    static float Cross(const Vector2& v1, const Vector2& v2) {
        return (v1.x * v2.y) - (v1.y * v2.x);
    }

    void normalize();
    void rotate(const float angle);

    Vector2 normalized() const;
    Vector2 rotated(double angleInRad) const;
    Vector2 rotated(const Vector2& point, double angleInRad) const;

    // inline
    float length() const
    {
        return sqrtf((x*x) + (y*y));
    }

    // inline
    float squaredLength() const
    {
        return ((x*x) + (y*y));
    }

    // inline
    float dotProduct(const Vector2 & v) const
    {
        return (x * v.x) + (y * v.y);
    }

    // inline
    Vector2 lerp(const Vector2& v, float factor) const
    {
        return ((*this) * (1.0f - factor)) + (v * factor);
    }

    static const Vector2 Zero;

    // inline
    void print() const
    {
        printf("(%.2f, %.2f)", x, y);
    }

}; // struct Vector2


#endif // VECTOR2_H
