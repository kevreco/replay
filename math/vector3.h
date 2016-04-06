#ifndef VECTOR3_H
#define VECTOR3_H

#include "vector2.h"
#include "maths.h"

#include <assert.h>
#include <limits>



struct Vector3
{
    typedef float float3[3];

    union {
        struct { float x; float y; float z; };
        float3 data;
    };

public:
    Vector3();
    Vector3(const Vector3& v);

    Vector3(float x, float y, float z);
    explicit Vector3(float val);

public:

    Vector2 xy() const { return Vector2(x, y); }
    //Operator
    Vector3 operator +(const Vector3& vec) const;
    Vector3 operator -(const Vector3& vec) const;
    Vector3 operator *(const Vector3& vec) const;
    Vector3 operator /(const Vector3& vec) const;
    Vector3 operator +(const float val) const;
    Vector3 operator -(const float val) const;
    Vector3 operator *(const float val) const;
    Vector3 operator /(const float val) const;

    //vector3  operator *(const float val) const;

    Vector3& operator +=(const Vector3& vec);
    Vector3& operator -=(const Vector3& vec);
    Vector3& operator *=(const Vector3& vec);
    Vector3& operator /=(const Vector3& vec);
    Vector3& operator +=(const float val);
    Vector3& operator -=(const float val);
    Vector3& operator *=(const float val);
    Vector3& operator /=(const float val);

    bool operator==(const Vector3& a) const;
    bool operator!=(const Vector3& v) const;

    Vector3 operator -() const;

    //Yep, strange operator for vector, it does not really make sense, but we need some for reason (sorting etc.)
    //TODO : choose an implementation among both -> (1,0,0) < (0,0,1)  == false ? (0,0,1) < (1,0,0)  == false ?
    bool operator <(const Vector3& vec) const;
    bool operator >(const Vector3& vec) const;
    bool operator <=(const Vector3& vec) const;
    bool operator >=(const Vector3& vec) const;

    float operator [](int i) const { return data[i];}
    float& operator [](int i) { return data[i];}

    inline operator const float* () const { return (const float*) this;}

    void set(float x, float y, float z);

    void normalize(); //noramlize "this" vector
    Vector3 normalized() const; //return a normalized copy

    void setLength(double length);
    Vector3 withLength(double length) const; //return a normalized copy with good length (normalize)
    Vector3 scaled(double length) const; //return a simple copy with good length

    double magnitude() const;
    inline double length() const { return magnitude(); } //overload

    Vector3& inverse(); //inverse "this" vector
    Vector3 inversed() const; //return a inversed copy

    double squaredLength() const { return squaredNorm(); }
    double squaredNorm() const { return x*x + y*y + z*z; }


    double dot(const Vector3& vec) const; //Dot product
    Vector3 cross( const Vector3&vec) const; //Cross product
    Vector3 crossNormalize(const Vector3& vec) const; //Cross and normalize

    inline bool isValid() const; //Returns true if all values are not "infinite" and are not "NaN" (Not a Number)

    inline void projectOn(const Vector3& axisDirection);
    Vector3 projectedOn(const Vector3& axisDirection) const;

    Vector3 getHorizontalAngle() const
    {
        Vector3 angle;

        angle.y = (float)(atan2(x, z) * M::RAD_TO_DEG);

        if (angle.y < 0.0f)
            angle.y += 360.0f;
        if (angle.y >= 360.0f)
            angle.y -= 360.0f;

        const double z1 = sqrt(x*x + z*z);

        angle.x = (float)(atan2((double)z1, (double)y) * M::RAD_TO_DEG - 90.0f);

        if (angle.x < 0.0f)
            angle.x += 360.0f;
        if (angle.x >= 360.0f)
            angle.x -= 360.0f;

        return angle;
    }


    bool equals(const Vector3& b)
    {
        return ( M::Equals(this->x, b.x) && M::Equals(this->y, b.y) && M::Equals(this->z, b.z) );
    }


    //Min max

    Vector3 min(float ceil) const
    {
        return Vector3(M::Min(x, ceil), M::Min(y, ceil), M::Min(z, ceil));
    }

    Vector3 min(const Vector3& ceil) const
    {
        return Vector3(M::Min(x, ceil.x), M::Min(y, ceil.y), M::Min(z, ceil.z));
    }

    Vector3 max(float floor) const
    {
        return Vector3(M::Max(x, floor), M::Max(y, floor), M::Max(z, floor));
    }

    Vector3 max(const Vector3& floor) const
    {
        return Vector3(M::Max(x, floor.x), M::Max(y, floor.y), M::Max(z, floor.z));
    }

    Vector3 clamp(float floor, float ceil) const
    {
        return min(ceil).max(floor);
    }

    Vector3 clamp(const Vector3& floor, const Vector3& ceil) const
    {
        return min(ceil).max(floor);
    }

    //Lerp
    static Vector3 Interpolate(const Vector3& a,const Vector3& b, const float ratio);


    /* union
    {
        struct { float x; float y;float z; };
        float v[3];
    };*/

    static const Vector3 Zero;
    static const Vector3 Identity;

    //TODO : choose the array, or the 6 axis separate ?

    static const Vector3 Unit[6];

    static const Vector3 Xpos;  //Positive X
    static const Vector3 Xneg; //Negative X
    static const Vector3 Ypos;
    static const Vector3 Yneg;
    static const Vector3 Zpos;
    static const Vector3 Zneg;

    void print()
    {
        printf("(%.2f, %.2f, %.2f)", x, y, z);
    }


}; // struct Vector3

enum {
    X = 0, NX, Y, NY, Z, NZ
};

#endif // VECTOR3_H

