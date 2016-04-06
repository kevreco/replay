#ifndef VECTOR4_H
#define VECTOR4_H

#include "vector3.h" // only used for function xyz()

class Vector4
{

public:

    union {
        float data[4];
        struct {
            float x, y, z, w;

        };
        struct {
            float r, g, b, a;
        };
        struct {
            float top, left, bottom, right;
        };
    };

public:

    Vector4();
    Vector4(float x, float y, float z, float w);
    Vector4(const float * values);
    Vector4(const Vector4& v);
    Vector4(const Vector3& v);
    Vector4(const Vector3& v, float value);

public:

public:

    Vector4 operator + (const Vector4& v) const;

    Vector4 operator - (const Vector4& v) const;

    Vector4 operator * (const float value) const;

    Vector4 operator / (const float value) const;

    bool operator == (const Vector4& v) const;

    bool operator != (const Vector4& v) const;
    void operator += (const Vector4& v);

    void operator -= (const Vector4& v);

    void operator *= (const float value);

    void operator /= (const float value);

    Vector4 operator - (void) const;

    Vector4 operator + (void) const;

    operator float* () const;

    operator const float* () const;

    float operator[](int i) const;

    float& operator[](int i);

    Vector3 xyz() {
        return Vector3(x, y, z);
    }

    void set(float x, float y, float z, float w);

    float dotProduct(const Vector4& v);

    Vector4 lerp(const Vector4& v, float factor) const;

    // inline
    void print() const
    {
        printf("(%.2f, %.2f, %.2f, %.2f)", x, y, z, w);
    }

}; // Vector4

#endif // VECTOR4_H
