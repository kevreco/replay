#include "vector4.h"

#include "vector3.h"

Vector4::Vector4() :
    x(0), y(0), z(0), w(0)
{

}

Vector4::Vector4(float x, float y, float z, float w) :
    x(x), y(y), z(z), w(w)
{

}

Vector4::Vector4(const float * values) :
    x(*(values)),
    y(*(values + 1)),
    z(*(values + 2)),
    w(*(values + 3))
{
}

Vector4::Vector4(const Vector4& v) :
    x(v.x),
    y(v.y),
    z(v.z),
    w(v.w)
{}

Vector4::Vector4(const Vector3& v) :
    x(v.x),
    y(v.y),
    z(v.z),
    w(1.0f)
{}

Vector4::Vector4(const Vector3& v, float value):
    x(v.x),
    y(v.y),
    z(v.z),
    w(value)
{}


Vector4 Vector4::operator +(const Vector4& v) const
{
    return Vector4(x + v.x, y + v.y, z + v.z, w + v.w);
}

Vector4 Vector4::operator - (const Vector4& v) const
{
    return Vector4(x - v.x, y - v.y, z - v.z, w - v.w);
}

Vector4 Vector4::operator * (const float value) const
{
    return Vector4(x * value, y * value, z * value, w * value);
}

Vector4 Vector4::operator / (const float value) const
{
    if(value == 0.0f) return Vector4(0.0f, 0.0f, 0.0f, 0.0f);
    return Vector4(x / value, y / value, z / value, w / value);
}



bool Vector4::operator == (const Vector4& v) const
{
    if(x == v.x && y == v.y && z == v.z && w == v.w){
        return true;
    }

    return false;
}

bool Vector4::operator != (const Vector4& v) const
{
    return !((*this) == v);
}

void Vector4::operator += (const Vector4& v)
{
    x += v.x;
    y += v.y;
    z += v.z;
    w += v.w;
}

void Vector4::operator -= (const Vector4& v)
{
    x -= v.x;
    y -= v.y;
    z -= v.z;
    w -= v.w;
}

void Vector4::operator *= (const float value)
{
    x *= value;
    y *= value;
    z *= value;
    w *= value;
}

void Vector4::operator /= (const float value)
{
    if(value == 0.0f){
        return;
    }
    else
    {
        x /= value;
        y /= value;
        z /= value;
        w /= value;
    }
}

Vector4 Vector4::operator - (void) const
{
    return Vector4(-x, -y, -z, -w);
}

Vector4 Vector4::operator + (void) const
{
    return (*this);
}

Vector4::operator float* () const
{
    return (float*) this;
}

Vector4::operator const float* () const
{
    return (const float*) this;
}

float Vector4::operator[](int i) const
{
    return data[i];
}

float& Vector4::operator[](int i)
{
    return data[i];
}

void Vector4::set(float x, float y, float z, float w)
{
    this->x = x;
    this->y = y;
    this->z = z;
    this->w = w;
}

float Vector4::dotProduct(const Vector4& v)
{
    return x*v.x + y*v.y + z*v.z + w*v.w;
}

Vector4 Vector4::lerp(const Vector4& v, float factor) const
{
    return ((*this) * (1.0f - factor)) + (v * factor);
}

Vector4 operator * (float factor, const Vector4& v)
{
    return v * factor;
}
