#include "vector2.h"

const Vector2 Vector2::Zero(0, 0);
const Vector2 Vector2::Xpos(1.0f, 0.0f);
const Vector2 Vector2::Xneg(-1.0f, 0.0f);
const Vector2 Vector2::Ypos(0.0f, 1.0f);
const Vector2 Vector2::Yneg(0.0f, -1.0f);

void Vector2::normalize()
{
    float length;
    float factor;

    length = this->length();

    if(length == 1 || length == 0){
        return;
    }

    factor = 1.0f / length;
    x *= factor;
    y *= factor;
}

void Vector2::rotate(const float angleInDeg)
{
    if(angleInDeg == 0.0){
        return;
    }

    float a = (float)(angleInDeg * M::DEG_TO_RAD);
    float sinAngle = sinf(a);
    float cosAngle = cosf(a);

    float rx = x * cosAngle - y * sinAngle;
    float ry = x * sinAngle + y * cosAngle;

    x = rx;
    y = ry;
}

Vector2 Vector2::rotated(double angleInRad) const
{
    if(angleInRad == 0.0){
        return (*this);
    }

    double sinAngle = sin(angleInRad);
    double cosAngle = cos(angleInRad);

    return Vector2(
        x * cosAngle - y * sinAngle,
        x * sinAngle + y * cosAngle
    );
}

Vector2 Vector2::rotated(const Vector2& point, double angleInRad) const
{
    if(angleInRad == 0.0)
        return (*this);


    double sinAngle = M::Sin(angleInRad);
    double cosAngle = M::Cos(angleInRad);
    Vector2 temp =   (*this) - point;

    return Vector2(
        point.x + (temp.x * cosAngle - temp.y * sinAngle),
        point.y + (temp.x * sinAngle + temp.y * cosAngle)
    );

}

Vector2 Vector2::normalized() const
{
    Vector2 result(*this);
    result.normalize();

    return result;
}
