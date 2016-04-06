#include "vector3.h"

const Vector3 Vector3::Zero(0,0,0);
const Vector3 Vector3::Identity(1,1,1);
const Vector3 Vector3::Unit[6] ={{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};
const Vector3 Vector3::Xpos(1.0f, 0.0f, 0.0f);
const Vector3 Vector3::Xneg(-1.0f, 0.0f, 0.0f); //Negative X
const Vector3 Vector3::Ypos(0.0f, 1.0f, 0);
const Vector3 Vector3::Yneg(0.0f, -1.0f, 0.0f);
const Vector3 Vector3::Zpos(0.0f, 0.0f, 1.0f);
const Vector3 Vector3::Zneg(0.0f, 0.0f, -1.0);

Vector3::Vector3() :
    x(0), y(0), z(0)
{

}

Vector3::Vector3(const Vector3& v) :
    x(v.x), y(v.y), z(v.z)
{
}

Vector3::Vector3(float x, float y, float z) :
    x(x), y(y), z(z)
{

}

Vector3::Vector3(float val) :
    x(val), y(val), z(val)
{

}

Vector3 Vector3::operator +(const Vector3& vec) const
{
  return Vector3(x + vec.x, y + vec.y, z + vec.z);
}


Vector3 Vector3::operator -(const Vector3& vec) const
{
  return Vector3(x - vec.x, y - vec.y, z - vec.z);
}


Vector3 Vector3::operator *(const Vector3& vec) const
{
  return Vector3(x * vec.x, y * vec.y, z * vec.z);
}


Vector3 Vector3::operator /(const Vector3& vec) const
{
  return Vector3(x / vec.x, y / vec.y, z / vec.z);
}


Vector3 Vector3::operator +(const float val) const
{
  return Vector3(x + val, y + val, z + val);
}


Vector3 Vector3::operator -(const float val) const
{
  return Vector3(x - val, y - val, z - val);
}


Vector3 Vector3::operator *(const float val) const
{
    return Vector3(x * val, y * val, z * val);
}


Vector3 Vector3::operator /(const float val) const
{
    assert(M::Abs(val) > 1.0E-10 && "vector3::operator / : dividing by a null value");
    float iVal=(float)1.0/val; //inverse value
    return Vector3(x * iVal, y * iVal, z * iVal);
}


Vector3& Vector3::operator +=(const Vector3& vec)
{
    x += vec.x; y += vec.y; z += vec.z; return *this;
}


Vector3& Vector3::operator -=(const Vector3& vec)
{
  x -= vec.x; y -= vec.y; z -= vec.z; return *this;
}


Vector3& Vector3::operator *=(const Vector3& vec)
{
  x *= vec.x; y *= vec.y; z *= vec.z; return *this;
}


Vector3& Vector3::operator /=(const Vector3& vec)
{
  x /= vec.x; y /= vec.y; z /= vec.z; return *this;
}


Vector3& Vector3::operator +=(const float val)
{
  x += val; y += val; z += val; return *this;
}


Vector3& Vector3::operator -=(const float val)
{
  x -= val; y -= val; z -= val; return *this;
}


Vector3& Vector3::operator *=(const float val)
{
  x *= val; y *= val; z *= val; return *this;
}




Vector3& Vector3::operator /=(const float val)
{
  assert(M::Abs(val) > 1.0E-10 && "vector3::operator / : dividing by a null value");
  float iVal=(float)1.0/val; //inverse value
  x *= iVal; y *= iVal; z *= iVal; return *this;
}


bool Vector3::operator==(const Vector3& vec) const
{
  return (M::Equals(x, vec.x) && M::Equals(y , vec.y)&& M::Equals(z , vec.z));
}


bool Vector3::operator!=(const Vector3& vec) const
{
  return !(*this==vec);
}


Vector3 Vector3::operator -() const
{
  return Vector3(-x, -y, -z);
}


bool Vector3::operator <(const Vector3&vec) const
{
  if (x > vec.x)
    return false;
  else if (x < vec.x)
    return true;
  else if (x == vec.x)
  {
    if((y > vec.y))
      return false;
    else if (y < vec.y)
      return true;
    else if (y == vec.y)
    {
      if((z > vec.z))
        return false;
      else if (z < vec.z)
        return true;
      else
        return false;
    }
  }
  return false;
}


bool Vector3::operator >(const Vector3& vec) const
{
  return (x > vec.x && y > vec.y && z > vec.z);
}


bool Vector3::operator <=(const Vector3& vec) const
{
  return (x <= vec.x && y <= vec.y && z <= vec.z);
}


bool Vector3::operator >=(const Vector3& vec) const
{
  return (x >= vec.x && y >= vec.y && z >= vec.z);
}



void Vector3::set(float x, float y, float z)
{
  this->x = x;
  this->y = y;
  this->z = z;
}


void Vector3::setLength(double length)
{
    double magn = magnitude();
    if (magn > 1.0E-10) {
        double len = (1/magn)*length;
        this->x *= len;
        this->y *= len;
        this->z *= len;
    } else {
        this->x = 0.0f;
        this->y = 0.0f;
        this->z = 0.0f;
    }
}


Vector3 Vector3::withLength(double length) const
{
Vector3 vec;
 double magn = magnitude();
  if (magn > 1.0E-10)
  {
   double len = (1.0/magn)*length;
    vec.x = this->x * len;
    vec.y = this->y * len;
    vec.z = this->z * len;
  }
  else
  {
    vec.x = 0.0f;
    vec.y = 0.0f;
    vec.z = 0.0f;
  }
  return vec;
}


Vector3 Vector3::scaled(double length) const
{
  return Vector3(this->x * length, this->y * length, this->z * length);
}


void Vector3::normalize()
{
 double magn = magnitude();
  if (magn > 1.0E-10)
  {
   double invMagn= 1.0/magn;
    this->x *= invMagn;
    this->y *= invMagn;
    this->z *= invMagn;
  }
  else
  {
    this->x = 0.0f;
    this->y = 0.0f;
    this->z = 0.0f;
  }
}


Vector3 Vector3::normalized() const
{
Vector3 v;
 double magn = magnitude();

  if (magn > 1.0E-10)
  {
   double invMagn= 1.0f / magn;
    v.x = this->x * invMagn;
    v.y = this->y * invMagn;
    v.z = this->z * invMagn;
  }
  else
  {
    v.x = 0.0f;
    v.y = 0.0f;
    v.z = 0.0f;
  }

  return v;
}


double Vector3::magnitude() const
{
  return sqrt(x*x + y*y + z*z);
}


Vector3& Vector3::inverse()
{
  x *= -1.0f; y *= -1.0f; z *= -1.0f; return *this;
}


Vector3 Vector3::inversed() const
{
  return Vector3(x* -1.0f, y * -1.0f, z * -1.0f);
}


double Vector3::dot(const Vector3& vec) const //Dot product
{
  return (x*vec.x + y*vec.y + z*vec.z);
}


Vector3 Vector3::cross( const Vector3 &vec) const //Cross product
{
 return Vector3(y*vec.z - z*vec.y,
           z*vec.x - x*vec.z,
           x*vec.y - y*vec.x);
}


Vector3 Vector3::crossNormalize(const Vector3& vec) const //Cross and normalize
{
Vector3 cross(y*vec.z - z*vec.y,
          z*vec.x - x*vec.z,
          x*vec.y - y*vec.x);
  cross.normalize();
  return cross;
}

bool Vector3::isValid() const
{
  return (M::IsFinite(x) && M::IsFinite(y) && M::IsFinite(z));
}

void Vector3::projectOn(const Vector3& axisDirection)
{
  assert(axisDirection.squaredNorm() > 1.0E-10 && "vector3::projectOn : axis direction is not normalized");

  *this = axisDirection * ( (*this)*axisDirection / axisDirection.squaredNorm() ) ;
}


Vector3 Vector3::projectedOn(const Vector3& axisDirection) const
{
  assert(axisDirection.squaredNorm() > 1.0E-10 && "vector3::projectedOn : axis direction is not normalized");

  return axisDirection * ( (*this).dot(axisDirection) / axisDirection.squaredNorm() ) ;
}

//Static

Vector3 Vector3::Interpolate(const Vector3 & a,const Vector3 & b, const float ratio)
{
  return Vector3(a.x + ( ratio * ( b.x - a.x )),
            a.y + ( ratio * ( b.y - a.y )),
            a.z + ( ratio * ( b.z - a.z )) );
}
