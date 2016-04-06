#ifndef QUATERNION_H
#define QUATERNION_H

#include "maths.h"
#include "vector3.h"

class Matrix4;
class Matrix3;

struct Quaternion
{


    Quaternion();

    Quaternion(float x, float y, float z, float w);
    Quaternion(double x, double y, double z, double w);

    Quaternion(const Vector3 & axis, const float & angleInRadian);
    Quaternion(const float & scalar, const Vector3 & v);
    Quaternion(const Quaternion & quaternion);

    //Quaternion(const float& pitchX,const float& rollY,const float& yawZ, bool useDegrees=true);

    void set_value(float x, float y=0.0f, float z=0.0f, float w=1.0f);
    void set_value_from_axis(const Vector3& axis, float angleInRadian);
    void set_value_from_euler(const float pitchX, const float rollY, const float yawZ, bool useDegrees=true);
    void set_value_from_euler_rad(const float pitchX, const float rollY, const float yawZ);


    inline void invert() { x=-x, y=-y, z=-z; }
    inline Quaternion conjugate() const { return Quaternion(-x, -y, -z, w); }
    inline Quaternion inversed() const { return Quaternion(-x, -y, -z, w); }
    //inline Quaternion inverse() const { return inversed(); }

    //dot product
    inline float dot(const Quaternion& q) const
    {
         return (x * q.x) + (y * q.y) + (z * q.z) + (w * q.w);
    }

    inline Vector3 vector() const { return Vector3(x, y, z);}

    Vector3 axis() const;
    Vector3 getAxis() const;

    // Get angle in radian
    float angle() const;
    // Get angle in radian
    inline float angle_in_rad() const { return angle(); } // overload
    inline float angle_in_deg() const { return angle() * M::RAD_TO_DEG; }

    void get_angle_axis(float& angle, Vector3& axis) const;

    void normalize();
    float  magnitude() const;

    Vector3 to_euler() const;

    inline Vector3 to_euler_in_rad() const {return to_euler();} // overload
    Vector3 to_euler_in_deg() const;

    Quaternion & operator =(const Quaternion& quaternionToCopy);
    Quaternion operator +( const Quaternion & quaternion) const;

    Quaternion operator *(const Quaternion & quaternion) const;
    Quaternion operator *(float s) const;
    Quaternion& operator *=(const Quaternion &q);
    Quaternion& operator *=(float s);

    inline const float& operator[](int i) const { return q[i];}

    inline float& operator[](int i)  { return q[i];}

    inline bool operator ==(const Quaternion & quat) const
	{
		return (x == quat.x && y == quat.y && z == quat.z && w == quat.w);
	}

	inline bool operator !=(const Quaternion & quat) const
	{
		return !((*this) == quat);
	}

    const float * matrix() const;

    void get_matrix(float m[4][4]) const;
    void get_matrix(float m[16]) const;

    void get_rotation_matrix(float m[3][3]) const;

    //vector3 operator* (const Vector3 & vectorR) const;
    Vector3 operator* (const Vector3& v) const;
    Vector3 rotate(const Vector3 & v) const;
    Vector3 inverse_rotate(const Vector3& v) const;
    //inline Vector3 inverseRotate(const Vector3& v) const{return inversedRotate(v);};

    void set_from_rotated_basis(const Vector3& x, const Vector3& y, const Vector3& z);
    void set_from_rotation_matrix(const double m[3][3]);

    Vector3 x_axis() const;
    Vector3 y_axis() const;
    Vector3 z_axis() const;
    Vector3 xy_axis() const;
    //Vector3 xz_axis() const;
    //Vector3 yz_axis() const;
    void get_xyz(Vector3 axis[3]) const;
    void get_xyz2(Vector3 axis[3]) const;
    void get_xy(Vector3 axis[2]) const;
    void get_xz(Vector3 axis[2]) const;
    void get_yz(Vector3 axis[2]) const;

    Quaternion mirrored(const Vector3& planeNorm);

    // inline
    void print() const {
        printf("(%f, %f, %f, %f)", x, y, z, w);
    }

public:

    static Quaternion FromEulerRadian(const float pitchX, const float rollY, const float yawZ);
    static Quaternion FromEulerDegree(const float pitchX, const float rollY, const float yawZ);

    static Quaternion FromAxisAngle(const Vector3& axis, float radianAngle);
    // Assuming Left Handed coordinate
    static Quaternion FromDirection(const Vector3& dir, const Vector3& up = Vector3::Ypos);
    static Quaternion FromMatrix(const Matrix4& m);
    static Quaternion BuildFromLookAt(const Vector3& pos, const Vector3& target);


public:
    static const Quaternion Identity;

    union {
        float q[4];
        struct {float x, y, z, w;} ;
    };

}; // end class Quaternion


inline Quaternion Interpolation(Quaternion q1, Quaternion q2, float ratio)
{
	float angle = q1.dot(q2);

    if (angle < 0.0f) {
		q1 *= -1.0f;
		angle *= -1.0f;
	}

	float scale;
	float invscale;

    if ((angle + 1.0f) > 0.05f) {
        if ((1.0f - angle) >= 0.05f) { // spherical interpolation
			const float theta = acosf(angle);
			const float invsintheta = 1.0f/sinf(theta);
			scale = sinf(theta * (1.0f-ratio)) * invsintheta;
			invscale = sinf(theta * ratio) * invsintheta;
		}
        else { // linear interploation
			scale = 1.0f - ratio;
			invscale = ratio;
		}
	}
    else {
        q2.set_value(-q1.y, q1.x, -q1.w, q1.z);
        scale = sinf(M::PI * (0.5f - ratio));
        invscale = sinf(M::PI * ratio);
	}

	return ((q1*scale) + (q2*invscale));
}

#endif // QUATERNION_H
