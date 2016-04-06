#ifndef MATRIX4_H
#define MATRIX4_H

#include <cmath>

#include <math/vector3.h>
#include <math/vector4.h>
#include <math/quaternion.h>

typedef float value_type;
typedef float float16[16];
typedef float float4x4[4][4];

struct Matrix4 {



    union {

        struct { float16 data; };
        // @TODO find a better name
        struct { float4x4 data4x4; };
        struct {
            float m11, m12, m13, m14,
            m21, m22, m23, m24,
            m31, m32, m33, m34,
            m41, m42, m43, m44;
        };
    };
public:

    //Operator
    Matrix4 operator +(const Matrix4& mat) const;
    Matrix4 operator -(const Matrix4& mat) const;
    Matrix4 operator *(const Matrix4& mat) const;

    Matrix4 operator *(float value) const;
    Matrix4 operator /(float value) const;

    Matrix4& operator =(const Matrix4& mat);
    Matrix4& operator =(const Matrix3& mat);
    Matrix4& operator =(const Vector3 axis[3]);


    Matrix4& operator +=(const Matrix4& mat);
    Matrix4& operator -=(const Matrix4& mat);
    Matrix4& operator *=(const Matrix4& mat);
    Matrix4& operator *=(float value);
    Matrix4& operator /=(float value);

    bool operator ==(const Matrix4& mat) const;
    bool operator !=(const Matrix4& mat) const;

    Matrix4 operator -() const;

    Vector4 operator *(const Vector4& vec) const;
    Vector3 operator *(const Vector3& vec) const;
    Vector2 multiplied(const Vector2& vec) const;

    inline float operator [](int i) const { return data[i];}
    inline float& operator [](int i) { return data[i]; }

public:
    typedef float* iterator;
    iterator begin() {
        return data;
    }
    iterator end() {
        return data + 16;
    }

public:

    static const Matrix4 Zero;
    static const Matrix4 Identity;

    bool isIdentity() const { return (*this == Identity); }

    float determinant() const;

    void scale(const Vector3& scaleFactor);
    inline void scale(float x, float y, float z) { scale(Vector3(x, y, z));} //overload
    void rotate(const Quaternion& quat);
    void rotate(const Vector3& axis, float angle);
    void translate(const Vector3& vec);

    void invert();
    void transpose();
    void affine_invert();
    void invert_transpose();
    void affine_invert_transpose();

    Matrix4 rotated(const Quaternion& quat) const;
    Matrix4 rotated(const Vector3& axis, float angle) const;
    Matrix4 translated(const Vector3& vec) const;

    Matrix4 inverted() const;
    Matrix4 transposed() const;
    Matrix4 inverse_transposed() const;
    Matrix4 affine_inverted() const;
    Matrix4 affine_inverse_transposed() const;

    Vector3 rotated_vector(const Vector3& vec) const;
    Vector3 inverse_rotated_vector(const Vector3& vec) const;

    void set_scale(const Vector3& scaleFactor); // set identity and set scale
    void set_scale(float size); // set identity and set scale
    void set_translation(const Vector3& translation); // set identity and set translation
    void set_translation_part(const Vector3& translation); // set translation
    void set_rotation_axis(float angle, const Vector3& axis);
    void set_rotation_x(float angle);
    void set_rotation_y(float angle);
    void set_rotation_z(float angle);
    void set_rotation_euler(float angleX, float angleY, float angleZ); //set identity and set translation
    void set_rotation_part_euler(float angleX, float angleY, float angleZ); //set identity and set translation
    void set_rotation_part_euler(const Vector3& rotations)
    {
        set_rotation_part_euler(rotations.x, rotations.y, rotations.z);
    } //overload

    Vector3 get_translation_part() const;

    Vector3 x_axis() const
    {
        return Vector3(data4x4[0][0], data4x4[0][1], data4x4[0][2]);
    }

    Vector3 y_axis() const
    {
        return  Vector3(data4x4[1][0], data4x4[1][1], data4x4[1][2]);
    }
    Vector3 z_axis() const
    {
        return  Vector3(data4x4[2][0], data4x4[2][1], data4x4[2][2]);
    }

    Vector3 get_euler_angles() const;

    void print() const
    {
        printf("\n(");
        printf("%g, %g, %g, %g\n", data[0], data[1], data[2], data[3]);
        printf("%g, %g, %g, %g\n", data[4], data[5], data[6], data[7]);
        printf("%g, %g, %g, %g\n", data[8], data[9], data[10], data[11]);
        printf("%g, %g, %g, %g\n", data[12], data[13], data[14], data[15]);
        printf(")");
    }

    // Construit une matrice de projection
    void SetProjectionMatrix (const float fov_, const float aspect_, const float nearPlane_, const float farPlane_)
    {
        float maxY = nearPlane_ * std::tan(fov_ * 3.14159256f / 360.0f);
        float minY = -maxY;
        float minX = minY * aspect_;
        float maxX = maxY * aspect_;

        //memset(data, 0, sizeof(Matrix4));

        data[0] = 2.0f * nearPlane_ / (maxX - minX);
        data[5] = 2.0f * nearPlane_ / (maxY - minY);
        data[8] = (maxX + minX) / (maxX - minX);
        data[9] = (maxY + minY) / (maxY - minY);
        data[10] = -(farPlane_ + nearPlane_) / (farPlane_ - nearPlane_);
        data[11] = -1.0f;
        data[14] = -(2.0f * farPlane_ * nearPlane_) / (farPlane_ - nearPlane_);
    }

    // Construit une matrice de transformation toute simple (uniquement translation)
    void SetTransformationMatrix (const float x_, const float y_, const float z_)
    {
        data[12] += x_;
        data[13] += y_;
        data[14] += z_;
    }

public:

    /// ASSUMING WE ARE IN LEFT HANDED AXIS

    // Left Handed
    static Matrix4 ModelViewFromDirectionLH(const Vector3& pos, const Vector3& dir, const Vector3& up = Vector3(0, 1.0, 0));
    static Matrix4 ModelViewFromLookAtLH(const Vector3& pos, const Vector3& target, const Vector3& up = Vector3(0, 1.0, 0));
    static Matrix4 ModelViewFromAxisLH(const Vector3& pos, const Vector3& x_axis, const Vector3& y_axis, const Vector3& z_axis);

    static Matrix4 OrthoLH(const Vector4& corners, const Vector2& zNearFar);
    static Matrix4 OrthoLH(const Vector2& min, const Vector2& max, const Vector2& zNearFar, bool flip_y_axis = false);
    static Matrix4 OrthoLHB(const Vector2& min, const Vector2& max, const Vector2& zNearFar, bool flip_y_axis = false);

    static Matrix4 OrthoLHWithWidth(const Vector4& corners, const Vector2& zNearFar);

    // Only use to make ImGui coordinate system working with Matrix from OrthoLH
    //static Matrix4 ImGuiOrthoLH(const Vector2& zNearFar, const Vector2& vpScaling, float screen_w, float screen_h);

    // @TODO Right Handed version


    // @TODO
    //static void FromRows(const Vector3& row1, const Vector3& row2, const Vector3& row3);
    // static void FromColumns(const Vector3& col1, const Vector3& col2, const Vector3& col3);

    //@TODO : check which of both function below is more efficient
    //matrix in parmater don't need to be initialized
    static void FromQuaternion(const Quaternion& quat, Matrix4& mat);
    //matrix in parmater don't need to be initialized
    static void FromQuaternion2(const Quaternion& quat, Matrix4& mat);

}; //end class Matrix4

#endif //MATRIX4_H
