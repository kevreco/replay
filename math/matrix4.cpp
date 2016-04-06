#include "matrix4.h"
#include "vector3.h"
//#include <Utils/Logger.h>

#include <memory.h>
#include <math/maths.h>

#include "matrix3.h"


///All fonction for matrix template


const Matrix4 Matrix4::Zero = {
    0.0f, 0.0f, 0.0f, 0.0f,
    0.0f, 0.0f, 0.0f, 0.0f,
    0.0f, 0.0f, 0.0f, 0.0f,
    0.0f, 0.0f, 0.0f, 0.0f
};


const Matrix4 Matrix4::Identity = {
    1.0f, 0.0f, 0.0f, 0.0f,
    0.0f, 1.0f, 0.0f, 0.0f,
    0.0f, 0.0f, 1.0f, 0.0f,
    0.0f, 0.0f, 0.0f, 1.0f
};

Matrix4 Matrix4::operator +(const Matrix4& mat) const
{
    return Matrix4{
        data[0] + mat.data[0],
                data[1] + mat.data[1],
                data[2] + mat.data[2],
                data[3] + mat.data[3],
                data[4] + mat.data[4],
                data[5] + mat.data[5],
                data[6] + mat.data[6],
                data[7] + mat.data[7],
                data[8] + mat.data[8],
                data[9] + mat.data[9],
                data[10] + mat.data[10],
                data[11] + mat.data[11],
                data[12] + mat.data[12],
                data[13] + mat.data[13],
                data[14] + mat.data[14],
                data[15] + mat.data[15]
    };
}


Matrix4 Matrix4::operator -(const Matrix4& mat) const
{
    return Matrix4{
        data[0] - mat.data[0],
                data[1] - mat.data[1],
                data[2] - mat.data[2],
                data[3] - mat.data[3],
                data[4] - mat.data[4],
                data[5] - mat.data[5],
                data[6] - mat.data[6],
                data[7] - mat.data[7],
                data[8] - mat.data[8],
                data[9] - mat.data[9],
                data[10] - mat.data[10],
                data[11] - mat.data[11],
                data[12] - mat.data[12],
                data[13] - mat.data[13],
                data[14] - mat.data[14],
                data[15] - mat.data[15]
    };
}


Matrix4 Matrix4::operator *(const Matrix4& m) const
{
    const Matrix4& a =(*this);
    Matrix4 r;

    r.data[0]  = a[0]*m.data[0] + a[4]*m.data[1] + a[8]*m.data[2] + a[12]*m.data[3];
    r.data[1]  = a[1]*m.data[0] + a[5]*m.data[1] + a[9]*m.data[2] + a[13]*m.data[3];
    r.data[2]  = a[2]*m.data[0] + a[6]*m.data[1] + a[10]*m.data[2] + a[14]*m.data[3];
    r.data[3]  = a[3]*m.data[0] + a[7]*m.data[1] + a[11]*m.data[2] + a[15]*m.data[3];

    r.data[4]  = a[0]*m.data[4] + a[4]*m.data[5] + a[8]*m.data[6] + a[12]*m.data[7];
    r.data[5]  = a[1]*m.data[4] + a[5]*m.data[5] + a[9]*m.data[6] + a[13]*m.data[7];
    r.data[6]  = a[2]*m.data[4] + a[6]*m.data[5] + a[10]*m.data[6] + a[14]*m.data[7];
    r.data[7]  = a[3]*m.data[4] + a[7]*m.data[5] + a[11]*m.data[6] + a[15]*m.data[7];

    r.data[8]  = a[0]*m.data[8] + a[4]*m.data[9] + a[8]*m.data[10] + a[12]*m.data[11];
    r.data[9]  = a[1]*m.data[8] + a[5]*m.data[9] + a[9]*m.data[10] + a[13]*m.data[11];
    r.data[10] = a[2]*m.data[8] + a[6]*m.data[9] + a[10]*m.data[10] + a[14]*m.data[11];
    r.data[11] = a[3]*m.data[8] + a[7]*m.data[9] + a[11]*m.data[10] + a[15]*m.data[11];

    r.data[12] = a[0]*m.data[12] + a[4]*m.data[13] + a[8]*m.data[14] + a[12]*m.data[15];
    r.data[13] = a[1]*m.data[12] + a[5]*m.data[13] + a[9]*m.data[14] + a[13]*m.data[15];
    r.data[14] = a[2]*m.data[12] + a[6]*m.data[13] + a[10]*m.data[14] + a[14]*m.data[15];
    r.data[15] = a[3]*m.data[12] + a[7]*m.data[13] + a[11]*m.data[14] + a[15]*m.data[15];

    return r;
    //    }

    /*
        //If bottom row of both matrices is (0, 0, 0, 1).
        if(	data[3]==0.0f && data[7]==0.0f &&
                data[11]==0.0f && data[15]==1.0f	&&
                mat.data[3]==0.0f && mat.data[7]==0.0f &&
                mat.data[11]==0.0f && mat.data[15]==1.0f)
        {
        return Matrix4(
                        data[0]*mat.data[0]+data[4]*mat.data[1]+data[8]*mat.data[2],
                        data[1]*mat.data[0]+data[5]*mat.data[1]+data[9]*mat.data[2],
                        data[2]*mat.data[0]+data[6]*mat.data[1]+data[10]*mat.data[2],
                        0.0f,
                        data[0]*mat.data[4]+data[4]*mat.data[5]+data[8]*mat.data[6],
                        data[1]*mat.data[4]+data[5]*mat.data[5]+data[9]*mat.data[6],
                        data[2]*mat.data[4]+data[6]*mat.data[5]+data[10]*mat.data[6],
                        0.0f,
                        data[0]*mat.data[8]+data[4]*mat.data[9]+data[8]*mat.data[10],
                        data[1]*mat.data[8]+data[5]*mat.data[9]+data[9]*mat.data[10],
                        data[2]*mat.data[8]+data[6]*mat.data[9]+data[10]*mat.data[10],
                        0.0f,
                        data[0]*mat.data[12]+data[4]*mat.data[13]+data[8]*mat.data[14]+data[12],
                        data[1]*mat.data[12]+data[5]*mat.data[13]+data[9]*mat.data[14]+data[13],
                        data[2]*mat.data[12]+data[6]*mat.data[13]+data[10]*mat.data[14]+data[14],
                        1.0f
                );
        }

        //If bottom row of first matrix is (0, 0, 0, 1).
        if(	data[3]==0.0f && data[7]==0.0f &&
                data[11]==0.0f && data[15]==1.0f)
        {
        return Matrix4(
                        data[0]*mat.data[0]+data[4]*mat.data[1]+data[8]*mat.data[2]+data[12]*mat.data[3],
                        data[1]*mat.data[0]+data[5]*mat.data[1]+data[9]*mat.data[2]+data[13]*mat.data[3],
                        data[2]*mat.data[0]+data[6]*mat.data[1]+data[10]*mat.data[2]+data[14]*mat.data[3],
                        mat.data[3],
                        data[0]*mat.data[4]+data[4]*mat.data[5]+data[8]*mat.data[6]+data[12]*mat.data[7],
                        data[1]*mat.data[4]+data[5]*mat.data[5]+data[9]*mat.data[6]+data[13]*mat.data[7],
                        data[2]*mat.data[4]+data[6]*mat.data[5]+data[10]*mat.data[6]+data[14]*mat.data[7],
                        mat.data[7],
                        data[0]*mat.data[8]+data[4]*mat.data[9]+data[8]*mat.data[10]+data[12]*mat.data[11],
                        data[1]*mat.data[8]+data[5]*mat.data[9]+data[9]*mat.data[10]+data[13]*mat.data[11],
                        data[2]*mat.data[8]+data[6]*mat.data[9]+data[10]*mat.data[10]+data[14]*mat.data[11],
                        mat.data[11],
                        data[0]*mat.data[12]+data[4]*mat.data[13]+data[8]*mat.data[14]+data[12]*mat.data[15],
                        data[1]*mat.data[12]+data[5]*mat.data[13]+data[9]*mat.data[14]+data[13]*mat.data[15],
                        data[2]*mat.data[12]+data[6]*mat.data[13]+data[10]*mat.data[14]+data[14]*mat.data[15],
                        mat.data[15]
                );
        }

        //If bottom row of 2nd matrix is (0, 0, 0, 1).
        if(	mat.data[3]==0.0f && mat.data[7]==0.0f &&
                mat.data[11]==0.0f && mat.data[15]==1.0f)
        {
        return Matrix4(
                        data[0]*mat.data[0]+data[4]*mat.data[1]+data[8]*mat.data[2],
                        data[1]*mat.data[0]+data[5]*mat.data[1]+data[9]*mat.data[2],
                        data[2]*mat.data[0]+data[6]*mat.data[1]+data[10]*mat.data[2],
                        data[3]*mat.data[0]+data[7]*mat.data[1]+data[11]*mat.data[2],
                        data[0]*mat.data[4]+data[4]*mat.data[5]+data[8]*mat.data[6],
                        data[1]*mat.data[4]+data[5]*mat.data[5]+data[9]*mat.data[6],
                        data[2]*mat.data[4]+data[6]*mat.data[5]+data[10]*mat.data[6],
                        data[3]*mat.data[4]+data[7]*mat.data[5]+data[11]*mat.data[6],
                        data[0]*mat.data[8]+data[4]*mat.data[9]+data[8]*mat.data[10],
                        data[1]*mat.data[8]+data[5]*mat.data[9]+data[9]*mat.data[10],
                        data[2]*mat.data[8]+data[6]*mat.data[9]+data[10]*mat.data[10],
                        data[3]*mat.data[8]+data[7]*mat.data[9]+data[11]*mat.data[10],
                        data[0]*mat.data[12]+data[4]*mat.data[13]+data[8]*mat.data[14]+data[12],
                        data[1]*mat.data[12]+data[5]*mat.data[13]+data[9]*mat.data[14]+data[13],
                        data[2]*mat.data[12]+data[6]*mat.data[13]+data[10]*mat.data[14]+data[14],
                        data[3]*mat.data[12]+data[7]*mat.data[13]+data[11]*mat.data[14]+data[15]
                );
        }

    return Matrix4(
                data[0]*mat.data[0]+data[4]*mat.data[1]+data[8]*mat.data[2]+data[12]*mat.data[3],
                data[1]*mat.data[0]+data[5]*mat.data[1]+data[9]*mat.data[2]+data[13]*mat.data[3],
                data[2]*mat.data[0]+data[6]*mat.data[1]+data[10]*mat.data[2]+data[14]*mat.data[3],
                data[3]*mat.data[0]+data[7]*mat.data[1]+data[11]*mat.data[2]+data[15]*mat.data[3],
                data[0]*mat.data[4]+data[4]*mat.data[5]+data[8]*mat.data[6]+data[12]*mat.data[7],
                data[1]*mat.data[4]+data[5]*mat.data[5]+data[9]*mat.data[6]+data[13]*mat.data[7],
                data[2]*mat.data[4]+data[6]*mat.data[5]+data[10]*mat.data[6]+data[14]*mat.data[7],
                data[3]*mat.data[4]+data[7]*mat.data[5]+data[11]*mat.data[6]+data[15]*mat.data[7],
                data[0]*mat.data[8]+data[4]*mat.data[9]+data[8]*mat.data[10]+data[12]*mat.data[11],
                data[1]*mat.data[8]+data[5]*mat.data[9]+data[9]*mat.data[10]+data[13]*mat.data[11],
                data[2]*mat.data[8]+data[6]*mat.data[9]+data[10]*mat.data[10]+data[14]*mat.data[11],
                data[3]*mat.data[8]+data[7]*mat.data[9]+data[11]*mat.data[10]+data[15]*mat.data[11],
                data[0]*mat.data[12]+data[4]*mat.data[13]+data[8]*mat.data[14]+data[12]*mat.data[15],
                data[1]*mat.data[12]+data[5]*mat.data[13]+data[9]*mat.data[14]+data[13]*mat.data[15],
                data[2]*mat.data[12]+data[6]*mat.data[13]+data[10]*mat.data[14]+data[14]*mat.data[15],
                data[3]*mat.data[12]+data[7]*mat.data[13]+data[11]*mat.data[14]+data[15]*mat.data[15]
        );*/
}


Matrix4 Matrix4::operator *(float value) const {
    return Matrix4{
        data[0] * value,
                data[1] * value,
                data[2] * value,
                data[3] * value,
                data[4] * value,
                data[5] * value,
                data[6] * value,
                data[7] * value,
                data[8] * value,
                data[9] * value,
                data[10] * value,
                data[11] * value,
                data[12] * value,
                data[13] * value,
                data[14] * value,
                data[15] * value
    };
}


Matrix4 Matrix4::operator /(float value) const {
    if(value == 0.0f || value == 1.0f)
        return (*this);

    return (*this) * (1.0f / value);
}

Matrix4& Matrix4::operator =(const Matrix4& mat) {
    data[0] = mat[0];
    data[1] = mat[1];
    data[2] = mat[2];
    data[3] = mat[3];
    data[4] = mat[4];
    data[5] = mat[5];
    data[6] = mat[6];
    data[7] = mat[7];
    data[8] = mat[8];
    data[9] = mat[9];
    data[10] = mat[10];
    data[11] = mat[11];
    data[12] = mat[12];
    data[13] = mat[13];
    data[14] = mat[14];
    data[15] = mat[15];

    return *this;
}

Matrix4& Matrix4::operator =(const Matrix3& m) {
    // @TODO optimize this function

    *this = Matrix4::Identity;

    data[0] = m.data3x3[0][0];
    data[1] = m.data3x3[0][1];
    data[2] = m.data3x3[0][2];

    data[4] = m.data3x3[1][0];
    data[5] = m.data3x3[1][1];
    data[6] = m.data3x3[1][2];

    data[8] = m.data3x3[2][0];
    data[9] = m.data3x3[2][1];
    data[10] = m.data3x3[2][2];

}

Matrix4& Matrix4::operator =(const Vector3 axis[3]) {
    // @TODO optimize this function

    *this = Matrix4::Identity;

    data[0] = axis[0].x;
    data[1] = axis[0].y;
    data[2] = axis[0].z;

    data[4] = axis[1].x;
    data[5] = axis[1].y;
    data[6] = axis[1].z;

    data[8] = axis[2].x;
    data[9] = axis[2].y;
    data[10] = axis[2].z;
}



Matrix4& Matrix4::operator +=(const Matrix4& mat) {
    data[0] += mat[0];
    data[1] += mat[1];
    data[2] += mat[2];
    data[3] += mat[3];
    data[4] += mat[4];
    data[5] += mat[5];
    data[6] += mat[6];
    data[7] += mat[7];
    data[8] += mat[8];
    data[9] += mat[9];
    data[10] += mat[10];
    data[11] += mat[11];
    data[12] += mat[12];
    data[13] += mat[13];
    data[14] += mat[14];
    data[15] += mat[15];

    return *this;
}


Matrix4& Matrix4::operator -=(const Matrix4& mat) {
    data[0] -= mat[0];
    data[1] -= mat[1];
    data[2] -= mat[2];
    data[3] -= mat[3];
    data[4] -= mat[4];
    data[5] -= mat[5];
    data[6] -= mat[6];
    data[7] -= mat[7];
    data[8] -= mat[8];
    data[9] -= mat[9];
    data[10] -= mat[10];
    data[11] -= mat[11];
    data[12] -= mat[12];
    data[13] -= mat[13];
    data[14] -= mat[14];
    data[15] -= mat[15];

    return *this;
}


Matrix4& Matrix4::operator *=(const Matrix4& mat)
{
    (*this) = (*this) * mat;
    return *this;
}


Matrix4& Matrix4::operator *=(float value)
{
    data[0] *= value;
    data[1] *= value;
    data[2] *= value;
    data[3] *= value;
    data[4] *= value;
    data[5] *= value;
    data[6] *= value;
    data[7] *= value;
    data[8] *= value;
    data[9] *= value;
    data[10] *= value;
    data[11] *= value;
    data[12] *= value;
    data[13] *= value;
    data[14] *= value;
    data[15] *= value;

    return *this;
}


Matrix4& Matrix4::operator /=(float value)
{
    (*this) = (*this) / value;
    return *this;
}


bool Matrix4::operator ==(const Matrix4& mat) const
{
    return (data[0] == mat.data[0] &&
            data[1] == mat.data[1] &&
            data[2] == mat.data[2] &&
            data[3] == mat.data[3] &&
            data[4] == mat.data[4] &&
            data[5] == mat.data[5] &&
            data[6] == mat.data[6] &&
            data[7] == mat.data[7] &&
            data[8] == mat.data[8] &&
            data[9] == mat.data[9] &&
            data[10] == mat.data[10] &&
            data[11] == mat.data[11] &&
            data[12] == mat.data[12] &&
            data[13] == mat.data[13] &&
            data[14] == mat.data[14] &&
            data[15] == mat.data[15]
            );
}


bool Matrix4::operator !=(const Matrix4& mat) const
{
    return !((*this) == mat);
}


Matrix4 Matrix4::operator -() const
{
    Matrix4 result(*this);

    result.data[0] *= -1.0f;
    result.data[1] *= -1.0f;
    result.data[2] *= -1.0f;
    result.data[3] *= -1.0f;

    result.data[4] *= -1.0f;
    result.data[5] *= -1.0f;
    result.data[6] *= -1.0f;
    result.data[7] *= -1.0f;

    result.data[8] *= -1.0f;
    result.data[9] *= -1.0f;
    result.data[10] *= -1.0f;
    result.data[11] *= -1.0f;

    result.data[12] *= -1.0f;
    result.data[13] *= -1.0f;
    result.data[14] *= -1.0f;
    result.data[15] *= -1.0f;

    return result;
}



Vector4 Matrix4::operator *(const Vector4& vec) const
{
    // if bottom row is (0, 0, 0, 1)
    if(data[3] == 0.0f && data[7] == 0.0f && data[11] == 0.0f && data[15] == 1.0f)
    {
        return Vector4(
                    data[0] * vec.x
                + data[4] * vec.y
                + data[8] * vec.z
                + data[12]* vec.w,

                data[1] * vec.x
                + data[5] * vec.y
                + data[9] * vec.z
                + data[13]* vec.w,

                data[2] * vec.x
                + data[6] * vec.y
                + data[10]* vec.z
                + data[14]* vec.w,

                vec.w
                );
    }

    return Vector4(
                data[0] * vec.x
            + data[4] * vec.y
            + data[8] * vec.z
            + data[12]* vec.w,

            data[1] * vec.x
            + data[5] * vec.y
            + data[9] * vec.z
            + data[13]* vec.w,

            data[2] * vec.x
            + data[6] * vec.y
            + data[10]* vec.z
            + data[14]* vec.w,

            data[3] * vec.x
            + data[7] * vec.y
            + data[11]* vec.z
            + data[15]* vec.w
            );
}



Vector3 Matrix4::operator *(const Vector3& vec) const
{
    return Vector3(
                data[0] * vec.x
            + data[4] * vec.y
            + data[8] * vec.z
            + data[12],

            data[1] * vec.x
            + data[5] * vec.y
            + data[9] * vec.z
            + data[13],

            data[2] * vec.x
            + data[6] * vec.y
            + data[10]* vec.z
            + data[14]
            );
}


Vector2 Matrix4::multiplied(const Vector2& vec) const
{
    return Vector2(
                data[0] * vec.x
            + data[4] * vec.y
            ,
            data[1] * vec.x
            + data[5] * vec.y
            );
}
float Matrix4::determinant() const
{
    float a = m22 * (m33 * m44 - m43 * m34) - m23 * (m32 * m44 - m34 * m42) + m24 * (m32 * m43 - m33 * m42);
    float b = m21 * (m33 * m44 - m43 * m34) - m23 * (m31 * m44 - m34 * m41) + m24* (m31 * m43 - m33 * m41);
    float c = m21 * (m32 * m44 - m34 * m42) - m22 * (m31 * m44 - m34 * m41) + m24 * (m31 * m42 - m32 * m41);
    float d = m21 * (m32 * m43 - m33 * m42) - m22 * (m31 * m43 - m33 * m41) + m23 * (m31 * m42 - m32 * m41);
    return m11 * a - m12 * b + m13 * c - m14 * d; //cofactors(a, b, c, d)
}



void Matrix4::scale(const Vector3& scaleFactor)
{
    data[0] *= scaleFactor.x;
    data[1] *= scaleFactor.x;
    data[2] *= scaleFactor.x;
    data[3] *= scaleFactor.x;
    data[4] *= scaleFactor.y;
    data[5] *= scaleFactor.y;
    data[6] *= scaleFactor.y;
    data[7] *= scaleFactor.y;
    data[8] *= scaleFactor.z;
    data[9] *= scaleFactor.z;
    data[10] *= scaleFactor.z;
    data[11] *= scaleFactor.z;
}




void Matrix4::rotate(const Quaternion& quat)
{
    Matrix4 matrix;
    FromQuaternion(quat, matrix);

    (*this) = (*this) * matrix;
}


void Matrix4::rotate(const Vector3& axis, float angle)
{
    Matrix4 matrix;
    matrix.set_rotation_axis(angle, axis);
    (*this) = (*this) * matrix;
}


void Matrix4::translate(const Vector3& translate)
{
    Matrix4 matrix = Matrix4::Identity;
    matrix.set_translation_part(translate);
    (*this) *= matrix;

}


void Matrix4::invert()
{
    *this = inverted();
}


void Matrix4::transpose()
{
    *this = transposed();
}


void Matrix4::invert_transpose()
{
    *this = inverse_transposed();
}


void Matrix4::affine_invert()
{
    (*this) = affine_inverted();
}


void Matrix4::affine_invert_transpose()
{
    (*this) = affine_inverse_transposed();
}


Matrix4 Matrix4::rotated(const Quaternion& quat) const
{
    Matrix4 matrix;
    FromQuaternion(quat, matrix);

    return (*this) * matrix;
}


Matrix4 Matrix4::rotated(const Vector3& axis, float angle) const
{
    Matrix4 matrix;
    matrix.set_rotation_axis(angle, axis);
    return (*this) * matrix;
}


Matrix4  Matrix4::translated(const Vector3& vec) const
{
    Matrix4 matrix = Matrix4::Identity;
    matrix.set_translation_part(vec);
    return ((*this) * matrix);

}


Matrix4 Matrix4::inverted() const
{
    float det = determinant();
    float invDet = 1.0f / det;

    Matrix4 invMat;

    if(M::Abs(det) > M::EPSILON)
    {
        invMat.m11 =  (m22 * (m33 * m44 - m43 * m34) - m23 * (m32 * m44 - m34 * m42) + m24 * (m32 * m43 - m33 * m42)) * invDet;
        invMat.m21 = -(m21 * (m33 * m44 - m43 * m34) - m23 * (m31 * m44 - m34 * m41) + m24* (m31 * m43 - m33 * m41)) * invDet;
        invMat.m31 =  (m21 * (m32 * m44 - m34 * m42) - m22 * (m31 * m44 - m34 * m41) + m24 * (m31 * m42 - m32 * m41)) * invDet;
        invMat.m41 = -(m21 * (m32 * m43 - m33 * m42) - m22 * (m31 * m43 - m33 * m41) + m23 * (m31 * m42 - m32 * m41)) * invDet;

        invMat.m12 = -(m12 * (m33 * m44 - m34 * m43) - m13 * (m32 * m44 - m34 * m42) + m14 * (m32 * m43 - m33 *  m42)) * invDet;
        invMat.m22 =  (m11 * (m33 * m44 - m34 * m43) - m13 * (m31 * m44 - m34 * m41) + m14 * (m31 * m43 - m33 *  m41)) * invDet;
        invMat.m32 = -(m11 * (m32 * m44 - m34 * m42) - m12 * (m31 * m44 - m34 * m41) + m14 * (m31 * m42 - m32 *  m41)) * invDet;
        invMat.m42 =  (m11 * (m32 * m43 - m33 * m42) - m12 * (m31 * m43 - m33 * m41) + m13 * (m31 * m42 - m32 *  m41)) * invDet;

        invMat.m13 =  (m12 * (m23 * m44 - m24 * m43) - m13 * (m22 * m44 - m24 * m42) + m14 * (m22 * m43 - m23 *  m42)) * invDet;
        invMat.m23 = -(m11 * (m23 * m44 - m24 * m43) - m13 * (m21 * m44 - m24 * m41) + m14 * (m21 * m43 - m23 *  m41)) * invDet;
        invMat.m33 =  (m11 * (m22 * m44 - m24 * m42) - m12 * (m21 * m44 - m24 * m41) + m14 * (m21 * m42 - m22 *  m41)) * invDet;
        invMat.m43 = -(m11 * (m22 * m43 - m23 * m42) - m12 * (m21 * m43 - m23 * m41) + m13 * (m21 * m42 - m22 *  m41)) * invDet;

        invMat.m14 = -(m12 * (m23 * m34 - m24 * m33) - m13 * (m22 * m34 - m24 * m32) + m14 * (m22 * m33 - m23 *  m32)) * invDet;
        invMat.m24 =  (m11 * (m23 * m34 - m24 * m33) - m13 * (m21 * m34 - m24 * m31) + m14 * (m21 * m33 - m23 *  m31)) * invDet;
        invMat.m34 = -(m11 * (m22 * m34 - m24 * m32) - m12 * (m21 * m34 - m24 * m31) + m14 * (m21 * m32 - m22 *  m31)) * invDet;
        invMat.m44 =  (m11 * (m22 * m33 - m23 * m32) - m12 * (m21 * m33 - m23 * m31) + m13 * (m21 * m32 - m22 *  m31)) * invDet;
    }
    return invMat;
}


Matrix4 Matrix4::transposed() const
{
    return Matrix4{
        data[ 0], data[ 4], data[ 8], data[12],
                data[ 1], data[ 5], data[ 9], data[13],
                data[ 2], data[ 6], data[10], data[14],
                data[ 3], data[ 7], data[11], data[15]
    };
}


Matrix4 Matrix4::inverse_transposed() const
{
    Matrix4 result;

    float tmp[12]; // pair storage
    float det; // determinant

    // calculate pairs for first 8 elements (cofactors)
    tmp[0] = data[10] * data[15];
    tmp[1] = data[11] * data[14];
    tmp[2] = data[9] * data[15];
    tmp[3] = data[11] * data[13];
    tmp[4] = data[9] * data[14];
    tmp[5] = data[10] * data[13];
    tmp[6] = data[8] * data[15];
    tmp[7] = data[11] * data[12];
    tmp[8] = data[8] * data[14];
    tmp[9] = data[10] * data[12];
    tmp[10] = data[8] * data[13];
    tmp[11] = data[9] * data[12];

    // calculate first 8 elements (cofactors)
    result[0] =	tmp[0]*data[5] + tmp[3]*data[6] + tmp[4]*data[7]
            -	tmp[1]*data[5] - tmp[2]*data[6] - tmp[5]*data[7];

    result[1] =	tmp[1]*data[4] + tmp[6]*data[6] + tmp[9]*data[7]
            -	tmp[0]*data[4] - tmp[7]*data[6] - tmp[8]*data[7];

    result[2] =	tmp[2]*data[4] + tmp[7]*data[5] + tmp[10]*data[7]
            -	tmp[3]*data[4] - tmp[6]*data[5] - tmp[11]*data[7];

    result[3] =	tmp[5]*data[4] + tmp[8]*data[5] + tmp[11]*data[6]
            -	tmp[4]*data[4] - tmp[9]*data[5] - tmp[10]*data[6];

    result[4] =	tmp[1]*data[1] + tmp[2]*data[2] + tmp[5]*data[3]
            -	tmp[0]*data[1] - tmp[3]*data[2] - tmp[4]*data[3];

    result[5] =	tmp[0]*data[0] + tmp[7]*data[2] + tmp[8]*data[3]
            -	tmp[1]*data[0] - tmp[6]*data[2] - tmp[9]*data[3];

    result[6] =	tmp[3]*data[0] + tmp[6]*data[1] + tmp[11]*data[3]
            -	tmp[2]*data[0] - tmp[7]*data[1] - tmp[10]*data[3];

    result[7] =	tmp[4]*data[0] + tmp[9]*data[1] + tmp[10]*data[2]
            -	tmp[5]*data[0] - tmp[8]*data[1] - tmp[11]*data[2];

    // calculate pairs for second 8 elements (cofactors)
    tmp[0] = data[2]*data[7];
    tmp[1] = data[3]*data[6];
    tmp[2] = data[1]*data[7];
    tmp[3] = data[3]*data[5];
    tmp[4] = data[1]*data[6];
    tmp[5] = data[2]*data[5];
    tmp[6] = data[0]*data[7];
    tmp[7] = data[3]*data[4];
    tmp[8] = data[0]*data[6];
    tmp[9] = data[2]*data[4];
    tmp[10] = data[0]*data[5];
    tmp[11] = data[1]*data[4];

    // calculate second 8 elements (cofactors)
    result[8] =	tmp[0]*data[13] + tmp[3]*data[14] + tmp[4]*data[15]
            -	tmp[1]*data[13] - tmp[2]*data[14] - tmp[5]*data[15];

    result[9] =	tmp[1]*data[12] + tmp[6]*data[14] + tmp[9]*data[15]
            -	tmp[0]*data[12] - tmp[7]*data[14] - tmp[8]*data[15];

    result[10] =	tmp[2]*data[12] + tmp[7]*data[13] + tmp[10]*data[15]
            -	tmp[3]*data[12] - tmp[6]*data[13] - tmp[11]*data[15];

    result[11] =	tmp[5]*data[12] + tmp[8]*data[13] + tmp[11]*data[14]
            -	tmp[4]*data[12] - tmp[9]*data[13] - tmp[10]*data[14];

    result[12] =	tmp[2]*data[10] + tmp[5]*data[11] + tmp[1]*data[9]
            -	tmp[4]*data[11] - tmp[0]*data[9] - tmp[3]*data[10];

    result[13] =	tmp[8]*data[11] + tmp[0]*data[8] + tmp[7]*data[10]
            -	tmp[6]*data[10] - tmp[9]*data[11] - tmp[1]*data[8];

    result[14] =	tmp[6]*data[9] + tmp[11]*data[11] + tmp[3]*data[8]
            -	tmp[10]*data[11] - tmp[2]*data[8] - tmp[7]*data[9];

    result[15] =	tmp[10]*data[10] + tmp[4]*data[8] + tmp[9]*data[9]
            -	tmp[8]*data[9] - tmp[11]*data[10] - tmp[5]*data[8];

    // calculate determinant
    det	= data[0]*result[0]
            + data[1]*result[1]
            + data[2]*result[2]
            + data[3]*result[3];

    if(det == 0.0f)
    {
        Matrix4 id;
        return id;
    }

    result = result/det;

    return result;
}


Matrix4 Matrix4::affine_inverted() const
{
    return Matrix4{
        data[0],
                data[4],
                data[8],
                0.0f,
                data[1],
                data[5],
                data[9],
                0.0f,
                data[2],
                data[6],
                data[10],
                0.0f,
                -(data[0]*data[12] + data[1]*data[13] + data[2]*data[14]),
                -(data[4]*data[12] + data[5]*data[13] + data[6]*data[14]),
                -(data[8]*data[12] + data[9]*data[13] + data[10]*data[14]),
                1.0f
    };
}


Matrix4 Matrix4::affine_inverse_transposed() const
{
    return Matrix4{
        data[0],
                data[1],
                data[2],
                -(data[0]*data[12] + data[1]*data[13] + data[2]*data[14]),
                data[4],
                data[5],
                data[6],
                -(data[4]*data[12] +data[5]*data[13] + data[6]*data[14]),
                data[8],
                data[9],
                data[10],
                -(data[8]*data[12] +data[9]*data[13] + data[10]*data[14]),
                0.0f, 0.0f, 0.0f, 1.0f
    };
}


Vector3 Matrix4::rotated_vector(const Vector3& vec) const
{
    return Vector3(
                data[0]*vec.x + data[4]*vec.y + data[8]*vec.z,
            data[1]*vec.x + data[5]*vec.y + data[9]*vec.z,
            data[2]*vec.x + data[6]*vec.y + data[10]*vec.z
            );
}


Vector3 Matrix4::inverse_rotated_vector(const Vector3& vec) const
{
    return Vector3(
                data[0]*vec.x + data[1]*vec.y + data[2]*vec.z,
            data[4]*vec.x + data[5]*vec.y + data[6]*vec.z,
            data[8]*vec.x + data[9]*vec.y + data[10]*vec.z
            );
}



void Matrix4::set_scale(const Vector3& scaleFactor)
{
    *this = Matrix4::Identity;

    data[0] = scaleFactor.x;
    data[5] = scaleFactor.y;
    data[10] = scaleFactor.z;
}


void Matrix4::set_scale(float size)
{
    *this = Matrix4::Identity;
    data[0] = data[5] = data[10] = size;
}


void Matrix4::set_translation(const Vector3& translation)
{
    *this = Matrix4::Identity;
    set_translation_part(translation);
}


void Matrix4::set_translation_part(const Vector3& translation)
{
    data[12] = translation.x;
    data[13] = translation.y;
    data[14] = translation.z;
}

Vector3 Matrix4::get_translation_part() const {
    return Vector3(data[12], data[13], data[14]);
}




void Matrix4::set_rotation_axis(float angle, const Vector3& axis)
{
    Vector3 u = axis.normalized();

    float sinAngle = (float)sin(angle /** DEG_TO_RAD*/);
    float cosAngle = (float)cos(angle /** DEG_TO_RAD*/);
    float oneMinusCosAngle = 1.0f - cosAngle;

    *this = Matrix4::Identity;

    data[0] = (u.x)*(u.x) + cosAngle * (1.0f-(u.x)*(u.x));
    data[4] = (u.x)*(u.y)*(oneMinusCosAngle) - sinAngle*u.z;
    data[8] = (u.x)*(u.z)*(oneMinusCosAngle) + sinAngle*u.y;

    data[1] = (u.x)*(u.y)*(oneMinusCosAngle) + sinAngle*u.z;
    data[5] = (u.y)*(u.y) + cosAngle * (1.0f-(u.y)*(u.y));
    data[9] = (u.y)*(u.z)*(oneMinusCosAngle) - sinAngle*u.x;

    data[2] = (u.x)*(u.z)*(oneMinusCosAngle) - sinAngle*u.y;
    data[6] = (u.y)*(u.z)*(oneMinusCosAngle) + sinAngle*u.x;
    data[10] = (u.z)*(u.z) + cosAngle * (1.0f-(u.z)*(u.z));
}


void Matrix4::set_rotation_x(float angle)
{
    *this = Matrix4::Identity;

    data[5] = (float)cos(angle * M::DEG_TO_RAD);
    data[6] = (float)sin(angle * M::DEG_TO_RAD);

    data[9] = -data[6];
    data[10] = data[5];
}


void Matrix4::set_rotation_y(float angle)
{
    *this = Matrix4::Identity;

    data[0] =  (float)cos(angle * M::DEG_TO_RAD);
    data[2] = -(float)sin(angle * M::DEG_TO_RAD);

    data[2] = -data[8];
    data[10] = data[0];
}


void Matrix4::set_rotation_z(float angle)
{
    *this = Matrix4::Identity;

    data[0] = (float)cos(angle * M::DEG_TO_RAD);
    data[1] = (float)sin(angle * M::DEG_TO_RAD);

    data[4] = -data[1];
    data[5] =  data[0];
}


void Matrix4::set_rotation_euler(float angleX, float angleY, float angleZ)
{
    *this = Matrix4::Identity;
    set_rotation_part_euler(angleX, angleY, angleZ);
}


void Matrix4::set_rotation_part_euler(float angleX, float angleY, float angleZ)
{
    double cr = cos( angleX * M::DEG_TO_RAD );
    double sr = sin( angleX * M::DEG_TO_RAD );
    double cp = cos( angleY * M::DEG_TO_RAD );
    double sp = sin( angleY * M::DEG_TO_RAD );
    double cy = cos( angleZ * M::DEG_TO_RAD );
    double sy = sin( angleZ * M::DEG_TO_RAD );

    data[0] = (float)(cp*cy);
    data[1] = (float)(cp*sy);
    data[2] = (float)(-sp );

    double srsp = sr*sp;
    double crsp = cr*sp;

    data[4] = (float)(srsp*cy-cr*sy);
    data[5] = (float)(srsp*sy+cr*cy);
    data[6] = (float)(sr*cp );

    data[8] = (float)(crsp*cy+sr*sy);
    data[9] = (float)(crsp*sy-sr*cy);
    data[10] = (float)(cr*cp);
}


Vector3 Matrix4::get_euler_angles() const
{
    //Get euler angles from matrix.
    float cy = sqrtf(data[0] * data[0] + data[1] * data[1]);
    if(cy > (16.0 * 1.192092896e-07F))
    {
        Vector3 euler;
        Vector3 euler1;
        Vector3 euler2;

        euler1.x = atan2f(data[6], data[10]);
        euler1.y = atan2f(- data[2], cy);
        euler1.z = atan2f(data[1], data[0]);

        euler2.x = atan2f(- data[6], - data[10]);
        euler2.y = atan2f(- data[2], - cy);
        euler2.z = atan2f(- data[1], - data[0]);

        if( (fabs(euler1.x) + fabs(euler1.y) + fabs(euler1.z)) >
                (fabs(euler2.x) + fabs(euler2.y) + fabs(euler2.z)))
        {
            euler = euler2;
        }
        else
        {
            euler = euler1;
        }

        return euler * M::RAD_TO_DEG;
    }
    else
    {
        Vector3 euler;

        euler.x = atan2f(- data[9], data[5]);
        euler.y = atan2f(- data[2], cy);
        euler.z = 0.0f;

        return euler * (float)M::RAD_TO_DEG;
    }

    return Vector3(0, 0, 0);
}
/*

Quaternion Matrix4::getQuaternion() const
{
    Quaternion q;
     // Algorithm in Ken Shoemake's article in 1987 SIGGRAPH course notes
                    // article "Quaternion Calculus and Fast Animation".
    Log()<<"TEST\n";

    const float trace = data2[0][0] + data2[1][1] + data2[2][2];

    if (trace>0)
    {
        // |w| > 1/2, may as well choose w > 1/2
        Log()<<"TEST\n";
        float root = M::Sqrt(trace + 1.0f);  // 2w
        q.w = 0.5f * root;
        root = 0.5f / root;  // 1/(4w)
        q.x = (data2[2][1]-data2[1][2]) * root;
        q.y = (data2[0][2]-data2[2][0]) * root;
        q.z = (data2[1][0]-data2[0][1]) * root;
    }
    else
    {
        // |w| = 1/2
        Log()<<"TEST\n";

        static int next[3] = { 1, 2, 0 };

        int i = 0;
        if (data2[1][1]>data2[0][0])  i = 1;
        if (data2[2][2]>data2[i][i]) i = 2;
        int j = next[i];
        int k = next[j];

        float root = M::Sqrt(data2[i][i]-data2[j][j]-data2[k][k] + 1.0f);
        float *quaternion[3] = { &q.x, &q.y, &q.z };
        *quaternion[i] = 0.5f * root;
        root = 0.5f / root;
        q.w = (data2[k][j]-data2[j][k])*root;
        *quaternion[j] = (data2[j][i]+data2[i][j])*root;
        *quaternion[k] = (data2[k][i]+data2[i][k])*root;
    }

    return q;

    const float diag = data2[0][0] + data2[1][1]+ data2[2][2] + 1;

        if( diag > 0.0f )
        {
                const float scale = sqrtf(diag) * 2.0f; // get scale from diagonal

                // TODO: speed this up
                q.x = ( data2[2][1] - data2[1][2]) / scale;
                q.y = ( data2[0][2] - data2[2][0]) / scale;
                q.z = ( data2[1][0] - data2[0][1]) / scale;
                q.w = 0.25f * scale;
        }
        else
        {
                if ( data2[0][0]> data2[1][1] && data2[0][0] > data2[2][2])
                {
                        // 1st element of diag is greatest value
                        // find scale according to 1st element, and double it
                        const float scale = sqrtf( 1.0f + data2[0][0] - data2[1][1] - data2[2][2]) * 2.0f;

                        // TODO: speed this up
                        q.x = 0.25f * scale;
                        q.y = (data2[0][1] + data2[1][0]) / scale;
                        q.z = (data2[2][0] + data2[0][2]) / scale;
                        q.w = (data2[2][1] - data2[1][2]) / scale;
                }
                else if ( data2[1][1] > data2[2][2])
                {
                        // 2nd element of diag is greatest value
                        // find scale according to 2nd element, and double it
                        const float scale = sqrtf( 1.0f + data2[1][1] - data2[0][0]- data2[2][2]) * 2.0f;

                        // TODO: speed this up
                        q.x = (data2[0][1]+ data2[1][0] ) / scale;
                        q.y = 0.25f * scale;
                        q.z = (data2[1][2] + data2[2][1] ) / scale;
                        q.w = (data2[0][2] - data2[2][0] ) / scale;
                }
                else
                {
                        // 3rd element of diag is greatest value
                        // find scale according to 3rd element, and double it
                        const float scale = sqrtf( 1.0f + data2[2][2] - data2[0][0] - data2[1][1]) * 2.0f;

                        // TODO: speed this up
                        q.x = (data2[0][2] + data2[2][0]) / scale;
                        q.y = (data2[1][2] + data2[2][1]) / scale;
                        q.z = 0.25f * scale;
                        q.w = (data2[1][0] - data2[0][1]) / scale;
                }
        }
    q.normalize();
        return q;

    // Algorithm in Ken Shoemake's article in 1987 SIGGRAPH course notes
    // article "Quaternion Calculus and Fast Animation".

    float trace = data2[0][0]+ data2[1][1]+ data2[2][2];
    float root;

    if ( trace > 0.0 )
    {
        // |w| > 1/2, may as well choose w > 1/2
        root = M::Sqrt(trace + 1.0f);  // 2w
        q.w = 0.5f*root;
        root = 0.5f/root;  // 1/(4w)
        q.x = (data2[2][1]-data2[1][2])*root;
        q.y = (data2[0][2]-data2[2][0])*root;
        q.z = (data2[1][0]-data2[0][1])*root;
    }
    else
    {
        // |w| <= 1/2
        static size_t s_iNext[3] = { 1, 2, 0 };
        size_t i = 0;
        if ( data2[1][1] > data2[0][0] )
            i = 1;
        if ( data2[2][2] > data2[i][i] )
            i = 2;
        size_t j = s_iNext[i];
        size_t k = s_iNext[j];

        root = M::Sqrt(data2[i][i]-data2[j][j]-data2[k][k] + 1.0f);
        float* apkQuat[3] = { &q.x, &q.y, &q.z };
        *apkQuat[i] = 0.5f*root;
        root = 0.5f/root;
        q.w = (data2[k][j]-data2[j][k])*root;
        *apkQuat[j] = (data2[j][i]+data2[i][j])*root;
        *apkQuat[k] = (data2[k][i]+data2[i][k])*root;
    }
    return q;

    //Compute one plus the trace of the matrix.
    const double onePlusTrace = 1.0 + data2[0][0] + data2[1][1] + data2[2][2];

    //Check the diagonal.
    if (onePlusTrace > 1E-5)
    {
      //Direct computation.
        const double s = sqrt(onePlusTrace) * 2.0;
        Quaternion q((data2[2][1] - data2[1][2]) / s,(data2[0][2] - data2[2][0]) / s,(data2[1][0] - data2[0][1]) / s,0.25 * s);
    }
    else
    {
        //Computation depends on major diagonal term.
        if ((data2[0][0] > data2[1][1])&(data2[0][0] > data2[2][2]))
        {
          const double s = sqrt(1.0 + data2[0][0] - data2[1][1] - data2[2][2]) * 2.0;
           q = Quaternion( 0.25 * s,(data2[0][1] + data2[1][0]) / s, (data2[0][2] + data2[2][0]) / s,(data2[1][2] - data2[2][1]) / s);
        }
        else if (data2[1][1] > data2[2][2])
        {
            const double s = sqrt(1.0 + data2[1][1] - data2[0][0] - data2[2][2]) * 2.0;
            q = Quaternion( (data2[0][1] + data2[1][0]) / s,
                            0.25 * s,
                            (data2[1][2] + data2[2][1]) / s,
                            (data2[0][2] - data2[2][0]) / s);
        }
        else
        {
            const double s = sqrt(1.0 + data2[2][2] - data2[0][0] - data2[1][1]) * 2.0;
            q = Quaternion( (data2[0][2] + data2[2][0]) / s,
                            (data2[1][2] + data2[2][1]) / s,
                            0.25 * s,
                            (data2[0][1] - data2[1][0]) / s);
        }
    }
    //Set inverse.
    return q.inversed();
}*/

// WORKS !

Matrix4 Matrix4::ModelViewFromDirectionLH(const Vector3& eyePos, const Vector3& dir, const Vector3& up) {

    Vector3 xAxis, yAxis, zAxis;

    zAxis = dir.normalized();
    xAxis = up.cross(zAxis).normalized();
    yAxis = zAxis.cross(xAxis);

    Matrix4 m;
    float4x4& data = m.data4x4;
    data[0][0] = xAxis.x;
    data[1][0] = xAxis.y;
    data[2][0] = xAxis.z;

    data[0][1] = yAxis.x;
    data[1][1] = yAxis.y;
    data[2][1] = yAxis.z;

    data[0][2] = zAxis.x;
    data[1][2] = zAxis.y;
    data[2][2] = zAxis.z;

    data[3][0] = -xAxis.dot(eyePos);
    data[3][1] = -yAxis.dot(eyePos);
    data[3][2] = -zAxis.dot(eyePos);

    data[0][3] = 0.0f;
    data[1][3] = 0.0f;
    data[2][3] = 0.0f;

    m[15] = 1.0;

    return m;

}

// Works !

Matrix4 Matrix4::ModelViewFromLookAtLH(const Vector3& eyePos, const Vector3& targetPos, const Vector3& up) {
    return ModelViewFromDirectionLH(eyePos, targetPos - eyePos, up);
}

// Works !

Matrix4 Matrix4::ModelViewFromAxisLH(const Vector3& pos, const Vector3& xAxis, const Vector3& yAxis, const Vector3& zAxis) {

    Matrix4 m;
    m[0] = xAxis.x;
    m[1] = yAxis.x;
    m[2] = zAxis.x;
    m[3] = 0.0;

    m[4] = xAxis.y;
    m[5] = yAxis.y;
    m[6] = zAxis.y;
    m[7] = 0.0;

    m[8] = xAxis.z;
    m[9] = yAxis.z;
    m[10] = zAxis.z;
    m[11] = 0.0;

    m[12] = -xAxis.dot(pos);
    m[13] = -yAxis.dot(pos);
    m[14] = -zAxis.dot(pos);
    m[15] = 1.0f;

    return m;
}

Matrix4 Matrix4::OrthoLH(const Vector4& corners, const Vector2& zNearFar) {


    Matrix4 m = Matrix4::Identity;

    const float top = corners.top;
    const float left = corners.left;
    const float bottom = corners.bottom;
    const float right = corners.right;


    const float zNear = zNearFar.x;
    const float zFar = zNearFar.y;

    const float inv_right_minus_left = 1.0f / (right - left);
    const float inv_top_minus_bottom = 1.0f / (top - bottom);

    const float inv_far_minus_near = 1.0f / (zFar - zNear);


    m.data4x4[0][0] = 2.0f * inv_right_minus_left;
    m.data4x4[1][1] = 2.0f * inv_top_minus_bottom;

    m.data4x4[3][0] = - (right + left) * inv_right_minus_left;
    m.data4x4[3][1] = - (top + bottom) * inv_top_minus_bottom;

    m.data4x4[2][2] = 2.0f * inv_far_minus_near;
    m.data4x4[3][2] = - (zFar + zNear) * inv_far_minus_near;

    return m;
}

Matrix4 Matrix4::OrthoLH(const Vector2& min, const Vector2& max, const Vector2& zNearFar, bool flip_y_axis) {


    Matrix4 m = Matrix4::Identity;

    const float left = min.x;
    const float right = max.x;
    const float bottom = min.y;
    const float top = max.y;



    const float zNear = zNearFar.x;
    const float zFar = zNearFar.y;

    const float inv_right_minus_left = 1.0f / (right - left);
    const float inv_top_minus_bottom = 1.0f / (top - bottom);

    const float inv_far_minus_near = 1.0f / (zFar - zNear);


    m.data4x4[0][0] = 2.0f * inv_right_minus_left;
    m.data4x4[1][1] = 2.0f * inv_top_minus_bottom;

    m.data4x4[3][0] = - (right + left) * inv_right_minus_left;
    m.data4x4[3][1] = - (top + bottom) * inv_top_minus_bottom;

    //if (flip_y_axis) {
    //m.data4x4[3][1] = -m.data4x4[3][1];
    // }

    m.data4x4[2][2] = 2.0f * inv_far_minus_near;
    m.data4x4[3][2] = - (zFar + zNear) * inv_far_minus_near;

    return m;
}

Matrix4 Matrix4::OrthoLHB(const Vector2& min, const Vector2& max, const Vector2& zNearFar, bool flip_y_axis) {

    const float& l = min.x;
    const float& r = max.x;
    const float& b = min.y;
    const float& t = max.y;

    const float& n = zNearFar.x;
    const float& f = zNearFar.y;

    Matrix4 m;
    m.data4x4[0][0] = 2.f/(r-l);
    m.data4x4[0][1] = m.data4x4[0][2] = m.data4x4[0][3] = 0.f;

    m.data4x4[1][1] = 2.f/(t-b);
    m.data4x4[1][0] = m.data4x4[1][2] = m.data4x4[1][3] = 0.f;

    m.data4x4[2][2] = -2.f/(f-n);
    m.data4x4[2][0] = m.data4x4[2][1] = m.data4x4[2][3] = 0.f;

    m.data4x4[3][0] = -(r+l)/(r-l);
    m.data4x4[3][1] = -(t+b)/(t-b);
    m.data4x4[3][2] = -(f+n)/(f-n);
    m.data4x4[3][3] = 1.f;
    return m;
}



//matrix4 Matrix4::ImGuiOrthoLH(const Vector2& zNearFar, const Vector2& vpScaling, float screen_w, float screen_h) {

//    const float ratio = screen_h / screen_w;

//    Vector4 corners;
//    corners.top = 0;
//    corners.right = 1.0 / ((screen_h) * (1.0f/ (vpScaling.x * 2.0f)));
//    corners.bottom = -1.0 /((screen_w * ratio ) * (1.0f/ (vpScaling.y * 2.0f)));
//    corners.left = 0.0f;
//    return OrthoLH(corners, zNearFar);
//}


/*
Matrix4 Matrix4::FromRows(const Vector3& row1, const Vector3& row2, const Vector3& row3) {
    Matrix4 m;
    m[0] = row1.x;
    m[1] = row1.x;
    m[2] = row1.x;
    m[3] = 0.0;

    m[4] = row2.y;
    m[5] = row2.y;
    m[6] = row2.y;
    m[7] = 0.0;

    m[8] = row3.z;
    m[9] = row3.z;
    m[10] = row3.z;
    m[11] = 0.0;

    m[12] = 0.0f;
    m[13] = 0.0f;
    m[14] = 0.0f;
    m[15] = 1.0f;

    return m;
}

Matrix4 Matrix4::FromColumns(const Vector3& col1, const Vector3& col2, const Vector3& col3) {
    Matrix4 m;
    m[0] = col1.x;
    m[1] = col2.x;
    m[2] = col3.x;
    m[3] = 0.0f;

    m[4] = col1.y;
    m[5] = col2.y;
    m[6] = col3.y;
    m[7] = 0.0f;

    m[8] = col1.z;
    m[9] = col2.z;
    m[10] = col3.z;
    m[11] = 0.0f;

    m[12] = 0.0f;
    m[13] = 0.0f;
    m[14] = 0.0f;
    m[15] = 1.0f;

    return m;
}
*/

void Matrix4::FromQuaternion(const Quaternion& quat, Matrix4& mat)
{

    mat[0] = 1.0f - 2.0f*quat.y*quat.y - 2.0f*quat.z*quat.z;
    mat[1] = 2.0f*quat.x*quat.y + 2.0f*quat.z*quat.w;
    mat[2] = 2.0f*quat.x*quat.z - 2.0f*quat.y*quat.w;
    mat[3] = 0.0f;

    mat[4] = 2.0f*quat.x*quat.y - 2.0f*quat.z*quat.w;
    mat[5] = 1.0f - 2.0f*quat.x*quat.x - 2.0f*quat.z*quat.z;
    mat[6] = 2.0f*quat.z*quat.y + 2.0f*quat.x*quat.w;
    mat[7] = 0.0f;

    mat[8] = 2.0f*quat.x*quat.z + 2.0f*quat.y*quat.w;
    mat[9] = 2.0f*quat.z*quat.y - 2.0f*quat.x*quat.w;
    mat[10] = 1.0f - 2.0f*quat.x*quat.x - 2.0f*quat.y*quat.y;
    mat[11] = 0.0f;

    mat[12] = 0.0f;
    mat[13] = 0.0f;
    mat[14] = 0.0f;
    mat[15] = 1.f;

    /*float fTX  = 2.0f * quat.x;
        float fTY  = 2.0f * quat.y;
        float fTZ  = 2.0f * quat.z;
        float fTWX = fTX * quat.w;
        float fTWY = fTY * quat.w;
        float fTWZ = fTZ * quat.w;
        float fTXX = fTX * quat.x;
        float fTXY = fTY * quat.x;
        float fTXZ = fTZ * quat.x;
        float fTYY = fTY * quat.y;
        float fTYZ = fTZ * quat.y;
        float fTZZ = fTZ * quat.z;

    //mat = Matrix4::Identity;

        mat[0]  = 1.0f - ( fTYY + fTZZ );
        mat[1]  = fTXY - fTWZ;
        mat[2]  = fTXZ + fTWY;

        mat[4]  = fTXY + fTWZ;
        mat[5]  = 1.0f - ( fTXX + fTZZ );
        mat[6]  = fTYZ - fTWX;

        mat[8]  = fTXZ - fTWY;
        mat[9]  = fTYZ + fTWX;
        mat[10] = 1.0f - ( fTXX + fTYY );

        mat[3] = mat[7] = mat[11] = mat[12] = mat[13] = mat[14] = 0;
        mat[15] = 1.0;*/
}



void Matrix4::FromQuaternion2(const Quaternion& quat, Matrix4& mat) {

    float fTX  = 2.0f * quat.x;
    float fTY  = 2.0f * quat.y;
    float fTZ  = 2.0f * quat.z;
    float fTWX = fTX * quat.w;
    float fTWY = fTY * quat.w;
    float fTWZ = fTZ * quat.w;
    float fTXX = fTX * quat.x;
    float fTXY = fTY * quat.x;
    float fTXZ = fTZ * quat.x;
    float fTYY = fTY * quat.y;
    float fTYZ = fTZ * quat.y;
    float fTZZ = fTZ * quat.z;

    mat[0] = 1.0f - fTYY - fTZZ;
    mat[1] = fTXY + fTWZ;
    mat[2] = fTXZ - fTWY;
    mat[3] = 0.0f;

    mat[4] = fTXY - fTWZ;
    mat[5] = 1.0f - fTXX - fTZZ;
    mat[6] = fTYZ + fTWX;
    mat[7] = 0.0f;

    mat[8] = fTXZ + fTWY;
    mat[9] = fTYZ - fTWX;
    mat[10] = 1.0f - fTXX - fTYY;
    mat[11] = 0.0f;

    mat[12] = 0.0f;
    mat[13] = 0.0f;
    mat[14] = 0.0f;
    mat[15] = 1.f;

    /*float fTX  = 2.0f * quat.x;
        float fTY  = 2.0f * quat.y;
        float fTZ  = 2.0f * quat.z;
        float fTWX = fTX * quat.w;
        float fTWY = fTY * quat.w;
        float fTWZ = fTZ * quat.w;
        float fTXX = fTX * quat.x;
        float fTXY = fTY * quat.x;
        float fTXZ = fTZ * quat.x;
        float fTYY = fTY * quat.y;
        float fTYZ = fTZ * quat.y;
        float fTZZ = fTZ * quat.z;

    //mat = Matrix4::Identity;

        mat[0]  = 1.0f - ( fTYY + fTZZ );
        mat[1]  = fTXY - fTWZ;
        mat[2]  = fTXZ + fTWY;

        mat[4]  = fTXY + fTWZ;
        mat[5]  = 1.0f - ( fTXX + fTZZ );
        mat[6]  = fTYZ - fTWX;

        mat[8]  = fTXZ - fTWY;
        mat[9]  = fTYZ + fTWX;
        mat[10] = 1.0f - ( fTXX + fTYY );

        mat[3] = mat[7] = mat[11] = mat[12] = mat[13] = mat[14] = 0;
        mat[15] = 1.0;*/
}
