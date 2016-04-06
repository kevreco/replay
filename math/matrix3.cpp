#include "matrix3.h"

#include "vector3.h"


#include <memory.h>
#include "maths.h"


const Matrix3 Matrix3::Identity = {
    1, 0, 0,
    0, 1, 0,
    0, 0, 1
};

/*
Matrix3::Matrix3(const Vector2& translation) {
    *this = Matrix3::Identity;
    this->translate(translation);
}*/

Matrix3 Matrix3::operator *(float value) const {
    if (value == 1.0f)
        return *this;

    return Matrix3{
        data[0] * value,
                data[1] * value,
                data[2] * value,
                data[3] * value,
                data[4] * value,
                data[5] * value,
                data[6] * value,
                data[7] * value,
                data[8] * value
    };
}


Matrix3 Matrix3::operator /(float value) const
{
    if(value == 0.0f || value == 1.0f)
        return (*this);

    return (*this) * (1.0f / value);
}

Matrix3& Matrix3::operator *=(float value) {
    if (value == 1.0)
        return *this;

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

    return *this;
}


Matrix3& Matrix3::operator /=(float value) {
    (*this) = (*this) / value;
    return *this;
}


bool Matrix3::operator ==(const Matrix3& mat) const {
    return (data[0] == mat.data[0] &&
            data[1] == mat.data[1] &&
            data[2] == mat.data[2] &&
            data[3] == mat.data[3] &&
            data[4] == mat.data[4] &&
            data[5] == mat.data[5] &&
            data[6] == mat.data[6] &&
            data[7] == mat.data[7] &&
            data[8] == mat.data[8]
            );
}

bool Matrix3::operator !=(const Matrix3& mat) const {
    return !((*this) == mat);
}

Matrix3 Matrix3::operator -() const {
    Matrix3 result(*this);

    result.data[0] *= -1.0f;
    result.data[1] *= -1.0f;
    result.data[2] *= -1.0f;

    result.data[3] *= -1.0f;
    result.data[4] *= -1.0f;
    result.data[5] *= -1.0f;

    result.data[6] *= -1.0f;
    result.data[7] *= -1.0f;
    result.data[8] *= -1.0f;

    return result;
}
/*

TVector2<T> Matrix3::operator *(const TVector2<T>& vec) const
{
        return TVector2<T>(
                  data[0] * vec.x
                + data[3] * vec.y
                + data[6],

                  data[1] * vec.x
                + data[4] * vec.y
                + data[7]
        );
}*/


void Matrix3::apply(float px, float py, float& rx, float& ry) const {
    rx = data[0] * px + data[3] * py + data[6];
    ry = data[1] * px + data[4] * py + data[7];
}

Vector2 Matrix3::apply(float px, float py) const {
    return Vector2(
                data[0] * px
            + data[3] * py
            + data[6],

            data[1] * px
            + data[4] * py
            + data[7]
            );
}


void Matrix3::toIdentity() {
    data[0] = 1;
    data[1] = 0;
    data[2] = 0;

    data[3] = 0;
    data[4] = 1;
    data[5] = 0;

    data[6] = 0;
    data[7] = 0;
    data[8] = 1;

}


float Matrix3::determinant() const {
    /*return m11*(m33 * m22 - m32 * m23) -
        m21*(m33 * m12 - m32* m13) + m31*(m23 * m12 - m22 * m13);*/
    // m13 and m23 == 0
    return (m11*m22*m33 /*+ m21*m32*m13*/ /*+ m31*m12*m23*/) -
            (m13*m22*m31 /*+ m23*m32*m11 + m33*m12*m21*/);
}


void Matrix3::setTranslation(float x, float y) {
    tx = x;
    ty = y;
}

Vector2 Matrix3::translation() const {
    return Vector2(tx, ty);
}


Vector2 Matrix3::scale() const {
    return Vector2(m11, m22);
}
/*
void Matrix3::rotate(float angle, bool inDegree) {
    (inDegree) ? rotateInDeg(angle): rotateInRad(angle);
}

void Matrix3::rotateInDeg(float angle) {
   rotateInRad(M::DegToRad(angle));
}


void Matrix3::invert() {
        *this = inverted();
}


Matrix3 Matrix3::inverted() const
{
        T det = determinant();
    T invDet = (float)1.0 / det;

    Matrix3 invMat;
        return invMat;
}


Transform3::Transform3(const Transform3& tr) {
    memcpy(data, tr.data, 9 * sizeof(float));
    mLvl = tr.mLvl;
}

Transform3 Transform3::operator *(float value) const {

    Matrix3::operator *(value);

    if (mLvl < Scaled)
        mLvl = Scaled;

    return  *this;
}

Transform3 Transform3::operator /(float value) const {
    if(value == 0.0f || value == 1.0f)
        return (*this);

    if (mLvl < Scaled)
        mLvl = Scaled;

    return (*this) * (1.0f / value);
}


Transform3& Transform3::operator *=(const Transform3& m) {
    if (m.isIdentity())
        return *this;

    if (isIdentity()) {
        *this = m;
        return *this;
    }


    if (m.mLvl > mLvl)
        mLvl = m.mLvl;

    switch(mLvl) {
        case Translated: {
            tx += m.tx;
            ty += m.ty;
            break;
        }
        case Scaled: {

            const T tmp11 = m11 * m.m11;
            const T tmp22 = m22 * m.m22;
            tx = m11 * m.tx + tx;
            ty = m22 * m.ty + ty;
            m11 = tmp11;
            m22 = tmp22;
            break;
        }
        case Rotated: {
            const float tmp11 = m11 * m.m11 + m21 * m.m12;
            const float tmp12 = m12 * m.m11 + m22 * m.m12;

            const float tmp21 = m11 * m.m21 + m.m22 * m21;
            const float tmp22 = m12 * m.m21 + m.m22 * m22;

            tx = m11 * m.m31 + m21 * m.m32 + m31;
            ty = m12 * m.m31 + m22 * m.m32 + m32;

            m11 = tmp11;
            m12 = tmp12;
            m21 = tmp21;
            m22 = tmp22;
            break;
        }
        default: break;
    }

    return *this;
}


Transform3& Transform3::operator *=(float value) {
    *this = Matrix3::operator *=(value);

    if (mLvl < Scaled)
        mLvl = Scaled;

    return *this;
}

bool Transform3::operator ==(const Transform3& tr) const {
    return (Matrix3::operator ==(*this, tr)
            && mLvl == tr.mLvl );
}

bool Transform3::operator !=(const Transform3& mat) const {
    return !((*this) == mat);
}


Transform3 Transform3::operator -() const {
    return Matrix3::operator -(*this);
}


void Matrix3::translate(float x, T y)
{
    if (x == 0 && y == 0)
        return;

    switch(mLvl) {
        case Identity:
            tx = x;
            ty = y;
            break;
        case Translated:
            tx += x;
            ty += y;
            break;
        case Scaled:
            tx += x * m11;
            ty += y * m22;
            break;
        case Rotated:
            tx += x * m11 + y * m21;
            ty += y * m22 + x * m12;
            break;
        default:break;
    }

    if (mLvl < Translated)
        mLvl = Translated;
}


void Matrix3::scale(float x, T y)
{
    if (x == 1.0 && y == 1.0)
        return;

    switch(mLvl) {
        case Identity:
            m11 = x;
            m22 = y;
            break;
        case Translated:
            m11 = x;
            m22 = y;
            break;
        case Scaled:
            m11 *= x;
            m22 *= y;
            break;
        case Rotated:
            m12 *= x;
            m21 *= y;
            m11 *= x;
            m22 *= y;
            break;
        default :break;
    }

    if (mLvl < Scaled)
        mLvl = Scaled;
}


void Matrix3::rotateInRad(float angle) {
    T sina = M::Sin(angle);
    T cosa = M::Cos(angle);

    switch(mLvl)
    {
        case Identity:
            m11 = cosa;
            m12 = sina;
            m21 = -sina;
            m22 = cosa;
            break;
        case Translated:
            m11 = cosa;
            m12 = sina;
            m21 = -sina;
            m22 = cosa;
            break;
        case Scaled: {
            T s11 = cosa*m11;
            T s12 = sina*m22;
            T s21 = -sina*m11;
            T s22 = cosa*m22;
            m11 = s11;
            m12 = s12;
            m21 = s21;
            m22 = s22;
            break;
        }
        case Rotated:
        {
            T r11 = cosa * m11 + sina * m21;
            T r12 = cosa * m12 + sina * m22;
            T r21 = -sina * m11 + cosa * m21;
            T r22 = -sina * m12 + cosa * m22;
            m11 = r11;
            m12 = r12;
            m21 = r21;
            m22 = r22;
            break;
        }
        default:break;
    }

    if (mLvl < Rotated)
        mLvl = Rotated;
}


Matrix3 Matrix3::translated(float x, T y) const {
    if (x == 0 && y == 0)
        return *this;

    TMatrix3 matrix();

    //TODO refactorize
    matrix.mLvl = mLvl;
    switch(matrix.mLvl) {
        case Identity:
            matrix.tx = x;
            matrix.ty = y;
            break;
        case Translated:
            matrix.tx = tx;
            matrix.ty = ty;
            matrix.tx += x;
            matrix.ty += y;
            break;
        case Scaled:
            matrix.tx = tx;
            matrix.ty = ty;
            matrix.tx += x * m11;
            matrix.ty += y * m22;
            break;
        case Rotated:
            matrix.tx = tx;
            matrix.ty = ty;
            matrix.tx += x*m11 + y*m21;
            matrix.ty += y*m22 + x*m12;
            break;
    }

    if (matrix.mLvl < Translated)
        matrix.mLvl = Translated;

    return matrix;

}


Matrix3 Matrix3::scaled(float x, T y) const
{
    if (x == 1.0 && y == 1.0)
        return *this;

    TMatrix3 matrix = *this;
    switch(matrix.mLvl) {
        case Identity:
            matrix.m11 = x;
            matrix.m22 = y;
            break;
        case Translated:
            matrix.m11 = x;
            matrix.m22 = y;
            break;
        case Scaled:
            matrix.m11 *= x;
            matrix.m22 *= y;
            break;
        case Rotated:
            matrix.m12 *= x;
            matrix.m21 *= y;
            matrix.m11 *= x;
            matrix.m22 *= y;
            break;
        default :break;
    }

    if (matrix.mLvl < Scaled)
        matrix.mLvl = Scaled;

    return matrix;
}


Matrix3 Matrix3::rotated(float angle) const
{
    T sina = M::Sin(angle);
    T cosa = M::Cos(angle);

    TMatrix3 matrix = TMatrix3::Identity;

    matrix.mLvl = mLvl;
    switch(matrix.mLvl)
    {
        case Identity:
            matrix.m11 = cosa;
            matrix.m12 = sina;
            matrix.m21 = -sina;
            matrix.m22 = cosa;
            break;
        case Translated:
            matrix.m11 = cosa;
            matrix.m12 = sina;
            matrix.m21 = -sina;
            matrix.m22 = cosa;
            break;
        case Scaled: {
            matrix.s11 = cosa * m11;
            matrix.s12 = sina * m22;
            matrix.s21 = -sina * m11;
            matrix.s22 = cosa * m22;
            break;
        }
        case Rotated:
        {
            matrix.m11 = cosa * m11 + sina * m21;
            matrix.m12 = cosa * m12 + sina * m22;
            matrix.m21 = -sina * m11 + cosa * m21;
            matrix.m22 = -sina * m12 + cosa * m22;
            break;
        }
    }

    if (matrix.mLvl < Rotated)
        matrix.mLvl = Rotated;

    return matrix;
}

*/
