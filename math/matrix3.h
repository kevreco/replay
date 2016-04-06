#ifndef MATRIX3_H
#define MATRIX3_H

#include <math/vector2.h>
#include <math/vector3.h>

// THIS CLASS DOES NOT WORK !!! CAREFULLY


typedef float float9[9];
typedef float float3x3[3][3];
typedef Vector3 vector3x3[3];

struct Matrix3 {

    union {
        float9 data;
        float3x3 data3x3;

        struct {
            float m11, m12, m13,
            m21, m22, m23;
            union {
                struct {float m31, m32;};
                struct {float tx, ty;}; //translation
            };
            float m33;
        };
    };

    //Enum to know the "complexity" of transformation
    //two matrix "only translated" can be translate more easly than two matrix scaled
    //two matrix "scaled" can be scale more easly than two matrix rotated etc.
    //for now, Rotated > Scaled > Translated > None
    //then greater matrix can be retrograded.
    //Ex : a rotated matrix can't be back to scaled or translated matrix.

public:

    //Operator
    Matrix3 operator *(const Matrix3& mat) const;
    //matrix3 operator *(const Matrix3& mat) const;
    Matrix3 operator *(float value) const;
    Matrix3 operator /(float value) const;

    Matrix3& operator *=(float value);
    Matrix3& operator /=(float value);

    bool operator ==(const Matrix3& mat) const;
    bool operator !=(const Matrix3& mat) const;

    Matrix3 operator -() const;

    //TVector3<T> operator *(const TVector3<T>& vec) const;

    inline float operator [](int i) const { return data[i];}
    // inline float* operator [] (size_t row) const { return (float*)data3x3[row];}
    inline float& operator [](int i) { return data[i]; }
    //operator const float* () const { return (const float*) this;}
    inline operator const float* () const { return (const float*) data;}

public:
    static const Matrix3 Identity;

public:
    void toIdentity();
    float determinant() const;

    void setTranslation(float x, float y);
    //void translate(float x, float y);
    //inline void translate(const Vector2& vec) {translate(vec.x, vec.y);}
    Vector2 translation() const;

    void scale(float sx, float sy);
    inline void scale(const Vector2& s) { scale(s.x, s.y);} //overload

    Vector2 scale() const;
    //void rotate(float angle, bool inDegree=true);
    //void rotateInDeg(float angle);
    //void rotateInRad(float angle);

    //matrix3 translated(float x, float y) const;
    //matrix3 translated(const Vector2& vec) const { return translated(vec.x, vec.y);}
    //matrix3 scaled(float sx, float sy) const;
    //matrix3 rotated(float angle) const;

    void invert();
    //matrix3 inverted() const;

    void apply(float px, float py, float& rx, float& ry) const;
    Vector2 apply(float px, float py) const;
    Vector2 apply( const Vector2& vec) const { return apply(vec.x, vec.y);}

    const float* dataPtr() const { return (const float*) this;}

    /*
    friend std::ostream& operator<<(std::ostream& os, const Matrix3& m) {
        os << "\n("
            << m.data[0]<<", "<<m.data[1]<<", "<<m.data[2] << "\n"
            << m.data[3]<<", "<<m.data[4]<<", "<<m.data[5] << "\n"
            << m.data[6]<<", "<<m.data[7]<<", "<<m.data[8] << ")";
        return os;
    }
    */

    inline const char* c_str()
    {
        static char str[512] = {};

        sprintf(str, "\n(%f, %f, %f\n %f, %f, %f,\n %f, %f, %f)",
                data[0], data[1], data[2],
                data[3], data[4], data[5],
                data[6], data[7], data[8]);

        return str;
    }


}; //end class Matrix3


/*
 *
// Transform3 is not commutative !
class Transform3 : public Matrix3 {

    TransformationLevel mLvl;
    //Enum to know the "complexity" of transformation
    //two matrix "only translated" can be translate more easly than two matrix scaled
    //two matrix "scaled" can be scale more easly than two matrix rotated etc.
    //for now, Rotated > Scaled > Translated > None
    //then greater matrix can be retrograded.
    //Ex : a rotated matrix can't be back to scaled or translated matrix.



public:
    enum TransformationLevel {
        Identity = 0,
        Translated,
        Scaled,
        Rotated,
        NumTransformation
    };

    void setLvl(TransformationLevel lvl) {mLvl=lvl;}

    inline Transform3() : Matrix3(),
                        mLvl(Identity) {}

    //Transform3(const Matrix3& mat);
    Transform3(const Vector2& translation);

public:

    //Operator
    Transform3 operator *(const Transform3& mat) const;
    //Transform3 operator *(const Transform3& mat) const;
    Transform3 operator *(float value) const;
    Transform3 operator /(float value) const;

    Transform3& operator *=(const Transform3& mat);
    Transform3& operator *=(float value);
    Transform3& operator /=(float value);

    bool operator ==(const Transform3& mat) const;
    bool operator !=(const Transform3& mat) const;

    Transform3 operator -() const;

    //TVector3<T> operator *(const TVector3<T>& vec) const;

    inline float operator [](int i) const { return data[i];}
    inline float* operator [] (size_t row) const { return (float*)data2[row];}
    inline float operator [](int i) { return data[i]; }
    //operator const float* () const { return (const float*) this;}
    inline operator const float* () const { return (const float*) data;}

public:
    static const Transform3 Zero;
    static const Transform3 Identity;

public:
    void toIdentity();
    inline TransformationLevel lvl() const { return mLvl;}
    inline bool isIdentity() const { return mLvl==Identity;}
    float determinant() const;

    void setTranslation(float x, float y);
    void translate(float x, float y);
    inline void translate(const Vector2& vec) {translate(vec.x, vec.y);}
    Vector2 translation() const;

    void scale(float sx, float sy);
    inline void scale(const Vector2& s) { scale(s.x, s.y);} //overload

    Vector2 scale() const;
    void rotate(float angle, bool inDegree=true);
    void rotateInDeg(float angle);
    void rotateInRad(float angle);

    void invert();

    Transform3 translated(float x, float y) const;
    Transform3 translated(const Vector2& vec) const { translated(vec.x, vec.y);}
    Transform3 scaled(float sx, float sy) const;
    Transform3 rotated(float angle) const;

    Transform3 inverted() const;

    void apply(float px, float py, floatrx, float ry) const;
    Vector2 apply(float px, float py) const;
    inline Vector2 apply( const Vector2& vec) const { apply(vec.x, vec.y);}

    const float* dataPtr() const { return (const float*) this;}

    friend std::ostream& operator<<(std::ostream& os, const Transform3& m)
    {
        os << "\n("
            << m.data[0]<<", "<<m.data[1]<<", "<<m.data[2] << "\n"
            << m.data[3]<<", "<<m.data[4]<<", "<<m.data[5] << "\n"
            << m.data[6]<<", "<<m.data[7]<<", "<<m.data[8] << ")";
        return os;
    }



}; //end class Transform3

*/
#endif // MATRIX3_H
