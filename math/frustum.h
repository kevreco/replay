#ifndef FRUSTUM_H
#define FRUSTUM_H

#include "vector3.h"
#include "plane.h"

struct FrustumInfo
{
    Vector3 position;
    Vector3 center;
    Vector3 fCorners[4];

}; // struct FrustumInfo

struct Frustum
{

private:
    enum {
        TOP,
        BOTTOM,
        LEFT,
        RIGHT,
        NEARP,
        FARP
    };

public:
    enum Type { PERSPECTIVE, ORTHOGRAPHIC};

    Frustum();
    Frustum(float aspectRatio, float znear, float zfar, float fov);
    Frustum(const Frustum& f);
    ~Frustum() {}

    Frustum& operator =(const Frustum& f);
    //Calculate with camera param
    void compute_planes(const Vector3& pos, const Vector3& zAxis, const Vector3& yAxis);

    Type type;
    float aspectRatio;
    float znear, zfar;
    float fieldOfView; // in radian
    float top, bottom, left, right;

    Plane planes[6];

    FrustumInfo info;
}; // class Frustum

#endif // FRUSTUM_H
