#include "frustum.h"

Frustum::Frustum() :
    type(PERSPECTIVE),
    aspectRatio(16.0f/9.0f),
    znear(0.1f), //0.067
    zfar(300.0f),
    fieldOfView(45.0f)
{
    if (type == PERSPECTIVE)
    {
        top = tanf((fieldOfView * 0.5) * M::DEG_TO_RAD) * znear;
        bottom  = -top;
        left = bottom * aspectRatio;
        right = top * aspectRatio;
    }else
    {
        //TODO: implement ORTHOGRAPHIC

    }

}

Frustum::Frustum(float ratio, float zn, float zf, float fov) :
    type(PERSPECTIVE),
    aspectRatio(ratio),
    znear(zn),
    zfar(zf),
    fieldOfView(fov)
{
    if (type == PERSPECTIVE)
    {
        top = tanf((fieldOfView * 0.5f) * M::DEG_TO_RAD) * znear;
        bottom  = -top;
        left = bottom * aspectRatio;
        right = top * aspectRatio;
    }else
    {
        assert("TODO: implement ORTHOGRAPHIC");
    }
}

Frustum::Frustum(const Frustum& f) :
    type(f.type),
    aspectRatio(f.aspectRatio),
    znear(f.znear),
    zfar(f.zfar),
    fieldOfView(f.fieldOfView), // in radians
    top(f.top),
    bottom(f.bottom),
    left(f.left),
    right(f.right)
{
    planes[TOP] = f.planes[TOP];
    planes[BOTTOM] = f.planes[BOTTOM];
    planes[LEFT] = f.planes[LEFT];
    planes[RIGHT] = f.planes[RIGHT];
    planes[NEARP] = f.planes[NEARP];
    planes[FARP] = f.planes[FARP];

    info = f.info;
}

Frustum& Frustum::operator =(const Frustum& f)
{
    type = (f.type);
    aspectRatio = (f.aspectRatio);
    znear = (f.znear);
    zfar = (f.zfar);
    fieldOfView = (f.fieldOfView); // in radian
    top = (f.top);
    bottom = (f.bottom);
    left = (f.left);
    right = (f.right);

    planes[TOP] = f.planes[TOP];
    planes[BOTTOM] = f.planes[BOTTOM];
    planes[LEFT] = f.planes[LEFT];
    planes[RIGHT] = f.planes[RIGHT];
    planes[NEARP] = f.planes[NEARP];
    planes[FARP] = f.planes[FARP];

    info = f.info;

    return *this;

}

void Frustum::compute_planes(const Vector3& pos, const Vector3& zAxis, const Vector3& yAxis)
{
    Vector3 nc,fc, //center of planes
            X,Y,Z, //axis
            ntl,ntr,nbl,nbr, // 4 near corners
            ftl,ftr,fbl,fbr; // 4 far corners

    X = yAxis.cross(zAxis).normalized();
    Z = zAxis.normalized();
    Y = yAxis.normalized();

    // compute the centers of the near and far planes
    nc = pos + (Z * znear);
    fc = pos + (Z * zfar);
    info.position = pos;

    // update center of frustrum
    info.center = (nc+fc)*0.5f;

    double tan = tanf(fieldOfView*M::DEG_TO_RAD*0.5 ) ;
    // compute near plane dimensions
    float nHeight = znear*tan;
    float nWidth = nHeight*aspectRatio;
    // compute far plane dimensions
    float fHeight = zfar*tan;
    float fWidth = fHeight*aspectRatio;

    // compute the 4 corners of the frustum on the near plane
    ntl = nc + Y * nHeight - X * nWidth;
    ntr = nc + Y * nHeight + X * nWidth;
    nbl = nc - Y * nHeight - X * nWidth;
    nbr = nc - Y * nHeight + X * nWidth;

    // compute the 4 corners of the frustum on the far plane
    info.fCorners[0] = ftl = fc + Y * fHeight - X * fWidth;
    info.fCorners[1] = ftr = fc + Y * fHeight + X * fWidth;
    info.fCorners[2] = fbl = fc - Y * fHeight - X * fWidth;
    info.fCorners[3] = fbr = fc - Y * fHeight + X * fWidth;

    // compute the six planes
    // the function set3Points assumes that the points
    // are given in counter clockwise order
    planes[TOP] = Plane(ntr,ntl,ftl);
    planes[BOTTOM] = Plane(nbl,nbr,fbr);
    planes[LEFT] = Plane(ntl,nbl,fbl);
    planes[RIGHT] = Plane(nbr,ntr,ftr);
    planes[NEARP] = Plane(ntl,ntr,nbr);
    planes[FARP] = Plane(ftr,ftl,fbl);

}
