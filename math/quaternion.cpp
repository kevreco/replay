#include "quaternion.h"

#include "matrix4.h"
#include "matrix3.h"

const Quaternion Quaternion::Identity(0.0f, 0.0f, 0.0f, 1.0f);

Quaternion::Quaternion() : x(0.0f), y(0.0f), z(0.0f), w(1.0f)
{

}
Quaternion::Quaternion(float x, float y, float z, float w)
{
    this->x = x;
    this->y = y;
    this->z = z;
    this->w = w;
}
Quaternion::Quaternion(double x, double y, double z, double w)
{
    this->x = (float)x;
    this->y = (float)y;
    this->z = (float)z;
    this->w = (float)w;
}

Quaternion::Quaternion(const Vector3 & axis, const float & angleInRadian)
{
    set_value_from_axis(axis,angleInRadian);
}

Quaternion::Quaternion(const float & scalar, const Vector3 & v)
{
    x = v.x;
    y = v.y;
    z = v.z;
    w = scalar;
}

/*
Quaternion::Quaternion(const float& pitchX,const float& rollY,const float& yawZ, bool useDegrees)
{
    // Basically we create 3 Quaternions, one for pitch, one for yaw, one for roll
    // and multiply those together.
    // the calculation below does the same, just shorter

    setValueFromEuler(pitchX, rollY,yawZ,useDegrees);
}
*/

Quaternion::Quaternion(const Quaternion & quat)
{
    x = quat.x;
    y = quat.y;
    z = quat.z;
    w = quat.w;
}

void Quaternion::set_value(float x, float y, float z, float w)
{
    this->x = x;
    this->y = y;
    this->z = z;
    this->w = w;
}

void Quaternion::set_value_from_axis(const Vector3 & axis, float angleInRadian)
{

    //float magn = axis.magnitude();

    float sinAngle;
    angleInRadian *=  0.5; /** DEG_TO_RAD*/
    sinAngle = sin(angleInRadian);

    Vector3 normVector(axis);
    normVector.normalize();

    x = (normVector.x * sinAngle);
    y = (normVector.y * sinAngle);
    z = (normVector.z * sinAngle);
    w = cos(angleInRadian);
}

void Quaternion::set_value_from_euler_rad(const float x, const float y, const float z)
{
    float angle;

    /*{
        angleX *= M::DEG_TO_RAD;
        angleY *= M::DEG_TO_RAD;
        angleZ *= M::DEG_TO_RAD;
    }*/
    angle = x * 0.5;
    const double sr = sin(angle);
    const double cr = cos(angle);

    angle = y * 0.5;
    const double sp = sin(angle);
    const double cp = cos(angle);

    angle = z * 0.5;
    const double sy = sin(angle);
    const double cy = cos(angle);

    const double cpcy = cp * cy;
    const double spcy = sp * cy;
    const double cpsy = cp * sy;
    const double spsy = sp * sy;

    this->x = (sr * cpcy - cr * spsy);
    this->y = (cr * spcy + sr * cpsy);
    this->z = (cr * cpsy - sr * spcy);
    this->w = (cr * cpcy + sr * spsy);

    return normalize();
    /*

    double angleX = x;
    double angleY = y;
    double angleZ = z;

    if (useDegrees)
    {
        angleX *= M::DEG_TO_RAD;
        angleY *= M::DEG_TO_RAD;
        angleZ *= M::DEG_TO_RAD;
    }

    angleX *= 0.50f;
    const double sr = sin(angleX);
    const double cr = cos(angleX);

    angleY *= 0.50f;
    const double sp = sin(angleY);
    const double cp = cos(angleY);

    angleZ *= 0.50f;
    const double sy = sin(angleZ);
    const double cy = cos(angleZ);

    const double cpcy = cp * cy;
    const double spcy = sp * cy;
    const double cpsy = cp * sy;
    const double spsy = sp * sy;

    this->x = (float)(sr * cpcy - cr * spsy);
    this->y = (float)(cr * spcy + sr * cpsy);
    this->z = (float)(cr * cpsy - sr * spcy);
    this->w = (float)(cr * cpcy + sr * spsy);

    normalize();*/
}

void Quaternion::set_value_from_euler(const float x, const float y, const float z, bool useDegrees)
{
    float angle;

    /*{
        angleX *= M::DEG_TO_RAD;
        angleY *= M::DEG_TO_RAD;
        angleZ *= M::DEG_TO_RAD;
    }*/
    angle = x * 0.5 * (useDegrees ? M::DEG_TO_RAD : 1.0);
	const double sr = sin(angle);
	const double cr = cos(angle);

    angle = y * 0.5* (useDegrees ? M::DEG_TO_RAD : 1.0);
	const double sp = sin(angle);
	const double cp = cos(angle);

    angle = z * 0.5* (useDegrees ? M::DEG_TO_RAD : 1.0);
	const double sy = sin(angle);
	const double cy = cos(angle);

	const double cpcy = cp * cy;
	const double spcy = sp * cy;
	const double cpsy = cp * sy;
	const double spsy = sp * sy;

	this->x = (sr * cpcy - cr * spsy);
	this->y = (cr * spcy + sr * cpsy);
	this->z = (cr * cpsy - sr * spcy);
	this->w = (cr * cpcy + sr * spsy);

	return normalize();
	/*

    double angleX = x;
    double angleY = y;
    double angleZ = z;

    if (useDegrees)
    {
        angleX *= M::DEG_TO_RAD;
        angleY *= M::DEG_TO_RAD;
        angleZ *= M::DEG_TO_RAD;
    }

    angleX *= 0.50f;
    const double sr = sin(angleX);
    const double cr = cos(angleX);

    angleY *= 0.50f;
    const double sp = sin(angleY);
    const double cp = cos(angleY);

    angleZ *= 0.50f;
    const double sy = sin(angleZ);
    const double cy = cos(angleZ);

    const double cpcy = cp * cy;
    const double spcy = sp * cy;
    const double cpsy = cp * sy;
    const double spsy = sp * sy;

    this->x = (float)(sr * cpcy - cr * spsy);
    this->y = (float)(cr * spcy + sr * cpsy);
    this->z = (float)(cr * cpsy - sr * spcy);
    this->w = (float)(cr * cpcy + sr * spsy);

    normalize();*/
}

Vector3 Quaternion::axis() const
{



     /*   float angle;
        Vector3 res;
		angle = acosf(w);

		// pre-compute to save time
		float sinf_theta_inv = 1.0/sinf(angle);

		// now the vector
		res.x = this->x*sinf_theta_inv;
		res.y = this->y*sinf_theta_inv;
		res.z = this->z*sinf_theta_inv;

		// multiply by 2
		angle*=2;

		return res;*/



/*
    Vector3 res =
            Vector3(q[0], q[1], q[2]);
     const float sinus = res.magnitude();
     if (sinus > 1E-8)
     {
       res[0] /= sinus;
       res[1] /= sinus;
       res[2] /= sinus;
    }
     return (acos(q[3]) <= M_PI/2.0) ? res : -res;*/


    float v = q[3];

	if(1.0f < v){
		v = 1.0f;
	}
	else if(-1.0f > v){
		v = -1.0f;
	}

	float angle = (2.0f * acosf(v));
	float scale = (float)sin(angle / 2);
	Vector3 result;

	if(scale != 0)
	{
		scale = 1.0f / scale;
		result.x = scale * q[0];
		result.y = scale * q[1];
		result.z = scale * q[2];
		return result;
	}

	return result;//vector3(0.0f, 1.0f, 0.0f);


}

Vector3 Quaternion::getAxis() const {

    float sinTheta = sqrtf(1.0f - w * w);

    assert(M::Abs(sinTheta) > M::EPSILON);

    if( M::Abs(sinTheta) > M::EPSILON ) {
        float inv = 1.0f / sinTheta;
        return Vector3(x * inv, y * inv, z * inv);
    } else {
        return Vector3();
    }
}

float Quaternion::angle() const
{
    //return acos(w) * 2.0f;

    float v = q[3];

	if(1.0f < v){
		v = 1.0f;
	}
	else if(-1.0f > v){
		v = -1.0f;
	}

	return (float)((2.0f * acosf(v)) /** RAD_TO_DEG*/);

}




void Quaternion::get_angle_axis(float& angle, Vector3& axis) const
{
	const float scale = sqrtf(x*x + y*y + z*z);

    if (M::Equals(scale, 0.0f) || w > 1.0f || w < -1.0f)
	{
		angle = 0.0f;
		axis.x = 0.0f;
        axis.y = 0.0f;
        axis.z = -1.0f;
	}
	else
	{
		const float invscale = 1.0f / sqrtf(scale);
		angle = 2.0f * acosf(w);
		axis.x = x * invscale;
		axis.y = y * invscale;
		axis.z = z * invscale;
	}
}

Vector3 Quaternion::to_euler() const
{
    Vector3 v;
    const double sw = w*w;
    const double sq = x*x;
    const double sy = y*y;
    const double sz = z*z;

    // heading = rotation about z-axis
    v[2] = (float) (atan2(2.0 * (x*y +z*w),(sq - sy - sz + sw)));

    // bank = rotation about x-axis
    v[0] = (float) (atan2(2.0 * (y*z +x*w),(-sq - sy + sz + sw)));

    // attitude = rotation about y-axis
    v[1] = asinf( M::Clamped( -1.0f,-2.0f * (x*z - y*w), 1.0f) );

    return v;
}

Vector3 Quaternion::to_euler_in_deg() const
{
    Vector3 v = to_euler();
    v[0] = v[0] * M::RAD_TO_DEG;
    v[1] = v[1] * M::RAD_TO_DEG;
    v[2] = v[2] * M::RAD_TO_DEG;
    return v;
}

void Quaternion::normalize()
{

    float magnitude = this->magnitude();
    if(magnitude > 1.0E-8)
    {
        float invDist = 1.0f / magnitude;
        q[0] *= invDist;
        q[1] *= invDist;
        q[2] *= invDist;
        q[3] *= invDist;
    }
    else
	{
        *this = Quaternion::Identity;
	}
}

float Quaternion::magnitude() const
{
    return sqrtf( x*x
                + y* y
                + z* z
                + w* w);
}


Quaternion & Quaternion::operator=(const Quaternion& quaternionToCopy)
{
    x = quaternionToCopy.x;
    y = quaternionToCopy.y;
    z = quaternionToCopy.z;
    w = quaternionToCopy.w;

    return *this;
}

Quaternion Quaternion::operator +( const Quaternion & quaternion ) const
{
    return Quaternion(x + quaternion.x,
                      y + quaternion.y,
                      z + quaternion.z,
                      w + quaternion.w);
}

Quaternion Quaternion::operator* (const Quaternion & quaternion) const
{
    /*double rx, ry, rz, rw;		// temp result

            rw	= quaternionR.w*w - quaternionR.q*q - quaternionR.y*y - quaternionR.z*z;

            rx	= quaternionR.w*q + quaternionR.q*w + quaternionR.y*z - quaternionR.z*y;
            ry	= quaternionR.w*y + quaternionR.y*w + quaternionR.z*q - quaternionR.q*z;
            rz	= quaternionR.w*z + quaternionR.z*w + quaternionR.q*y - quaternionR.y*q;

            return(Quaternion((float)rx, (float)ry, (float)rz, (float)rw));*/
    return Quaternion(w * quaternion.x + x * quaternion.w + y * quaternion.z - z * quaternion.y,
                      w * quaternion.y + y * quaternion.w + z * quaternion.x - x * quaternion.z,
                      w * quaternion.z + z * quaternion.w + x * quaternion.y - y * quaternion.x,
                      w * quaternion.w - x * quaternion.x - y * quaternion.y - z * quaternion.z);
}

Quaternion Quaternion::operator *(float s) const
{
    return Quaternion(x*s,y*s,z*s,w*s);
}

Quaternion& Quaternion::operator *=(const Quaternion &q)
{
    *this = (*this)*q;
    return *this;
}

Quaternion& Quaternion::operator *=(float s)
{
    x*=s;
    y*=s;
    z*=s;
    w*=s;
    return *this;
}

/*! Returns the Quaternion associated 4x4 OpenGL rotation matrix.

 Use \c glMultMatrixd(q.matrix()) to apply the rotation represented by Quaternion \c q to the
 current OpenGL matrix.

 See also getMatrix(), getRotationMatrix() and inverseMatrix().

 \attention The result is only valid until the next call to matrix(). Use it immediately (as shown
 above) or consider using getMatrix() instead.

 \attention The matrix is given in OpenGL format (row-major order) and is the transpose of the
 actual mathematical European representation. Consider using getRotationMatrix() instead. */
const float * Quaternion::matrix() const
{
    static float m[4][4];
    get_matrix(m);
    return (const float*)(m);
}
/*void Quaternion::getMatrix(Matrix4 & mat) const
{
    float fTX  = 2.0f * qx;
	float fTY  = 2.0f * qy;
	float fTZ  = 2.0f * qz;
	float fTWX = fTX * qw;
	float fTWY = fTY * qw;
	float fTWZ = fTZ * qw;
	float fTXX = fTX * qx;
	float fTXY = fTY * qx;
	float fTXZ = fTZ * qx;
	float fTYY = fTY * qy;
	float fTYZ = fTZ * qy;
	float fTZZ = fTZ * qz;

    mat = Matrix4::IDENTITY;

	mat[0]  = 1.0f - ( fTYY + fTZZ );
	mat[1]  = fTXY - fTWZ;
	mat[2]  = fTXZ + fTWY;

	mat[4]  = fTXY + fTWZ;
	mat[5]  = 1.0f - ( fTXX + fTZZ );
	mat[6]  = fTYZ - fTWX;

	mat[8]  = fTXZ - fTWY;
	mat[9]  = fTYZ + fTWX;
	mat[10] = 1.0f - ( fTXX + fTYY );
}*/
void Quaternion::get_matrix(float m[4][4]) const
{
    const double q00 = 2.0l * x * x;
    const double q11 = 2.0l * y * y;
    const double q22 = 2.0l * z * z;

    const double q01 = 2.0l * x * y;
    const double q02 = 2.0l * x * z;
    const double q03 = 2.0l * x * w;

    const double q12 = 2.0l * y * z;
    const double q13 = 2.0l * y * w;

    const double q23 = 2.0l * z * w;

    m[0][0] = 1.0l - q11 - q22;
    m[1][0] =        q01 - q23;
    m[2][0] =        q02 + q13;

    m[0][1] =        q01 + q23;
    m[1][1] = 1.0l - q22 - q00;
    m[2][1] =        q12 - q03;

    m[0][2] =        q02 - q13;
    m[1][2] =        q12 + q03;
    m[2][2] = 1.0l - q11 - q00;

    m[0][3] = 0.0l;
    m[1][3] = 0.0l;
    m[2][3] = 0.0l;

    m[3][0] = 0.0l;
    m[3][1] = 0.0l;
    m[3][2] = 0.0l;
    m[3][3] = 1.0l;
}

/*! Same as getMatrix(), but with a \c GLdouble[16] parameter. See also getInverseMatrix() and Frame::getMatrix(). */
void Quaternion::get_matrix(float m[16]) const
{
  static float mat[4][4];
  get_matrix(mat);
  int count = 0;
  for (int i=0; i<4; ++i)
    for (int j=0; j<4; ++j)
      m[count++] = mat[i][j];
}

/*! Fills \p m with the 3x3 rotation matrix associated with the Quaternion.

  See also getInverseRotationMatrix().

  \attention \p m uses the European mathematical representation of the rotation matrix. Use matrix()
  and getMatrix() to retrieve the OpenGL transposed version. */
void Quaternion::get_rotation_matrix(float m[3][3]) const
{
  static float mat[4][4];
  get_matrix(mat);
  for (int i=0; i<3; ++i)
    for (int j=0; j<3; ++j)
      // Beware of transposition
      m[i][j] = mat[j][i];
}


Vector3 Quaternion::operator* (const Vector3& v) const {
    // nVidia SDK implementation
    Vector3 uv, uuv;
    Vector3 qvec(x, y, z);
    uv = qvec.cross(v);
    uuv = qvec.cross(uv);
    uv *= (2.0f * w);
    uuv *= 2.0f;

    return v + uv + uuv;
}

/*! Returns the image of \p v by the Quaternion rotation.

See also inverseRotate() and operator*(const Quaternion&, const Vec&). */
Vector3 Quaternion::rotate(const Vector3 & v) const
{
    return (*this * Quaternion(0,v) * conjugate()).vector();

    //qDebug() << "Other World :"  << v.x() << " " <<v.y() << " " << v.z();
    //return (*this * Quaternion(v,0) * conjugate()).vector();

    const double q00 = 2.00 * x * x;
    const double q11 = 2.00 * y * y;
    const double q22 = 2.00 * z * z;

    const double q01 = 2.00 * x * y;
    const double q02 = 2.00 * x * z;
    const double q03 = 2.00 * x * w;

    const double q12 = 2.00 * y * z;
    const double q13 = 2.00 * y * w;

    const double q23 = 2.00 * z * w;

    return Vector3((1.0 - q11 - q22)*v.x + (      q01 - q23)*v.y + (      q02 + q13)*v.z,
                    (      q01 + q23)*v.x + (1.0 - q22 - q00)*v.y + (      q12 - q03)*v.z,
                    (      q02 - q13)*v.x + (      q12 + q03)*v.y + (1.0 - q11 - q00)*v.z);
}

/*! Returns the image of \p v by the Quaternion inverse() rotation.

rotate() performs an inverse transformation. Same as inverse().rotate(v). */
Vector3 Quaternion::inverse_rotate(const Vector3& v) const
{
  return inversed().rotate(v);
}

/*! Sets the Quaternion from the three rotated vectors of an orthogonal basis.

  The three vectors do not have to be normalized but must be orthogonal and direct (X^Y=k*Z, with k>0).

  \code
  Quaternion q;
  q.setFromRotatedBasis(X, Y, Z);
  // Now wotate(Vec(1,0,0)) == X and xnverseRotate(X) == Vec(1,0,0)
  // Same goes for Y and Z with Vec(0,1,0) and Vec(0,0,1).
  \endcode

  See also setFromRotationMatrix() and Quaternion(const Vec&, const Vec&). */
void Quaternion::set_from_rotated_basis(const Vector3& x, const Vector3& y, const Vector3& z)
{
  double m[3][3];
  double normX = x.magnitude();
  double normY = y.magnitude();
  double normZ = z.magnitude();

  for (int i=0; i<3; ++i)
    {
      m[i][0] = x[i] / normX;
      m[i][1] = y[i] / normY;
      m[i][2] = z[i] / normZ;
    }

  set_from_rotation_matrix(m);
}

/*! Set the Quaternion from a (supposedly correct) 3x3 rotation matrix.

  The matrix is expressed in European format: its three \e columns are the images by the rotation of
  the three vectors of an orthogonal basis. Note that OpenGL uses a symmetric representation for its
  matrices.

  setFromRotatedBasis() sets a Quaternion from the three axis of a rotated frame. It actually fills
  the three columns of a matrix with these rotated basis vectors and calls this method. */
void Quaternion::set_from_rotation_matrix(const double m[3][3])
{
  // Compute one plus the trace of the matrix
  const double onePlusTrace = 1.0 + m[0][0] + m[1][1] + m[2][2];

  //Check the diagonal
  if (onePlusTrace > 1E-5)
    {
      // Direct computation
      const double s = sqrt(onePlusTrace) * 2.0;
      q[0] = (m[2][1] - m[1][2]) / s;
      q[1] = (m[0][2] - m[2][0]) / s;
      q[2] = (m[1][0] - m[0][1]) / s;
      q[3] = 0.25 * s;
    }
  else
    {
      // Computation depends on major diagonal term
      if ((m[0][0] > m[1][1])&(m[0][0] > m[2][2]))
        {
          const double s = sqrt(1.0 + m[0][0] - m[1][1] - m[2][2]) * 2.0;
          q[0] = 0.25 * s;
          q[1] = (m[0][1] + m[1][0]) / s;
          q[2] = (m[0][2] + m[2][0]) / s;
          q[3] = (m[1][2] - m[2][1]) / s;
        }
      else
        if (m[1][1] > m[2][2])
          {
            const double s = sqrt(1.0 + m[1][1] - m[0][0] - m[2][2]) * 2.0;
            q[0] = (m[0][1] + m[1][0]) / s;
            q[1] = 0.25 * s;
            q[2] = (m[1][2] + m[2][1]) / s;
            q[3] = (m[0][2] - m[2][0]) / s;
          }
        else
          {
            const double s = sqrt(1.0 + m[2][2] - m[0][0] - m[1][1]) * 2.0;
            q[0] = (m[0][2] + m[2][0]) / s;
            q[1] = (m[1][2] + m[2][1]) / s;
            q[2] = 0.25 * s;
            q[3] = (m[0][1] - m[1][0]) / s;
          }
    }
    normalize();
}
Vector3 Quaternion::x_axis() const
{
    float fTy  = 2.0f*y;
    float fTz  = 2.0f*z;
    float fTwy = fTy*w;
    float fTwz = fTz*w;
    float fTxy = fTy*x;
    float fTxz = fTz*x;
    float fTyy = fTy*y;
    float fTzz = fTz*z;

    return Vector3(1.0f-(fTyy+fTzz), fTxy+fTwz, fTxz-fTwy);
}

Vector3 Quaternion::y_axis() const
{
    float fTx  = 2.0f*x;
    float fTy  = 2.0f*y;
    float fTz  = 2.0f*z;
    float fTwx = fTx*w;
    float fTwz = fTz*w;
    float fTxx = fTx*x;
    float fTxy = fTy*x;
    float fTyz = fTz*y;
    float fTzz = fTz*z;

    return Vector3(fTxy-fTwz, 1.0f-(fTxx+fTzz), fTyz+fTwx);
}

Vector3 Quaternion::z_axis() const
{
    float fTx  = 2.0f*x;
    float fTy  = 2.0f*y;
    float fTz  = 2.0f*z;
    float fTwx = fTx*w;
    float fTwy = fTy*w;
    float fTxx = fTx*x;
    float fTxz = fTz*x;
    float fTyy = fTy*y;
    float fTyz = fTz*y;

    return Vector3(fTxz+fTwy, fTyz-fTwx, 1.0f-(fTxx+fTyy));
}
Vector3 Quaternion::xy_axis() const
{
    Vector3 xy[2];
    get_xy(xy);
    return (xy[0] + xy[1]).normalized();
}

void Quaternion::get_xyz(Vector3 axis[3]) const
{
     //x
    float fTy  = 2.0f*y;
    float fTz  = 2.0f*z;
    float fTwy = fTy*w;
    float fTwz = fTz*w;
    float fTxy = fTy*x;
    float fTxz = fTz*x;
    float fTyy = fTy*y;
    float fTzz = fTz*z;

    axis[0] = Vector3(1.0f-(fTyy+fTzz), fTxy+fTwz, fTxz-fTwy);
    //y
    float fTx  = 2.0f*x;
    //float fTy  = 2.0f*y;
    //float fTz  = 2.0f*z;
    float fTwx = fTx*w;
    //float fTwz = fTz*w;
    float fTxx = fTx*x;
    //float fTxy = fTy*x;
    float fTyz = fTz*y;
    //float fTzz = fTz*z;

    axis[1] = Vector3(fTxy-fTwz, 1.0f-(fTxx+fTzz), fTyz+fTwx);

    axis[2] = Vector3(fTxz+fTwy, fTyz-fTwx, 1.0f-(fTxx+fTyy));
}
void Quaternion::get_xyz2(Vector3 axis[3]) const
{
    axis[0].x = 1 - (2*y*y + 2*z*z);
    axis[0].y = 2*x*y + 2*z*w;
    axis[0].z = 2*x*z - 2*y*w;
    axis[1].x = 2*x*y - 2*z*w;
    axis[1].y = 1 - (2*x*x  + 2*z*z);
    axis[1].z = 2*y*z + 2*x*w;
    axis[2].x = 2*x*z + 2*y*w;
    axis[2].y = 2*y*z - 2*x*w;
    axis[2].z = 1 - (2*x*x  + 2*y*y);
}

void Quaternion::get_xy(Vector3 axis[2]) const
{
    //x
    float fTy  = 2.0f*y;
    float fTz  = 2.0f*z;
    float fTwy = fTy*w;
    float fTwz = fTz*w;
    float fTxy = fTy*x;
    float fTxz = fTz*x;
    float fTyy = fTy*y;
    float fTzz = fTz*z;

    axis[0] = Vector3(1.0f-(fTyy+fTzz), fTxy+fTwz, fTxz-fTwy);
    //y
    float fTx  = 2.0f*x;
    float fTwx = fTx*w;
    float fTxx = fTx*x;
    float fTyz = fTz*y;

    axis[1] = Vector3(fTxy-fTwz, 1.0f-(fTxx+fTzz), fTyz+fTwx);
}
void Quaternion::get_xz(Vector3 axis[2]) const
{
     //x
    float fTy  = 2.0f*y;
    float fTz  = 2.0f*z;
    float fTwy = fTy*w;
    float fTwz = fTz*w;
    float fTxy = fTy*x;
    float fTxz = fTz*x;
    float fTyy = fTy*y;
    float fTzz = fTz*z;

    axis[0] = Vector3(1.0f-(fTyy+fTzz), fTxy+fTwz, fTxz-fTwy);

    //z
    float fTx  = 2.0f*x;
    float fTwx = fTx*w;
    float fTxx = fTx*x;

    float fTyz = fTz*y;

    axis[1] = Vector3(fTxz+fTwy, fTyz-fTwx, 1.0f-(fTxx+fTyy));
}

void Quaternion::get_yz(Vector3 axis[2]) const
{
     //y
    float fTx  = 2.0f*x;
    float fTy  = 2.0f*y;
    float fTz  = 2.0f*z;
    float fTwx = fTx*w;
    float fTwz = fTz*w;
    float fTxx = fTx*x;
    float fTxy = fTy*x;
    float fTyz = fTz*y;
    float fTzz = fTz*z;

    axis[0] = Vector3(fTxy-fTwz, 1.0f-(fTxx+fTzz), fTyz+fTwx);

    //z
    float fTwy = fTy*w;
    float fTxz = fTz*x;
    float fTyy = fTy*y;

    axis[1] = Vector3(fTxz+fTwy, fTyz-fTwx, 1.0f-(fTxx+fTyy));
}

Quaternion Quaternion::mirrored(const Vector3& in)
{
    double ow= 0;
    double ox= ( -x*x +y*y +z*z)*in.x +( -2*x*y)*in.y +( -2*x*z)*in.z;
    double oy= ( -2*x*y)*in.x +( x*x -y*y +z*z)*in.y +( -2*y*z)*in.z;
    double oz= ( -2*x*z)*in.x +( -2*y*z)*in.y +( x*x +y*y -z*z)*in.z;

    return Quaternion(ox, oy, oz, ow);
}

Quaternion Quaternion::FromEulerRadian(const float pitchX,const float rollY,const float yawZ)
{
    Quaternion q;
    q.set_value_from_euler_rad(pitchX, rollY, yawZ);
    return q;
}

Quaternion Quaternion::FromEulerDegree(const float pitchX,const float rollY,const float yawZ)
{
    Quaternion q;
    q.set_value_from_euler_rad(pitchX * M::DEG_TO_RAD, rollY * M::DEG_TO_RAD, yawZ * M::DEG_TO_RAD) ;
    return q;
}

Quaternion Quaternion::FromAxisAngle(const Vector3& axis, float degree) {
    Quaternion q;
    const float halfAngle = 0.5f * degree * M::DEG_TO_RAD;
    const float sin = M::Sin(halfAngle);
    q.w = M::Cos(halfAngle);
    q.x = axis.x * sin;
    q.y = axis.y * sin;
    q.z = axis.z * sin;
    return q;
}

Quaternion Quaternion::FromDirection(const Vector3& dir, const Vector3& up) {

    Vector3 zAxis = dir.normalized();
    Vector3 xAxis = up.cross(zAxis).normalized();
    Vector3 yAxis = zAxis.cross(xAxis);

    Matrix4 m;

    m[0] = xAxis.x;
    m[1] = xAxis.y;
    m[2] = xAxis.z;

    m[4] = yAxis.x;
    m[5] = yAxis.y;
    m[6] = yAxis.z;

    m[8] = zAxis.x;
    m[9] = zAxis.y;
    m[10] = zAxis.z;

    m.transpose();
    return FromMatrix(m);
}

Quaternion Quaternion::FromMatrix(const Matrix4& m) {
    // Compute one plus the trace of the matrix

    const float4x4& data = m.data4x4;

    const double onePlusTrace = 1.0 + data[0][0] + data[1][1] + data[2][2];

    Quaternion q;
    //Check the diagonal
    if (onePlusTrace > 1E-5) {
      // Direct computation
      const double s = sqrt(onePlusTrace) * 2.0;
      q[0] = (data[2][1] - data[1][2]) / s;
      q[1] = (data[0][2] - data[2][0]) / s;
      q[2] = (data[1][0] - data[0][1]) / s;
      q[3] = 0.25 * s;
    } else {
      // Computation depends on major diagonal term
      if ((data[0][0] > data[1][1])&(data[0][0] > data[2][2])) {
          const double s = sqrt(1.0 + data[0][0] - data[1][1] - data[2][2]) * 2.0;
          q[0] = 0.25 * s;
          q[1] = (data[0][1] + data[1][0]) / s;
          q[2] = (data[0][2] + data[2][0]) / s;
          q[3] = (data[1][2] - data[2][1]) / s;
        } else if (data[1][1] > data[2][2]) {
            const double s = sqrt(1.0 + data[1][1] - data[0][0] - data[2][2]) * 2.0;
            q[0] = (data[0][1] + data[1][0]) / s;
            q[1] = 0.25 * s;
            q[2] = (data[1][2] + data[2][1]) / s;
            q[3] = (data[0][2] - data[2][0]) / s;
          } else {
            const double s = sqrt(1.0 + data[2][2] - data[0][0] - data[1][1]) * 2.0;
            q[0] = (data[0][2] + data[2][0]) / s;
            q[1] = (data[1][2] + data[2][1]) / s;
            q[2] = 0.25 * s;
            q[3] = (data[0][1] - data[1][0]) / s;
          }
    }
    q.normalize();
    return q;
}

Quaternion Quaternion::BuildFromLookAt(const Vector3& pos, const Vector3& target) {
    const Matrix4 m = Matrix4::ModelViewFromLookAtLH(pos, target);
    return FromMatrix(m);
}

