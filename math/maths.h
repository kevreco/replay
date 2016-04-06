#ifndef MATHS_H
#define MATHS_H


#include <cmath>

namespace M {

const double PI= 3.1415926535897932384626433832795028841971693993751;
const double DEG_TO_RAD = (PI/180.0);
const double RAD_TO_DEG = (180.0/PI);

 const double EPSILON = (1.0e-08);

//Returns if a equals b, with rounding errors
inline bool Equals(const float a, const float b, const float tolerance = 0.000001f)
{
    return (a + tolerance >= b) && (a - tolerance <= b);
}

template <class T>
inline const T & Min(const T & a,const T & b)
{
    return (a < b) ? a : b;
}


template <class T>
inline const T & Max(const T& a,const T & b)
{
    return (a > b) ? a : b;
}

// Returns absolute value
template <class T>
inline T Abs(const T & a)
{
    return (a >= 0) ? a :-a;
}

//Return type for simple or double precision
template<class T>
struct PrecisionTrait { typedef float type; };
template<>
struct PrecisionTrait<double> { typedef double type; };

// Sin
template<class T>
inline typename PrecisionTrait<T>::type Sin(const T& a) //use simple precision by default
{
    return sinf(a);
}

template<>
inline PrecisionTrait<double>::type Sin(const double& a) //use double precision
{
    return sin(a);
}

// Cos

template<class T>
inline typename PrecisionTrait<T>::type Cos(const T& a) //use simple precision by default
{
    return cosf(a);
}

template<>
inline PrecisionTrait<double>::type Cos(const double& a) //use double precision
{
    return cos(a);
}

// Sqrt

template<class T>
inline typename PrecisionTrait<T>::type Sqrt(const T& a) //use simple precision by default
{
    return sqrtf(a);
}

template<>
inline PrecisionTrait<double>::type Sqrt(const double& a) //use double precision
{
    return sqrt(a);
}

template <typename T>
inline bool IsFinite(T a)
{
    return std::isfinite(a);
}

template <class T>
inline const T & Clamped(const T & low,const T & value, const T & high)
{
    if (value < low)
        return low;
    else if(value > high)
        return high;
    else
        return value;
}

template <class T>
inline void Clamp(const T & low,T & value, const T & high)
{
    if (value < low)
        value = low;
    else if(value > high)
        value = high;
}

template <class T>
inline int RoundInt(const T & d)
{
    return d >= 0.0 ? int(d + 0.5) : int(d - int(d-1) + 0.5) + int(d-1);
}

}
/*
//#include <Utils/Logger.h>

#include <limits>
#include <cmath>
namespace Math
{
    /// CONSTANTS

    const float PI= 3.14159265358979323846f;
    const double PId= 3.1415926535897932384626433832795028841971693993751;

    const float TWO_PI = 3.14159265358979323846f*2.0f;
    const double TWO_PId= 3.1415926535897932384626433832795028841971693993751*2.0;


    const float DEG_TO_RAD = (PI/180.0f);
    const double DEG_TO_RADd = (PId/180.0);

    const float RAD_TO_DEG = (180.0f/PI);
    const double RAD_TO_DEGd = (180.0/PId);

    const float Tolerance = 1e-06f;
    const double dTolerance = 1e-08;

    const float EPSILON = (1.0e-08f);
    const double dEPSILON = (1.0e-08);

    /// COMMON STUFF

    template <class T>
    inline double DegToRad(T a)
    {
        return a*DEG_TO_RADd;
    }
    template <class T>
    inline double RadToDeg(T a)
    {
        return a*RAD_TO_DEGd;
    }

    template <class T>
    inline const T & Min(const T & a,const T & b)
    {
        return (a < b) ? a : b;
    }


    template <class T>
    inline const T & Max(const T& a,const T & b)
    {
        return (a > b) ? a : b;
    }

    //Return type for simple or double precision
    template<class T>
    struct PrecisionTrait { typedef float type; };
    template<>
    struct PrecisionTrait<double> { typedef double type; };

    //Sin
    template<class T>
    inline typename PrecisionTrait<T>::type Sin(const T& a) //use simple precision by default
    {
        return sinf(a);
    }

    template<>
    inline PrecisionTrait<double>::type Sin(const double& a) //use double precision
    {
        return sin(a);
    }

    //Cos

    template<class T>
    inline typename PrecisionTrait<T>::type Cos(const T& a) //use simple precision by default
    {
        return cosf(a);
    }

    template<>
    inline PrecisionTrait<double>::type Cos(const double& a) //use double precision
    {
        return cos(a);
    }

    //Sqrt

    template<class T>
    inline typename PrecisionTrait<T>::type Sqrt(const T& a) //use simple precision by default
    {
        return sqrtf(a);
    }

    template<>
    inline PrecisionTrait<double>::type Sqrt(const double& a) //use double precision
    {
        return sqrt(a);
    }

    template <class T>
    inline void Threshold(T & value, const T & min,const T & max)
    {
        if (value < min) value = min;
        else if (value > max) value = max;
    }

    template <class T>
    inline void Swap(T & a, T & b)
    {
        T c = a;
        a = b;
        b = c;
    }

    //Returns the absolute value
    template <class T>
    inline T Abs(const T & a)
    {
        return (a>=0)?a:-a;
    }


    template <class T>
    inline int Round(const T & d)
    {
        return d >= 0.0 ? int(d + 0.5) : int(d - int(d-1) + 0.5) + int(d-1);
    }


    //Returns if a equals b, with rounding errors
    inline bool Equals(const float a, const float b, const float tolerance = 0.000001f)
    {
        return (a + tolerance >= b) && (a - tolerance <= b);
    }



    inline void LimitRange(const float & low,float &  value,const float &  high)
    {
        if (value < low)
            value = low;
        else if(value > high)
            value = high;
    }

    template <class T>
    inline const T & Clamped(const T & low,const T & value, const T & high)
    {
        if (value < low)
            return low;
        else if(value > high)
            return high;
        else
            return value;
    }

    template <class T>
    inline void Clamp(const T & low,T & value, const T & high)
    {
        if (value < low)
            value = low;
        else if(value > high)
            value = high;
    }

    // if a = 68 and step = 10; return 7
    // if a = 63 and step = 10; return 6
    inline int RoundStep(int a, int step)
    {
        int half = step/2;
        return (a >= 0.0 ? (a + half) / step : (a - half) / step );
    }

    // T1 should be a float or double, T2 should be a short, int, long etc.
    template <typename T1, typename T2>
    inline T2 RoundInt(T1 a)
    {
        return (a >= 0.0 ? T2(a + 0.5) : T2(a - T2(a-1) + 0.5) + T2(a-1));
    }

    template <typename T>
    inline bool IsFinite(T a)
    {
        return std::isfinite(a);
    }

    template <typename T>
    inline T Average(const T& a, const T& b)
    {
        return (a+b)/static_cast<T>(2);
    }

    template <typename T>
    inline T Average(const T& a, const T& b, const T& c)
    {
        return (a+b+c)/static_cast<T>(3);
    }

    template <typename T>
    inline T Average(const T& a, const T& b, const T& c, const T& d)
    {
        return (a+b+c+d)/static_cast<T>(4);
    }

    inline int NextPowerOf2(int n)
    {
        --n;
        n |= n >> 1;
        n |= n >> 2;
        n |= n >> 4;
        n |= n >> 8;
        n |= n >> 16;
        ++n;
        return n;
    }


    inline int NextMultipleOf4(int n)
    {
        return (n + 3) & ~0x03;
    }

    //PI over 4
    const float PIo4 = M::PI/4.0f;

    //found on http://en.wikipedia.org/wiki/File:Circle_and_cubic_bezier.svg
     const float Kappa = (4.0f/3.0f)*((1.0f-M::Cos(PIo4))/(M::Sin(PIo4)))*1.0f;

    template<int min_range_src, int max_range_src, int min_range_dst, int max_range_dst>
    inline float IntervalTranspose(float v)
    {

         return ( (v - min_range_src) / (max_range_src - min_range_src) ) * (max_range_dst - min_range_dst);
    }

    inline float IntervalTranspose2(float min_range_src, float max_range_src, float min_range_dst, float max_range_dst, float v)
    {
         return ( (v - min_range_src) / (max_range_src - min_range_src) ) * (max_range_dst - min_range_dst);
    }
}


*/
#endif // MATHS_H
