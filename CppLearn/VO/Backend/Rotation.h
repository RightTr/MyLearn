#ifndef __ROTATION_H__
#define __ROTATION_H__

#include <algorithm>
#include <cmath>
#include <limits>


template<typename T>
inline T dot(const T x[3], const T y[3])
{
    return (x[0] * y[0] + x[1] * y[1] + x[2] * y[2]);
}

template<typename T>
inline void cross(const T x[3], const T y[3], T res[3])
{
    res[0] = x[1] * y[2] - x[2] * y[1];
    res[1] = x[2] * y[0] - x[0] * y[2];
    res[2] = x[0] * y[1] - x[1] * y[0];
}

template<typename T>
inline void Rot2Point(const T axis[3], const T pt[3], T res[3]) // Rodrigues' rotation formula
{
    const T thetasquare = dot(axis, axis);
    if (thetasquare > T(std::numeric_limits<double>::epsilon()))
    {
        const T theta = sqrt(thetasquare);
        const T costheta = cos(theta);
        const T sintheta = sin(theta);

        const T theta_inv = 1.0 / theta;

        const T axis_[3] = {
            axis[0] * theta_inv, 
            axis[1] * theta_inv, 
            axis[2] * theta_inv, 
        };

        const T temp = (T(1.0) - costheta) * dot(axis_, pt);

        T w[3];
        cross(axis_, pt, w);

        res[0] = pt[0] * costheta + temp * axis_[0] + sintheta * w[0];
        res[1] = pt[1] * costheta + temp * axis_[1] + sintheta * w[1];
        res[2] = pt[2] * costheta + temp * axis_[2] + sintheta * w[2];
        
    }
    else //when theta is near 0
    {
        T w[3];
        cross(axis, pt, w);

        res[0] = pt[0] + w[0];
        res[1] = pt[1] + w[1];
        res[2] = pt[2] + w[2];

    }

}


#endif