#ifndef __COSTFUN_H__
#define __COSTFUN_H_

#include <ceres/ceres.h>
#include <iostream>
#include "Rotation.h"
#include <chrono>

class MyCostFunction
{
    public:
        MyCostFunction(double obs_x, double obs_y) : obs_x_(obs_x), obs_y_(obs_y) 
        {

        }

        template<typename T>
        bool operator()(const T* camera, const T* points, T* residual) const
        {
            T pre[2];
            CamProjectionWithDistortion(camera, points, pre);

            residual[0] = pre[0] - T(obs_x_);
            residual[1] = pre[1] - T(obs_y_);

            return true;
        }

        template<typename T>    
        static inline bool CamProjectionWithDistortion(const T* camera, const T* points, T* pre) // projection
        {
            T p[3];
            Rot2Point(camera, points, p);

            p[0] += camera[3];
            p[1] += camera[4];
            p[2] += camera[5];


            T x_norm = -p[0] / p[2];
            T y_norm = -p[1] / p[2];

            const T& l1 = camera[7];
            const T& l2 = camera[8];

            T r2 = x_norm * x_norm + y_norm * y_norm;
            T distortion = T(1.0) + r2 * (l1 + l2 * r2);

            const T& f = camera[6];

            pre[0] = f * distortion * x_norm;
            pre[1] = f * distortion * y_norm;
            
            return true;
        }

        static ceres::CostFunction* Create(const double obs_x, const double obs_y)
        {
            return (new ceres::AutoDiffCostFunction<MyCostFunction, 2, 9, 3>(new MyCostFunction(obs_x, obs_y)));
        }

    private:
        double obs_x_;
        double obs_y_;
};

#endif