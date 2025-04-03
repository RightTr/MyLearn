#ifndef __BAL_H__
#define __BAL_H__

#include <cstdio>
#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <string>
#include <fstream>
#include <ostream>

#include "Rotation.h"
#include "Random.h"

using namespace std;
using namespace Eigen;

typedef Map<Vector3d> VectorRef;
typedef Map<const Vector3d> VectorConstRef;

class BALProblem
{
    public:
        explicit BALProblem(const string& path);

        ~BALProblem()
        {
            delete[] index_cam_;
            delete[] index_pts_;
            delete[] param_;
            delete[] obs_;
        }

        void normalize();

        void perturb(const double point_sigma, const double rot_sigma, const double tran_sigma);

        void Cam2AxisCenter(const double* camera, double* axis, double* center) const;

        void AxisCenter2Cam(double* camera, double* axis, const double* center) const;

        void WriteToPLYFile(const std::string &filename) const;

        double* mutable_points() const
        {
            return param_ + 9 * num_cam_;
        }

        double* mutable_cameras() const
        {
            return param_;
        }

        const double* observations() const
        {
            return obs_;
        }

        const int* cam_index() const
        {
            return index_cam_;
        }

        const int* pts_index() const
        {
            return index_pts_;
        }

        const int num_obs() const
        {
            return num_obs_;
        }

        const double* cameras() const
        {
            return param_;
        }

    private:
        int num_cam_;
        int num_pts_;
        int num_obs_;
        int num_param_;

        int* index_cam_;
        int* index_pts_;
        double* param_;
        double* obs_;

};

#endif