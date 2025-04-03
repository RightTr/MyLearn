#include <opencv2/core.hpp>
#include <iostream>
#include <ceres/ceres.h>
#include <chrono>

using namespace std;
using namespace ceres;

struct COSTFUNCTOR
{
    COSTFUNCTOR(double x, double y) : _x(x), _y(y)
    {

    }

    const double _x, _y;

    template<typename T>
    bool operator()(const T *const abc, T *residual) const
    {
        residual[0] = T(_y) - ceres::exp(abc[0] * T(_x) * T(_x) + abc[1] * T(_x) + abc[2]);
        return true;
    }
};


int main(int argc, char const *argv[])
{
    double a = 1.0, b = 2.0, c = 1.0;
    double a_ = 2.0, b_ = -1.0, c_ = 5.0;
    cv::RNG rng;

    vector<double> arr_y;
    vector<double> arr_x;

    for(int i = 0; i < 100; i++)
    {
        double x = i / 100.0;
        arr_x.push_back(x);
        arr_y.push_back(ceres::exp(a * x * x + b * x + c) + rng.gaussian(1.0 * 1.0));
    }

    double abc[3] = {a_, b_, c_};

    Problem problem;
    for(int i = 0; i < 100; i++)
    {
        CostFunction* costfunction = new AutoDiffCostFunction<COSTFUNCTOR, 1, 3>(new COSTFUNCTOR(arr_x[i], arr_y[i]));
        problem.AddResidualBlock(costfunction, nullptr, abc);
    }

    Solver::Options options;
    options.linear_solver_type = DENSE_NORMAL_CHOLESKY;
    options.minimizer_progress_to_stdout = true;

    Solver::Summary summary;
    chrono::steady_clock::time_point t1 = chrono::steady_clock::now();
    Solve(options, &problem, &summary);
    chrono::steady_clock::time_point t2 = chrono::steady_clock::now();
    chrono::duration<double> duration = chrono::duration_cast<chrono::duration<double>>(t2 - t1);
    cout << "Duration: " << duration.count() << endl;

    cout << summary.BriefReport() << endl;
    cout << abc[0] << "," << abc[1] << "," << abc[2] << endl;

    return 0;
}
