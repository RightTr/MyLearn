#include <iostream>
#include <chrono>
#include <opencv2/core.hpp>
#include <Eigen/Core>
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

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
        arr_y.push_back(exp(a * x * x + b * x + c) + rng.gaussian(1.0 * 1.0));
    }

    double last_cost = 0;

    chrono::steady_clock::time_point t1 = chrono::steady_clock::now();
    for(int i = 0; i < 100; i++)
    {
        Matrix3d H = Matrix3d::Zero();
        Vector3d g = Vector3d::Zero();
        double cost = 0;

        for(int j = 0; j < 100; j++)
        {
            double x = arr_x[j];
            double y = arr_y[j];

            double e = y - exp(a_ * x * x + b_ * x + c_);
            cost += e * e;
            Vector3d J;
    
            J[0] = -x * x * exp(a_ * x * x + b_ * x + c_);
            J[1] = -x * exp(a_ * x * x + b_ * x + c_);
            J[2] = -exp(a_ * x * x + b_ * x + c_);

            H += J * J.transpose();
            g += -e * J ;
        }

        Vector3d dp = H.ldlt().solve(g);

        cout << "cost: " << cost << endl;

        if(cost >= last_cost && i > 0)
        {
            break;
        }

        a_ += dp[0];
        b_ += dp[1];
        c_ += dp[2];
        last_cost = cost;

        cout << "iteration:" <<  i << ",dp:" << dp[0] << "," << dp[1] << "," << dp[2] << endl; 
        
        if(dp.norm() < 1e-6)
        {
            break;
        }
    }
    chrono::steady_clock::time_point t2 = chrono::steady_clock::now();
    chrono::duration<double> duration = chrono::duration_cast<chrono::duration<double>>(t2 - t1);

    cout << a_ << "," << b_ << "," << c_ << endl;
    cout << "Duration: " << duration.count() << endl;

    return 0;
}
