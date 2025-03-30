#include <iostream>
#include <chrono>
#include <opencv2/core.hpp>
#include <g2o/core/g2o_core_api.h>
#include <g2o/core/base_vertex.h>
#include <g2o/core/base_unary_edge.h>
#include <g2o/core/block_solver.h>
#include <g2o/core/optimization_algorithm_levenberg.h>
#include <g2o/core/optimization_algorithm_gauss_newton.h>
#include <g2o/core/optimization_algorithm_dogleg.h>
#include <g2o/solvers/dense/linear_solver_dense.h>
#include <cmath>
#include <Eigen/Core>

using namespace std;
using namespace g2o;
using namespace Eigen;

class CurveFittingVertex : public BaseVertex<3, Vector3d>
{
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW //内存对齐

        virtual void setToOriginImpl() override
        {
            _estimate << 0, 0, 0;
        }

        virtual void oplusImpl(const double* update) override
        {
            _estimate += Vector3d(update);
        }

        virtual bool read(istream &in)
        {
            return true;
        }

        virtual bool write(ostream &out) const
        {
            return true;
        }
};

class CurveFittingEdge : public BaseUnaryEdge<1, double, CurveFittingVertex>
{
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    
        double _x;

        CurveFittingEdge(double x) : BaseUnaryEdge() , _x(x)
        {

        }

        virtual void computeError() override
        {
            const CurveFittingVertex *v = static_cast<const CurveFittingVertex *> (_vertices[0]);
            const Eigen::Vector3d abc = v->estimate();
            _error[0] = _measurement - std::exp(abc[0] * _x * _x + abc[1] * _x + abc[2]);
        }

        virtual void linearizeOplus() override
        {
            const CurveFittingVertex* v = static_cast<const CurveFittingVertex*>(_vertices[0]);
            const Vector3d abc = v->estimate();
            double y = exp(abc[0] * _x * _x + abc[1] * _x + abc[2]);
            _jacobianOplusXi[0] = -_x * _x * y;
            _jacobianOplusXi[1] = -_x * y;
            _jacobianOplusXi[2] = -y;
        }

        virtual bool read(istream &in)
        {
            return true;
        }
    
        virtual bool write(ostream &out) const
        {
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
        arr_y.push_back(std::exp(a * x * x + b * x + c) + rng.gaussian(1.0 * 1.0));
    }

    typedef BlockSolver<BlockSolverTraits<3, 1>> BlockSolverType;
    typedef LinearSolverDense<BlockSolverType::PoseMatrixType> LinearSolverType;

    auto solver = new OptimizationAlgorithmGaussNewton(std::make_unique<BlockSolverType>
                    (std::make_unique<LinearSolverType>()));

    SparseOptimizer opt;
    opt.setAlgorithm(solver);
    opt.setVerbose(true);

    CurveFittingVertex *v = new CurveFittingVertex();
    v->setEstimate(Eigen::Vector3d(a_, b_, c_));
    v->setId(0);
    opt.addVertex(v);

    for(int i = 0; i < 100; i++)
    {
        CurveFittingEdge* edge = new CurveFittingEdge(arr_x[i]);
        edge->setId(i);
        edge->setVertex(0, v);
        edge->setMeasurement(arr_y[i]);
        edge->setInformation(Matrix<double, 1, 1>::Identity() * 1);
        opt.addEdge(edge);
    }

    opt.initializeOptimization();
    opt.optimize(10);

    Eigen::Vector3d abc_ = v->estimate();
    cout << "estimated: " << abc_.transpose() << endl;

    return 0;
}
