#include "BAL.h"
#include <g2o/core/base_vertex.h>
#include <g2o/core/base_binary_edge.h>
#include <g2o/core/block_solver.h>
#include <g2o/core/optimization_algorithm_levenberg.h>
#include <g2o/solvers/csparse/linear_solver_csparse.h>
#include <g2o/core/robust_kernel_impl.h>
#include <iostream>

#include "sophus/se3.hpp"

using namespace Sophus;
using namespace std;
using namespace Eigen;

class PoseAndIntrinsics
{
    public:
        PoseAndIntrinsics(){}

        explicit PoseAndIntrinsics(const double* data) //显示调用
        {
            rotation = SO3d::exp(Vector3d(data[0], data[1], data[2]));
            transform = Vector3d(data[3], data[4], data[5]);
            focal = data[6];
            k1 = data[7];
            k2 = data[8];
        }

        void to_data(double* data)
        {
            auto fai = rotation.log();
            for(int i = 0; i < 3; i++) data[i] = fai[i];
            for(int i = 0; i < 3; i++) data[i + 3] = transform[i];
            data[6] = focal;
            data[7] = k1;
            data[8] = k2;
        }
        
        SO3d rotation;
        Vector3d transform = Vector3d::Zero();
        double focal = 0;
        double k1 = 0, k2 = 0;
};

class VertexPoseAndIntrinsics : public g2o::BaseVertex<9, PoseAndIntrinsics> //位姿+内参
{
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW; //内存对齐，加速运算

        VertexPoseAndIntrinsics(){}

        virtual void setToOriginImpl() override //重置优化变量
        {
            _estimate = PoseAndIntrinsics();
        }

        virtual void oplusImpl(const double* update) override //更新优化变量
        {
            _estimate.rotation = SO3d::exp(Vector3d(update[0], update[1], update[2])) * _estimate.rotation;
            _estimate.transform += Vector3d(update[3], update[4], update[5]);
            _estimate.focal += update[6];
            _estimate.k1 += update[7];
            _estimate.k2 += update[8];
        }

        Vector2d project(const Vector3d& point) //投影函数
        {
            Vector3d point_ = _estimate.rotation * point + _estimate.transform;
            point_ = -point_ / point_[2];
            double r2 = sqrt(point_[0] * point_[0] + point_[1] * point_[1]);
            double distortion = 1 + r2 * (_estimate.k1 + _estimate.k2 * r2);
            return Vector2d(_estimate.focal * distortion * point_[0],
                            _estimate.focal * distortion * point_[1]);
        }

        virtual bool read(istream &in){return true;}

        virtual bool write(ostream &out) const {return true;}

};

class VertexPoint : public g2o::BaseVertex<3, Vector3d> //空间点
{
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    
        virtual void setToOriginImpl() override //重置优化变量
        {
            _estimate = Vector3d(0, 0, 0);
        }

        virtual void oplusImpl(const double* update) override //更新优化变量
        {
            _estimate += Vector3d(update[0], update[1], update[2]);
        }

        virtual bool read(istream &in){return true;}

        virtual bool write(ostream &out) const {return true;}
};

class EdgeProjection : public g2o:: BaseBinaryEdge<2, Vector2d, //二元边
                    VertexPoseAndIntrinsics, VertexPoint>
{
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

    virtual void computeError() override //计算误差
    {
        auto v0 = (VertexPoseAndIntrinsics*)_vertices[0];
        auto v1 = (VertexPoint*)_vertices[1];
        auto proj = v0->project(v1->estimate()); //投影
        _error = proj - _measurement; //误差
    }

        virtual bool read(istream &in) {return true;}

        virtual bool write(ostream &out) const {return true;}
};

void SolveBA(BALProblem &bal_problem);

string path = "../../../data/BAL.txt";

int main(int argc, char const *argv[])
{
    BALProblem balproblem(path);
    balproblem.normalize();
    balproblem.perturb(0.5, 0.1, 0.5);
    balproblem.WriteToPLYFile("../../../data/initial_g2o.ply");
    SolveBA(balproblem);
    balproblem.WriteToPLYFile("../../../data/final_g2o.ply");
    return 0;
}

void SolveBA(BALProblem &balproblem)
{
    const double *obs = balproblem.observations();
    double* cameras = balproblem.mutable_cameras();
    double* points = balproblem.mutable_points();

    typedef g2o::BlockSolver<g2o::BlockSolverTraits<9, 3>> BlockSolverType; //优化问题维度
    typedef g2o::LinearSolverCSparse<BlockSolverType::PoseMatrixType> LinearSolverType; //矩阵的数据类型
    auto solver = new g2o::OptimizationAlgorithmLevenberg(
        std::make_unique<BlockSolverType>(std::make_unique<LinearSolverType>())); //求解器类型，LM法
    g2o::SparseOptimizer optimizer;
    optimizer.setAlgorithm(solver);
    optimizer.setVerbose(true);

    vector<VertexPoseAndIntrinsics*> pose_intrinsics_vertex;
    vector<VertexPoint*> points_vertex;
    for(int i = 0; i < balproblem.num_cam(); i++) //设置位姿内参顶点
    {
        VertexPoseAndIntrinsics* v = new VertexPoseAndIntrinsics();
        double* camera = cameras + 9 * i;
        v->setId(i);
        v->setEstimate(PoseAndIntrinsics(camera));
        optimizer.addVertex(v);
        pose_intrinsics_vertex.push_back(v);
    }

    for(int i = 0; i < balproblem.num_pts(); i++) //设置三维点顶点
    {
        VertexPoint* v = new VertexPoint();
        double* point = points + 3 * i;
        v->setId(i + balproblem.num_cam());
        v->setEstimate(Vector3d(point[0], point[1], point[2]));
        v->setMarginalized(true); //边缘化
        optimizer.addVertex(v);
        points_vertex.push_back(v);
    }

    for(int i = 0; i < balproblem.num_obs(); i++) //边连接顶点
    {
        EdgeProjection* e = new EdgeProjection;
        e->setVertex(0, pose_intrinsics_vertex[balproblem.cam_index()[i]]); //设置二元边的起点
        e->setVertex(1, points_vertex[balproblem.pts_index()[i]]); //设置二元边的1终点
        e->setMeasurement(Vector2d(obs[2 * i], obs[2 * i + 1])); //设置测量值
        e->setInformation(Matrix2d::Identity()); //设置信息矩阵
        e->setRobustKernel(new g2o::RobustKernelHuber()); //Huber核
        optimizer.addEdge(e);
    }

    optimizer.initializeOptimization();
    optimizer.optimize(40);

    for(int i = 0; i < balproblem.num_cam(); i++)
    {
        double* camera = cameras + 9 * i;
        auto vertex = pose_intrinsics_vertex[i];
        auto estimate = vertex->estimate(); //返回 VertexPoseAndIntrinsics
        estimate.to_data(camera);
    }    

    for(int i = 0; i < balproblem.num_pts(); i++)
    {
        double* point = points + 3 * i;
        auto vertex = points_vertex[i];
        for(int j = 0; j < 3; j++)
        {
            point[j] = vertex->estimate()[j];
        }
    }    


}