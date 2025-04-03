#include "CostFun.h"
#include "BAL.h"

using namespace ceres;

void SolveBA(BALProblem& balproblem);

string path = "../../data/BAL.txt";

int main(int argc, char const *argv[])
{
    BALProblem balproblem(path);
    balproblem.normalize();
    balproblem.perturb(0.5, 0.1, 0.5);
    balproblem.WriteToPLYFile("initial.ply");
    SolveBA(balproblem);
    balproblem.WriteToPLYFile("final.ply");

    return 0;
}

void SolveBA(BALProblem& balproblem)
{
    const double *obs = balproblem.observations();
    double* cameras = balproblem.mutable_cameras();
    double* points = balproblem.mutable_points();

    Problem problem;
    for(int i = 0; i < balproblem.num_obs(); i++)
    {
        CostFunction* costfunction = MyCostFunction::Create(obs[2 * i], obs[2 * i + 1]);

        LossFunction* lossfunction = new HuberLoss(1.0);
        
        double* camera = cameras + balproblem.cam_index()[i] * 9;
        double* point = points + balproblem.pts_index()[i] * 3;

        problem.AddResidualBlock(costfunction, lossfunction, camera, point);
    }
    
    Solver::Options options;
    options.linear_solver_type = ceres::LinearSolverType::SPARSE_SCHUR;
    options.minimizer_progress_to_stdout = true;
    options.logging_type = ceres::PER_MINIMIZER_ITERATION;

    Solver::Summary summary;

    ceres::Solve(options, &problem, &summary);

    cout << summary.FullReport() << "\n";
}
