#include "BAL.h"

template<typename T>
void GetData(FILE* fptr, const char* format, T* value) //value必须是指针
{
    int num_scanned = fscanf(fptr, format, value);
    if(num_scanned != 1)
    {
        cerr << "Get Data Failed!" << endl;
    }
}

double Median(std::vector<double> *data) 
{
    int n = data->size();
    std::vector<double>::iterator mid_point = data->begin() + n / 2;
    std::nth_element(data->begin(), mid_point, data->end());
    return *mid_point;
}

void PerturbPoint3(const double sigma, double *point) 
{
    for (int i = 0; i < 3; ++i)
    {
       point[i] += RandNormal() * sigma;
    }
}

void BALProblem::Cam2AxisCenter(const double* camera, double* axis, double* center) const
{
    axis[0] = camera[0];
    axis[1] = camera[1];
    axis[2] = camera[2];

    double axis_[3] = {0}; // rotation inv
    axis_[0] = -axis[0];
    axis_[1] = -axis[1];
    axis_[2] = -axis[2];

    Rot2Point(axis_, camera + 3, center);

    VectorRef(center, 3) *= -1.0;
}

void BALProblem::AxisCenter2Cam(double* camera, double* axis, const double* center) const
{
    VectorConstRef axis_(axis, 3);
    VectorRef(camera, 3) = axis_;

    Rot2Point(axis, center, camera + 3);
    VectorRef(camera + 3, 3) *= -1.0;
}

BALProblem::BALProblem(const string& path)
{
    FILE* fptr = fopen(path.c_str(), "r");

    if(fptr == NULL)
    {
        cerr << "Open File Failed" << endl;
        return ;
    }

    GetData(fptr, "%d", &num_cam_);
    GetData(fptr, "%d", &num_pts_);
    GetData(fptr, "%d", &num_obs_);

    cout << "Load:" << num_cam_ << "," << 
                    num_obs_ << "," <<
                    num_pts_ << endl;

    index_cam_ = new int[num_obs_];
    index_pts_ = new int[num_obs_];
    obs_ = new double[num_obs_ * 2];
    
    num_param_ = 9 * num_cam_ + 3 * num_pts_;
    param_ = new double[num_param_];

    for(int i = 0; i < num_obs_; i++)
    {
        GetData(fptr, "%d", index_cam_ + i);
        GetData(fptr, "%d", index_pts_ + i);
        GetData(fptr, "%lf", obs_ + 2 * i);
        GetData(fptr, "%lf", obs_ + 2 * i + 1);
    }

    for(int i = 0; i < num_param_; i++)
    {
        GetData(fptr, "%lf", param_ + i);
    }

    fclose(fptr);
}

void BALProblem::normalize()
{
    vector<double> temp(num_pts_);
    double* points = mutable_points();
    Vector3d median;
    for(int i = 0; i < 3; i++)
    {
        for(int j = 0; j < num_pts_; j++)
        {
            temp[j] = points[3 * j + i];
        }
        median(i) = Median(&temp);
    }

    for(int i = 0; i < num_pts_; i++)
    {
        VectorRef point(points + 3 * i, 3);
        temp[i] = (point - median).lpNorm<1>();
    }

    const double median_absolute_deviation = Median(&temp);

    const double scale = 100.0 / median_absolute_deviation;

    for(int i = 0; i < num_pts_; i++)
    {
        VectorRef point(points + 3 * i, 3);
        point = (point - median) * scale;
    }

    double* cameras = mutable_cameras();
    double center[3];
    double axis[3];

    for(int i = 0; i < num_cam_; i++)
    {
        Cam2AxisCenter(cameras + 9 * i, axis, center);
        VectorRef(center, 3) -= median;
        VectorRef(center, 3) *= scale;
        AxisCenter2Cam(cameras + 9 * i, axis, center);

    }
}

void BALProblem::perturb(const double point_sigma, const double rot_sigma, const double tran_sigma)
{
    assert(point_sigma >= 0);
    assert(rot_sigma >= 0);
    assert(tran_sigma >= 0);

    double *points = mutable_points();
    if(point_sigma > 0)
    {
        for(int i = 0; i < num_pts_; i++)
        {
            PerturbPoint3(point_sigma, points + 3 * i);
        }
    }

    double* cameras = mutable_cameras();
    double center[3];
    double axis[3];

    for(int i = 0; i < num_cam_; i++)
    {
        if(rot_sigma > 0)
        {   
            Cam2AxisCenter(cameras + 9 * i, axis, center);
            PerturbPoint3(rot_sigma, axis);
            AxisCenter2Cam(cameras + 9 * i, axis, center);
        }
        if(tran_sigma > 0)
        {
            PerturbPoint3(tran_sigma, cameras + 9 * i + 3);
        }
    }
}

void BALProblem::WriteToPLYFile(const std::string &filename) const 
{
    std::ofstream of(filename.c_str());
    of << "ply"
       << '\n' << "format ascii 1.0"
       << '\n' << "element vertex " << num_cam_ + num_pts_
       << '\n' << "property float x"
       << '\n' << "property float y"
       << '\n' << "property float z"
       << '\n' << "property uchar red"
       << '\n' << "property uchar green"
       << '\n' << "property uchar blue"
       << '\n' << "end_header" << std::endl;

    // Export extrinsic data (i.e. camera centers) as green points.
    double angle_axis[3];
    double center[3];
    for (int i = 0; i < num_cam_; i++)
    {
        const double* camera = cameras() + 9 * i;
        Cam2AxisCenter(camera, angle_axis, center);
        of << center[0] << ' ' << center[1] << ' ' << center[2]
           << " 0 255 0" << '\n';
    }

    // Export the structure (i.e. 3D Points) as white points.
    const double *points = param_ + 9 * num_cam_;
    for (int i = 0; i < num_pts_; ++i) 
    {
        const double *point = points + i * 3;
        for (int j = 0; j < 3; ++j) 
        {
            of << point[j] << ' ';
        }
        of << " 255 255 255\n";
    }
    of.close();
}