#include <chrono>
#include <fstream>
#include <random>
#include "P3P_HD.h"
#include "EPNP/EPNP.h"
#include "cost_function/cost_function.h"
#include <ceres/ceres.h>
using namespace std;

//Update the number of iterations of Ransac function
int updateNumIters(double p, double ep, int niters,int sample_nums) {
    p = max(p, 0.);
    p = min(p, 1.);
    ep = max(ep, 0.);
    ep = min(ep, 1.);

    //置信度
    double num = log(1-p);
    double denom = 1 - pow(ep,sample_nums);

    denom = std::log(denom);
    return denom >= 0 || -num >= niters * (-denom) ? niters : round(num / denom);
}

//generate random numbers
int getRand(int min, int max) {

    std::mt19937_64 generator(static_cast<uint64_t>(std::chrono::system_clock::now().time_since_epoch().count()));
    std::uniform_int_distribution<int> rn(min, max);
    return rn(generator);
}
// load datas
vector<double> loadData(string filePath) {
    ifstream fin(filePath);
    vector<double> points;
    if (!fin.is_open()) {
        cerr <<filePath + " open failed\n";
        exit(EXIT_FAILURE);
    }
    string line;
    getline(fin, line);

    istringstream linestream1(line);

    double number;
    while (linestream1 >> number) {
        points.push_back(number);
    }
    fin.close();
    return points;
}

// This section is adapted from Colmap's pose estimation process.
// We replaced the original P3P solver with our custom P3P solver for better accuracy.
// Reference: Colmap repository (https://github.com/colmap/colmap)
// Our algorithm for Lo_Ransac pipline,where 4 correspondences are sampled in one iteration of Ransac.
// Replace different methods(Mul4, Mul3, Uni-3-1V, Uni-3-3V ) to solve the P3P problems in the Ransac pipline.
int demo_Ransac(hd::P3P_methods methods_id, bool EPNP, bool Ceres) {

    vector<double> worldPoints;
    vector<double> imagePoints;
    vector<double> Intr;
    vector<double> Pose;

    worldPoints=loadData("../Cambridge/extract_data_ShopFacade_1/world.txt");
    imagePoints=loadData("../Cambridge/extract_data_ShopFacade_1/image.txt");
    Intr=loadData("../Cambridge/extract_data_ShopFacade_1/intr.txt");

    double Intrin[3]={Intr[0],Intr[1],Intr[2]};

    const int N = worldPoints.size()/3;

    vector<bool> inliers_flag(N,false);
    int niter = 2000;
    int maxGoodCount = 0;
    int iter=0;

    vector<Eigen::Vector2d> Points2D_in_cam;
    vector<Eigen::Vector2d> Points2D;
    vector<Eigen::Vector3d> Points3D;

    double threshold = 12 / Intrin[0];

    for (size_t k = 0; k < imagePoints.size(); k += 3)
    {
        Eigen::Vector2d vec1((imagePoints[k] - Intrin[1]) / Intrin[0], (imagePoints[k + 1] - Intrin[2]) / Intrin[0]);
        Points2D_in_cam.push_back(vec1);

        Eigen::Vector2d vec2(imagePoints[k], imagePoints[k + 1]);
        Points2D.push_back(vec2);

        Eigen::Vector3d vec3(worldPoints[k], worldPoints[k + 1], worldPoints[k + 2]);
        Points3D.push_back(vec3);
    }

    std::vector<EPNPEstimator::X_t> X_inlier;
    std::vector<EPNPEstimator::Y_t> Y_inlier;
    std::vector<EPNPEstimator::M_t> local_models;
    EPNPEstimator local_estimator;

    Eigen::Matrix<double, 3, 3> R;
    Eigen::Vector3<double> t;

    for (iter = 0; iter < niter; iter++) {
        int r1 = getRand(0, N - 1);
        int r2 = getRand(0, N - 1);
        int r3 = getRand(0, N - 1);
        int r4 = getRand(0, N - 1);

        while (r1 == r2) {
            r2 = getRand(0, N - 1);
        }
        while (r2 == r3 || r1 == r3) {
            r3 = getRand(0, N - 1);
        }
        while (r1 == r4 || r2 == r4 || r3==r4) {
            r4 = getRand(0, N - 1);
        }
        r1=17,r2=4,r3=5,r4=38;
        vector<int> index = {r1, r2, r3, r4};
        double worldPoints_4[12] = {worldPoints[3 * r1], worldPoints[3 * r2], worldPoints[3 * r3], worldPoints[3 * r4],
                                    worldPoints[3 * r1 + 1], worldPoints[3 * r2 + 1], worldPoints[3 * r3 + 1], worldPoints[3 * r4 + 1],
                                    worldPoints[3 * r1 + 2], worldPoints[3 * r2 + 2], worldPoints[3 * r3 + 2], worldPoints[3 * r4 + 2]};

        double imagePoints_4[12] = {imagePoints[3 * r1], imagePoints[3 * r2], imagePoints[3 * r3], imagePoints[3 * r4],
                                    imagePoints[3 * r1 + 1], imagePoints[3 * r2 + 1], imagePoints[3 * r3 + 1], imagePoints[3 * r4 + 1],
                                    imagePoints[3 * r1 + 2], imagePoints[3 * r2 + 2], imagePoints[3 * r3 + 2], imagePoints[3 * r4 + 2]};

        std::vector<Eigen::Matrix<double, 3, 3>> Rs(4);
        std::vector<Eigen::Vector3<double>> Ts(4);

        int num_sols=0;
        switch (methods_id) {
            case hd::P3P_methods::hd_mul4: {
                num_sols=hd::P3P_HD(imagePoints_4,worldPoints_4,Intrin,Rs,Ts,5,hd::P3P_methods::hd_mul4);
                break;
            }
            case hd::P3P_methods::hd_mul3: {
                num_sols=hd::P3P_HD(imagePoints_4,worldPoints_4,Intrin,Rs,Ts,5,hd::P3P_methods::hd_mul3);
                break;
            }
            case hd::P3P_methods::hd_uni3_1v: {
                num_sols=hd::P3P_HD(imagePoints_4,worldPoints_4,Intrin,Rs,Ts,5,hd::P3P_methods::hd_uni3_1v);
                break;
            }
            case hd::P3P_methods::hd_uni3_3v: {
                num_sols=hd::P3P_HD(imagePoints_4,worldPoints_4,Intrin,Rs,Ts,5,hd::P3P_methods::hd_uni3_3v);
                break;
            }
            default:
                break;
        }

        double err = 0;
        for (int ii = 0; ii < num_sols; ii++) {
            int num_inlier = 0;
            vector<bool> temp_flag(N, false);
            for (int j = 0; j < N; j++) {
                Eigen::Vector3<double> x0 = {worldPoints[3*j], worldPoints[3*j+1], worldPoints[3*j+2]};
                Eigen::Vector3<double> x_c = Rs[ii] * x0 + Ts[ii];
                double x_tmp = 1 / x_c[2];
                x_c[0] = x_c[0] * x_tmp;
                x_c[1] = x_c[1] * x_tmp;
                err = (x_c[0] - Points2D_in_cam[j][0]) * (x_c[0] - Points2D_in_cam[j][0]) +
                          (x_c[1] - Points2D_in_cam[j][1]) * (x_c[1] - Points2D_in_cam[j][1]);
                if ( err < threshold*threshold) {
                    num_inlier++;
                    temp_flag[j] = true;
                }
            }
            if (num_inlier > maxGoodCount) {
                maxGoodCount = num_inlier;
                copy(temp_flag.begin(), temp_flag.end(), inliers_flag.begin());

                R = Rs[ii];
                t = Ts[ii];

                //EPNP
                if (num_inlier >= EPNPEstimator::kMinNumSamples && EPNP)
                    {
                        const size_t kMaxNumLocalTrials = 10;

                        for (size_t local_num_trials = 0; local_num_trials < kMaxNumLocalTrials; ++local_num_trials)
                        {
                            X_inlier.clear();
                            Y_inlier.clear();
                            X_inlier.reserve(N);
                            Y_inlier.reserve(N);

                            for (size_t m = 0; m < inliers_flag.size(); ++m)
                            {
                                if (inliers_flag[m])
                                {
                                    X_inlier.push_back(Points2D_in_cam[m]);
                                    Y_inlier.push_back(Points3D[m]);
                                }
                            }
                            local_estimator.Estimate(X_inlier, Y_inlier, &local_models);

                            const size_t prev_best_num_inliers = maxGoodCount;

                            for (const auto &local_model : local_models)
                            {
                                Eigen::Matrix<double, 3, 3> local_Rs;
                                Eigen::Vector<double, 3> local_Ts;
                                local_Rs <<local_model(0, 0), local_model(0, 1), local_model(0, 2),
                                                                      local_model(1, 0), local_model(1, 1), local_model(1, 2),
                                                                      local_model(2, 0), local_model(2, 1), local_model(2, 2);
                                local_Ts <<local_model(0, 3), local_model(1, 3), local_model(2, 3);

                                int local_inlier = 0;
                                vector<bool> local_flag(N, false);
                                for (int k = 0; k < N; k++)
                                {
                                    Eigen::Vector3<double> x0;
                                    x0 << worldPoints[3 * k], worldPoints[3 * k + 1], worldPoints[3 * k + 2];
                                    Eigen::Vector3<double> x_c = local_Rs * x0 + local_Ts;
                                    double x_tmp = 1 / x_c[2];
                                    x_c[0] = x_c[0] * x_tmp;
                                    x_c[1] = x_c[1] * x_tmp;
                                    double err = (x_c[0] - Points2D_in_cam[k][0]) * (x_c[0] - Points2D_in_cam[k][0]) +
                                          (x_c[1] - Points2D_in_cam[k][1]) * (x_c[1] - Points2D_in_cam[k][1]);
                                    if (err < threshold * threshold)
                                    {
                                        local_inlier++;
                                        local_flag[k] = true;
                                    }
                                }

                                if (local_inlier > maxGoodCount)
                                {
                                    maxGoodCount = local_inlier;
                                    copy(local_flag.begin(), local_flag.end(), inliers_flag.begin());
                                    R = local_Rs;
                                    t = local_Ts;
                                }
                            }

                            if (maxGoodCount <= prev_best_num_inliers)
                                break;
                        }
                    }

                if(methods_id==hd::P3P_methods::hd_mul3||methods_id==hd::P3P_methods::hd_mul4)
                    niter = updateNumIters(0.99, (double)maxGoodCount / N, niter,3);
                else
                    niter = updateNumIters(0.99, (double)maxGoodCount / N, niter,4);

            }
        }

    }
    // Ceres----------------------------------------------------------------------------------------------------

    if (Ceres)
    {
        Eigen::Matrix3x4d cam_from_world_matrix{{R(0, 0), R(0, 1), R(0, 2), t(0, 0)},
                                                {R(1, 0), R(1, 1), R(1, 2), t(1, 0)},
                                                {R(2, 0), R(2, 1), R(2, 2), t(2, 0)}};

        Rigid3d *cam_from_world = new Rigid3d(Eigen::Quaterniond(cam_from_world_matrix.leftCols<3>()),
                                              cam_from_world_matrix.col(3));

        const auto loss_function = std::make_unique<ceres::CauchyLoss>(1);
        vector<double> cam_param = vector<double>{Intrin[0], Intrin[1], Intrin[2], 0.0};
        double *camera_params = cam_param.data();
        double *rig_from_world_rotation = cam_from_world->rotation.coeffs().data();
        double *rig_from_world_translation = cam_from_world->translation.data();

        ceres::Problem::Options problem_options;
        problem_options.loss_function_ownership = ceres::DO_NOT_TAKE_OWNERSHIP;
        ceres::Problem problem(problem_options);

        for (int k = 0; k < N; ++k)
        {
            if (!inliers_flag[k])
                continue;
            problem.AddResidualBlock(
                CameraCostFunction<ReprojErrorConstantPoint3DCostFunction>(
                    CameraModelId::kSimpleRadial, Points2D[k], Points3D[k]),
                loss_function.get(),
                rig_from_world_rotation,
                rig_from_world_translation,
                camera_params);
        }

        if (problem.NumResiduals() > 0)
        {
            SetQuaternionManifold(&problem, rig_from_world_rotation);
            // Camera parameterization.
            problem.SetParameterBlockConstant(camera_params);
        }

        ceres::Solver::Options solver_options;
        solver_options.gradient_tolerance = 1.0;
        solver_options.max_num_iterations = 100;
        solver_options.linear_solver_type = ceres::DENSE_QR;
        solver_options.logging_type = ceres::LoggingType::SILENT;

        // The overhead of creating threads is too large.
                solver_options.num_threads = 1;
#if CERES_VERSION_MAJOR < 2
                solver_options.num_linear_solver_threads = 1;
#endif // CERES_VERSION_MAJOR

       ceres::Solver::Summary summary;
       //            auto start_Ceres=std::chrono::high_resolution_clock ::now();
       ceres::Solve(solver_options, &problem, &summary);
       //            auto end_Ceres=std::chrono::high_resolution_clock ::now();
       //            auto duration_Ceres=std::chrono::duration_cast<std::chrono::nanoseconds>(end_Ceres-start_Ceres);
       //            Ceres_time_array.push_back(duration_Ceres.count());

       Eigen::Matrix3x4d transform_matrix = cam_from_world->ToMatrix();
       Eigen::Matrix<double, 3, 3> ceres_Rs;
       Eigen::Vector<double, 3> ceres_Ts;
       ceres_Rs << transform_matrix(0, 0), transform_matrix(0, 1), transform_matrix(0, 2),
                                             transform_matrix(1, 0), transform_matrix(1, 1), transform_matrix(1, 2),
                                             transform_matrix(2, 0), transform_matrix(2, 1), transform_matrix(2, 2);
       ceres_Ts << transform_matrix(0, 3), transform_matrix(1, 3), transform_matrix(2, 3);

       R = ceres_Rs;
       t = ceres_Ts;
    }

    //--------------------------------------------------------------------------------------------------------------

    //print the result
    switch (methods_id) {
        case hd::P3P_methods::hd_mul4: {
            cout<<"Ransac(Mul-4):"<<endl;
            break;
        }
        case hd::P3P_methods::hd_mul3: {
            cout<<"Ransac(Mul-3):"<<endl;
            break;
        }
        case hd::P3P_methods::hd_uni3_1v: {
            cout<<"Ransac(Uni-3-1V):"<<endl;
            break;
        }
        case hd::P3P_methods::hd_uni3_3v: {
            cout<<"Ransac(Uni-3-3V):"<<endl;
            break;
        }
        default:
            break;
    }

    cout << "The num of inlier ratio is :" << double(maxGoodCount)/N << endl;
    cout << "The iteration is :" << iter << endl;
    cout << "The matrix of the P3P problems is: " << endl;

    cout << "R matrix: " << endl;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            cout << R(i * 3 + j) << " ";
        }
        cout << endl;
    }
    cout << "t matrix: " << endl;
    for (int i = 0; i < 3; i++) {
        cout << t(i) << " ";

    }
    cout << endl << endl;
    return 0;
}