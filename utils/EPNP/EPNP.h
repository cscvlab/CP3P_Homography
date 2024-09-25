#ifndef P3P_REAL_WORLD_EPNP_H
#define P3P_REAL_WORLD_EPNP_H

#include <array>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Geometry>

namespace Eigen {

    typedef Eigen::Matrix<float, 3, 4> Matrix3x4f;
    typedef Eigen::Matrix<double, 3, 4> Matrix3x4d;
    typedef Eigen::Matrix<double, 6, 6> Matrix6d;
    typedef Eigen::Matrix<uint8_t, 3, 1> Vector3ub;
    typedef Eigen::Matrix<uint8_t, 4, 1> Vector4ub;
    typedef Eigen::Matrix<double, 6, 1> Vector6d;

}

void ComputeSquaredReprojectionError(
        const std::vector<Eigen::Vector2d>& points2D,
        const std::vector<Eigen::Vector3d>& points3D,
        const Eigen::Matrix3x4d& cam_from_world,
        std::vector<double>* residuals);

class EPNPEstimator {
public:
    // The 2D image feature observations.
    typedef Eigen::Vector2d X_t;
    // The observed 3D features in the world frame.
    typedef Eigen::Vector3d Y_t;
    // The transformation from the world to the camera frame.
    typedef Eigen::Matrix3x4d M_t;

    // The minimum number of samples needed to estimate a model.
    static const int kMinNumSamples = 4;

    // Estimate the most probable solution of the P3P problem from a set of
    // three 2D-3D point correspondences.
    //
    // @param points2D   Normalized 2D image points as 3x2 matrix.
    // @param points3D   3D world points as 3x3 matrix.
    //
    // @return           Most probable pose as length-1 vector of a 3x4 matrix.
    static void Estimate(const std::vector<X_t>& points2D,
                         const std::vector<Y_t>& points3D,
                         std::vector<M_t>* models);

    // Calculate the squared reprojection error given a set of 2D-3D point
    // correspondences and a projection matrix.
    //
    // @param points2D     Normalized 2D image points as Nx2 matrix.
    // @param points3D     3D world points as Nx3 matrix.
    // @param proj_matrix  3x4 projection matrix.
    // @param residuals    Output vector of residuals.
    static void Residuals(const std::vector<X_t>& points2D,
                          const std::vector<Y_t>& points3D,
                          const M_t& proj_matrix,
                          std::vector<double>* residuals);

private:
    bool ComputePose(const std::vector<Eigen::Vector2d>& points2D,
                     const std::vector<Eigen::Vector3d>& points3D,
                     Eigen::Matrix3x4d* proj_matrix);

    void ChooseControlPoints();
    bool ComputeBarycentricCoordinates();

    Eigen::Matrix<double, Eigen::Dynamic, 12> ComputeM();
    Eigen::Matrix<double, 6, 10> ComputeL6x10(
            const Eigen::Matrix<double, 12, 12>& Ut);
    Eigen::Matrix<double, 6, 1> ComputeRho();

    void FindBetasApprox1(const Eigen::Matrix<double, 6, 10>& L_6x10,
                          const Eigen::Matrix<double, 6, 1>& rho,
                          Eigen::Vector4d* betas);
    void FindBetasApprox2(const Eigen::Matrix<double, 6, 10>& L_6x10,
                          const Eigen::Matrix<double, 6, 1>& rho,
                          Eigen::Vector4d* betas);
    void FindBetasApprox3(const Eigen::Matrix<double, 6, 10>& L_6x10,
                          const Eigen::Matrix<double, 6, 1>& rho,
                          Eigen::Vector4d* betas);

    void RunGaussNewton(const Eigen::Matrix<double, 6, 10>& L_6x10,
                        const Eigen::Matrix<double, 6, 1>& rho,
                        Eigen::Vector4d* betas);

    double ComputeRT(const Eigen::Matrix<double, 12, 12>& Ut,
                     const Eigen::Vector4d& betas,
                     Eigen::Matrix3d* R,
                     Eigen::Vector3d* t);

    void ComputeCcs(const Eigen::Vector4d& betas,
                    const Eigen::Matrix<double, 12, 12>& Ut);
    void ComputePcs();

    void SolveForSign();

    void EstimateRT(Eigen::Matrix3d* R, Eigen::Vector3d* t);

    double ComputeTotalReprojectionError(const Eigen::Matrix3d& R,
                                         const Eigen::Vector3d& t);

    const std::vector<Eigen::Vector2d>* points2D_ = nullptr;
    const std::vector<Eigen::Vector3d>* points3D_ = nullptr;
    std::vector<Eigen::Vector3d> pcs_;
    std::vector<Eigen::Vector4d> alphas_;
    std::array<Eigen::Vector3d, 4> cws_;
    std::array<Eigen::Vector3d, 4> ccs_;
};

struct Rigid3d {
public:
    Eigen::Quaterniond rotation = Eigen::Quaterniond::Identity();
    Eigen::Vector3d translation = Eigen::Vector3d::Zero();

    Rigid3d() = default;
    Rigid3d(const Eigen::Quaterniond& rotation,
            const Eigen::Vector3d& translation)
            : rotation(rotation), translation(translation) {}

    inline Eigen::Matrix3x4d ToMatrix() const {
        Eigen::Matrix3x4d matrix;
        matrix.leftCols<3>() = rotation.toRotationMatrix();
        matrix.col(3) = translation;
        return matrix;
    }
};


#endif //P3P_REAL_WORLD_EPNP_H
