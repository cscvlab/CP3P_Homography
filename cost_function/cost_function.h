#ifndef P3P_REAL_WORLD_COST_FUNCTION_H
#define P3P_REAL_WORLD_COST_FUNCTION_H



#include "../colmap_util/eigen_alignment.h"
#include "../colmap_util/types.h"
#include "../colmap_util/models.h"
#include <Eigen/Core>
#include <ceres/ceres.h>
#include <ceres/conditioned_cost_function.h>
#include <ceres/rotation.h>

template <typename T>
using EigenVector3Map = Eigen::Map<const Eigen::Matrix<T, 3, 1>>;
template <typename T>
using EigenQuaternionMap = Eigen::Map<const Eigen::Quaternion<T>>;
using EigenMatrix6d = Eigen::Matrix<double, 6, 6>;

template <typename CameraModel>
class ReprojErrorCostFunction {
public:
    explicit ReprojErrorCostFunction(const Eigen::Vector2d& point2D)
            : observed_x_(point2D(0)), observed_y_(point2D(1)) {}

    static ceres::CostFunction* Create(const Eigen::Vector2d& point2D) {
        return (
                new ceres::AutoDiffCostFunction<ReprojErrorCostFunction<CameraModel>,
                        2,
                        4,
                        3,
                        3,
                        CameraModel::num_params>(
                        new ReprojErrorCostFunction(point2D)));
    }

    template <typename T>
    bool operator()(const T* const cam_from_world_rotation,
                    const T* const cam_from_world_translation,
                    const T* const point3D,
                    const T* const camera_params,
                    T* residuals) const {
        const Eigen::Matrix<T, 3, 1> point3D_in_cam =
                EigenQuaternionMap<T>(cam_from_world_rotation) *
                EigenVector3Map<T>(point3D) +
                EigenVector3Map<T>(cam_from_world_translation);
        CameraModel::ImgFromCam(camera_params,
                                point3D_in_cam[0],
                                point3D_in_cam[1],
                                point3D_in_cam[2],
                                &residuals[0],
                                &residuals[1]);
        residuals[0] -= T(observed_x_);
        residuals[1] -= T(observed_y_);
        return true;
    }

private:
    const double observed_x_;
    const double observed_y_;
};


template <typename CameraModel>
class ReprojErrorConstantPoint3DCostFunction
        : public ReprojErrorCostFunction<CameraModel> {
    using Parent = ReprojErrorCostFunction<CameraModel>;

public:
    ReprojErrorConstantPoint3DCostFunction(const Eigen::Vector2d& point2D,
                                           const Eigen::Vector3d& point3D)
            : Parent(point2D),
              point3D_x_(point3D(0)),
              point3D_y_(point3D(1)),
              point3D_z_(point3D(2)) {}

    static ceres::CostFunction* Create(const Eigen::Vector2d& point2D,
                                       const Eigen::Vector3d& point3D) {
        return (new ceres::AutoDiffCostFunction<
                ReprojErrorConstantPoint3DCostFunction<CameraModel>,
                2,
                4,
                3,
                CameraModel::num_params>(
                new ReprojErrorConstantPoint3DCostFunction(point2D, point3D)));
    }

    template <typename T>
    bool operator()(const T* const cam_from_world_rotation,
                    const T* const cam_from_world_translation,
                    const T* const camera_params,
                    T* residuals) const {
        const T point3D[3] = {T(point3D_x_), T(point3D_y_), T(point3D_z_)};
        return Parent::operator()(cam_from_world_rotation,
                                  cam_from_world_translation,
                                  point3D,
                                  camera_params,
                                  residuals);
    }

private:
    const double point3D_x_;
    const double point3D_y_;
    const double point3D_z_;
};

template <template <typename> class CostFunction, typename... Args>
ceres::CostFunction* CameraCostFunction(const CameraModelId camera_model_id,
                                        Args&&... args) {
    switch (camera_model_id) {
#define CAMERA_MODEL_CASE(CameraModel)                                     \
  case CameraModel::model_id:                                              \
    return CostFunction<CameraModel>::Create(std::forward<Args>(args)...); \
    break;

        CAMERA_MODEL_SWITCH_CASES

#undef CAMERA_MODEL_CASE
    }
}

inline void SetQuaternionManifold(ceres::Problem* problem, double* quat_xyzw) {
#if CERES_VERSION_MAJOR >= 3 || \
    (CERES_VERSION_MAJOR == 2 && CERES_VERSION_MINOR >= 1)
    problem->SetManifold(quat_xyzw, new ceres::EigenQuaternionManifold);
#else
    problem->SetParameterization(quat_xyzw,
                               new ceres::EigenQuaternionParameterization);
#endif
}

#endif //P3P_REAL_WORLD_COST_FUNCTION_H
