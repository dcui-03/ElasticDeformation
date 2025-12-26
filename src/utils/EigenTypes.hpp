// EigenTypes.hpp
#pragma once

#include <Eigen/Core>

namespace Utils {

    using Matrix3d = Eigen::Matrix3d;
    using Vector3d = Eigen::Vector3d;

    using Matrix9d = Eigen::Matrix<double, 9, 9>;
    using Vector9d = Eigen::Matrix<double, 9, 1>;

}