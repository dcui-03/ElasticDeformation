// snh.hpp
#pragma once

#include "energy.hpp"
#include <Eigen/Core>
#include <string>
#include "utils/EigenTypes.hpp"

namespace energy {

class SNH : public Energy {
    public:
        SNH(double _mu, double _lambda);

        // Compute energy itself
        double computeEnergy(Utils::Matrix3d& F) override;
        // Compute First Piola Kirchoff stress tensor
        Utils::Matrix3d PK1(Utils::Matrix3d& F) override;

        // Compute Energy Hessian
        Utils::Matrix9d computeHessian(Utils::Matrix3d& F) override;
        // Compute Clamped Hessian of the Energy 
        Utils::Matrix9d computeClampedHessian(Utils::Matrix3d& F) override;
    protected:

        // Compute energy derivatives of the invariants
        double EderivativeI1_1() override; 
        double EderivativeI1_2() override; 
        double EderivativeI2_1() override;
        double EderivativeI2_2() override; 
        double EderivativeI3_1(Utils::Matrix3d F);  // Already overrides parent
        double EderivativeI3_2(Utils::Matrix3d F);
    private:
        // "Chubbiness"
        double alpha;
        std::string energy_type = "SNH";
};

}   // namespace energy