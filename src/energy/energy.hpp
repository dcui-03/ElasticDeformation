// Energy.hpp
#pragma once

#include "utils/EigenTypes.hpp"
#include <Eigen/Core>
#include <string>

// A generic class for energies. When unspecified, uses Dirichlet energy
namespace energy {

class Energy {
    public:
        Energy(double _mu, double _lambda);
        virtual ~Energy();

        // Compute First Piola Kirchoff stress tensor
        virtual Utils::Matrix3d PK1(Utils::Matrix3d& F);

        // Compute Energy Hessian (generic)
        virtual Utils::Matrix9d computeHessian(Utils::Matrix3d& F);
        // Compute Clamped Hessian (for PSD)
        virtual Utils::Matrix9d computeClampedHessian(Utils::Matrix3d& F);
    protected:
        // Compute energy itself
        virtual double computeEnergy(Utils::Matrix3d& F);

        // Compute energy derivatives of the invariants
        virtual double EderivativeI1_1(); 
        virtual double EderivativeI1_2(); 
        virtual double EderivativeI2_1(); 
        virtual double EderivativeI2_2(); 
        virtual double EderivativeI3_1();
        virtual double EderivativeI3_2();
        
        // The three mixed double derivatives
        virtual double EderivativeI1I2();
        virtual double EderivativeI2I3();
        virtual double EderivativeI3I1();

        // material parameters
        double mu; double lambda;

    private:
        std::string energy_type = "Dirichlet";
};

} // namespace energy