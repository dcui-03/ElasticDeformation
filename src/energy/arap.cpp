#include "arap.hpp"

#include "utils/EigenTypes.hpp"
#include <Eigen/Core>
#include <iostream>
#include <algorithm>
#include "utils/matrixUtils.hpp"

namespace energy {

// Compute energy itself
// We're using mu as the coeff instead of mu/2
double ARAP::computeEnergy(Utils::Matrix3d& F) {
    double I1 = Utils::computeI1(F);
    double I2 = Utils::computeI2(F);

    return (mu) * (I2 - 2*I1 + 3);
}

// Compute First Piola Kirchoff stress tensor (dPsi/dF)
Utils::Matrix3d ARAP::PK1(Utils::Matrix3d& F) {
    Utils::Matrix3d R; Utils::Matrix3d S; 
    Utils::polarDecomposition(F, R, S);
    return 2*mu*(F - R);
}

// Compute Energy Hessian
Utils::Matrix9d ARAP::computeHessian(Utils::Matrix3d& F) {
    Utils::Matrix9d H1 = Utils::computeH1(F);
    Utils::Matrix9d I; I.setIdentity();
    return 2*mu*(I - H1);
}

// Compute Clamped/Abs Hessian of the Energy 
Utils::Matrix9d ARAP::computePSDHessian(Utils::Matrix3d& F) {
    Utils::Vector9d eigenvals;
    Utils::Matrix9d eigenvecs;

    // Compute rotation variant SVD decomp.
    Utils::Matrix3d U;
    Utils::Vector3d Sigma;
    Utils::Matrix3d V;
    Utils::rotationVariantSVD(F, U, Sigma, V);

    // Compute the scaling and flip eigenvalues
    for (int i = 0; i < 3; i++) {
        eigenvals(i) = 2*mu;
    } for (int i = 6; i < 9; i++) {
        eigenvals(i) = 2*mu;
    }

    // Compute twist eigenvalues (ARAP is a jackpot energy)
    eigenvals(3) = mu * (2 - 4/(Sigma(0) + Sigma(1)));
    eigenvals(4) = mu * (2 - 4/(Sigma(1) + Sigma(2)));
    eigenvals(5) = mu * (2 - 4/(Sigma(2) + Sigma(1)));

    // Get twist and flip eigenvectors
    Utils::buildTwistAndFlipEigenvectors(U, V, eigenvecs);
    // Get scaling eigenvectors
    Utils::buildScalingEigenvectors(U, V, eigenvecs);

    // Zero out the bad eigenvalues
    for (int i = 0; i < eigenvals.size(); i++) {
        if (use_abs) {
            eigenvals(i) = std::abs(eigenvals(i));
        } else {
            eigenvals(i) = std::max(eigenvals(i), 0.0);
        }
    }

    return (eigenvecs * eigenvals.asDiagonal() * eigenvecs.transpose());
}

// Compute energy derivatives of the invariants
double ARAP::EderivativeI1_1() {
    return -2*mu;
} 
double ARAP::EderivativeI1_2() {
    return 0;
}
double ARAP::EderivativeI2_1() {
    return mu;
} 
double ARAP::EderivativeI2_2() {
    return 0;
}
double ARAP::EderivativeI3_1() {
    return 0;
}
double ARAP::EderivativeI3_2() {
    return 0;
}

} // namespace energy