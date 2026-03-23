#include "energy.hpp"

#include "utils/EigenTypes.hpp"
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <iostream>
#include <algorithm>
#include "utils/matrixUtils.hpp"

namespace energy {

Energy::Energy(double _mu, double _lambda, bool abs): mu(_mu), lambda(_lambda), use_abs(abs) {
}

Energy::~Energy() = default;

// COMPUTE ENERGY
double Energy::computeEnergy(Utils::Matrix3d& F) {   // Dirichlet placeholder
    return (F.transpose() * F).trace();
}

// COMPUTE PK1 OF ENERGY
Utils::Matrix3d Energy::PK1(Utils::Matrix3d& F) {
    // Dirichlet Placeholder
    return 2*F;
}

// COMPUTE ENERGY HESSIAN
Utils::Matrix9d Energy::computeHessian(Utils::Matrix3d& F) {
    // Compute invariants
    double I1 = Utils::computeI1(F);
    double I2 = Utils::computeI2(F);
    double I3 = Utils::computeI3(F);

    Utils::Vector9d g1 = Utils::computeg1(F);
    Utils::Vector9d g2 = Utils::computeg2(F);
    Utils::Vector9d gj = Utils::computegj(F);

    Utils::Matrix9d H1 = Utils::computeRotationGradient(F);
    Utils::Matrix9d H2 = Utils::computeH2(F);
    Utils::Matrix9d Hj = Utils::computeHj(F);

    // Compute energy derivatives
    double eI1_1 = EderivativeI1_1();
    double eI1_2 = EderivativeI1_2();
    double eI2_1 = EderivativeI2_1();
    double eI2_2 = EderivativeI2_2();
    double eI3_1 = EderivativeI3_1();
    double eI3_2 = EderivativeI3_2();

    // Assemble Hessian
    Utils::Matrix9d Hessian;
    Hessian.setZero();
    Hessian += (eI1_2 * (g1 * g1.transpose())) + (eI1_1 * H1);
    Hessian += (eI2_2 * (g2 * g2.transpose())) + (eI2_1 * H2);
    Hessian += (eI3_2 * (gj * gj.transpose())) + (eI3_1 * Hj);

    return Hessian;
}

// COMPUTE CLAMPED ENERGY HESSIAN (For PSD)
Utils::Matrix9d Energy::computePSDHessian(Utils::Matrix3d& F) {
    Utils::Vector9d eigenvals;
    Utils::Matrix9d eigenvecs;

    // Compute rotation variant SVD decomp.
    Utils::Matrix3d U;
    Utils::Vector3d Sigma;
    Utils::Matrix3d V;
    Utils::rotationVariantSVD(F, U, Sigma, V);
    // Compute I1, I2, I3
    double I1 = Utils::computeI1(F);
    double I2 = Utils::computeI2(F);
    double I3 = Utils::computeI3(F);

    // Compute the twist and flip eigenvalues
    eigenvals(3) = (2/(Sigma(0) + Sigma(1))) * EderivativeI1_1() + 2 * EderivativeI2_1() + (Sigma(2) * EderivativeI3_1());
    eigenvals(4) = (2/(Sigma(1) + Sigma(2))) * EderivativeI1_1() + 2 * EderivativeI2_1() + (Sigma(0) * EderivativeI3_1());
    eigenvals(5) = (2/(Sigma(0) + Sigma(2))) * EderivativeI1_1() + 2 * EderivativeI2_1() + (Sigma(1) * EderivativeI3_1());
    eigenvals(6) = 2 * EderivativeI2_1() - (Sigma(2) * EderivativeI3_1());
    eigenvals(7) = 2 * EderivativeI2_1() - (Sigma(0) * EderivativeI3_1());
    eigenvals(8) = 2 * EderivativeI2_1() - (Sigma(1) * EderivativeI3_1());

    // Compute the scaling eigenvalue matrix
    Utils::Matrix3d A; A.setZero();
    for (int i = 0; i < 3; i++) {
        // Gross, copied this from the notes
        A(i,i) = 2*EderivativeI2_1() + EderivativeI1_2() + 4*Sigma(i)*Sigma(i)*EderivativeI2_2() +
                 (I3*I3/(Sigma(i)*Sigma(i)))*EderivativeI3_2() + 4*Sigma(i)*EderivativeI1I2() +
                 4*I3*EderivativeI2I3() + 2*(I3/Sigma(i))*EderivativeI3I1();
        for (int j = i+1; j < 3; j++) {
            int k;
            if (i == 0 && j == 1) {k = 2;};
            if (i == 0 && j == 2) {k = 1;};
            if (i == 1 && j == 2) {k = 0;};
            double temp = Sigma(k)*EderivativeI3_1() + EderivativeI1_2() + (4*I3/Sigma(k))*EderivativeI2_2() +
                          Sigma(k)*I3*EderivativeI3_2() + 2*Sigma(k)*(I2 - Sigma(k)*Sigma(k))*EderivativeI2I3() + 
                          (I1 - Sigma(k))*(Sigma(k)*EderivativeI3I1() + 2*EderivativeI1I2());
            A(i, j) = temp;
            A(j, i) = temp;
        }
    }
    // Compute scaling eigenvalues
    Eigen::SelfAdjointEigenSolver<Utils::Matrix3d> eigensys(A);
    eigenvals(0) = eigensys.eigenvalues()[0];
    eigenvals(1) = eigensys.eigenvalues()[1];
    eigenvals(2) = eigensys.eigenvalues()[2];

    // Get twist and flip eigenvectors
    Utils::buildTwistAndFlipEigenvectors(U, V, eigenvecs);
    // Get scaling eigenvectors
    Utils::buildScalingEigenvectors(U, eigensys.eigenvectors(), V, eigenvecs);

    // Zero out or absolute value the bad eigenvalues
    for (int i = 0; i < eigenvals.size(); i++) {
        if (use_abs) {
            eigenvals(i) = std::abs(eigenvals(i));
        } else {
            eigenvals(i) = std::max(eigenvals(i), 0.0);
        }
    }

    return (eigenvecs * eigenvals.asDiagonal() * eigenvecs.transpose());
}

// COMPUTE DERIVATIVES OF ENERGY W.R.T. INVARIANTS
double Energy::EderivativeI1_1() {
    return 0;
} 
double Energy::EderivativeI1_2() {
    return 0;
} 
double Energy::EderivativeI2_1() {  // For dirichlet this is 1
    return 1;
}
double Energy::EderivativeI2_2() {
    return 0;
}
double Energy::EderivativeI3_1() {
    return 0;
}
double Energy::EderivativeI3_2() {
    return 0;
}

// Note that for Dirichlet these are all 0
double Energy::EderivativeI1I2() {
    return 0;
}
double Energy::EderivativeI2I3() {
    return 0;
}
double Energy::EderivativeI3I1() {
    return 0;
}

} // namespace energy