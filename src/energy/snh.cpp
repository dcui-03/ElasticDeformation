#include "snh.hpp"

#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <iostream>
#include <algorithm>
#include "utils/matrixUtils.hpp"
#include "utils/EigenTypes.hpp"


namespace energy {

// Constructor
SNH::SNH(double _mu, double _lambda): Energy(_mu, _lambda) {
    alpha = 1 + (mu/lambda);
}

// Compute energy itself
double SNH::computeEnergy(Utils::Matrix3d& F) {
    double I2 = Utils::computeI2(F);
    double J = Utils::computeI3(F);

    return (mu/2) * (I2 - 3) - (mu * (J - 1)) + ((lambda/2) * (J-1)*(J-1));
}

// Compute First Piola Kirchoff stress tensor (dPsi/dF)
Utils::Matrix3d SNH::PK1(Utils::Matrix3d& F) {
    Utils::Matrix3d dJdF = Utils::partialJpartialF(F);
    double J = Utils::computeI3(F);
    return (mu * F) + (lambda * (J - alpha) * dJdF);
}

// Compute Energy Hessian
Utils::Matrix9d SNH::computeHessian(Utils::Matrix3d& F) {
    // First compute I3, its gradient, and its Hessian
    double J = Utils::computeI3(F);
    Utils::Matrix3d dJdF = Utils::partialJpartialF(F);
    Utils::Matrix9d Hj = Utils::computeHj(F);
    Utils::Matrix9d H2 = Eigen::MatrixXd::Identity(9, 9);

    return (mu*H2) + lambda*(Utils::vectorizeMatrix(dJdF) * Utils::vectorizeMatrix(dJdF).transpose()) + lambda*(J - alpha)*Hj;
}

// Compute Clamped Hessian of the Energy 
Utils::Matrix9d SNH::computeClampedHessian(Utils::Matrix3d& F) {
    Utils::Vector9d eigenvals;
    Utils::Matrix9d eigenvecs;

    // Compute rotation variant SVD decomp.
    Utils::Matrix3d U;
    Utils::Vector3d Sigma;
    Utils::Matrix3d V;
    Utils::rotationVariantSVD(F, U, Sigma, V);
    // Compute I3
    double I3 = Utils::computeI3(F);

    // Compute the twist and flip eigenvalues
    double eigen_temp = (lambda*(I3 - 1) - mu);
    eigenvals(3) = mu + Sigma(2)*eigen_temp;
    eigenvals(4) = mu + Sigma(0)*eigen_temp;
    eigenvals(5) = mu + Sigma(1)*eigen_temp;
    eigenvals(6) = mu - Sigma(2)*eigen_temp;
    eigenvals(7) = mu - Sigma(0)*eigen_temp;
    eigenvals(8) = mu - Sigma(1)*eigen_temp;

    // Compute the scaling eigenvalue matrix
    Utils::Matrix3d A; A.setZero();
    for (int i = 0; i < 3; i++) {
        A(i,i) = mu + lambda * (I3*I3/(Sigma(i)*Sigma(i)));
        for (int j = i+1; j < 3; j++) {
            int k;
            if (i == 0 && j == 1) {k = 2;};
            if (i == 0 && j == 2) {k = 1;};
            if (i == 1 && j == 2) {k = 0;};
            double temp = Sigma(k) * (lambda * (2*I3 - 1) - mu);
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

    // Zero out the bad eigenvalues
    for (int i = 0; i < eigenvals.size(); i++) {
        eigenvals(i) = std::max(eigenvals(i), 0.0);
    }

    return (eigenvecs * eigenvals.asDiagonal() * eigenvecs.transpose());
}

// Compute energy derivatives of the invariants
double SNH::EderivativeI1_1() {
    return 0;
} 
double SNH::EderivativeI1_2() {
    return 0;
}
double SNH::EderivativeI2_1() {
    return mu/2;
} 
double SNH::EderivativeI2_2() {
    return 0;
}
double SNH::EderivativeI3_1(Utils::Matrix3d F) {
    double J = Utils::computeI3(F);
    return (lambda * J) - (lambda + mu);
}
double SNH::EderivativeI3_2(Utils::Matrix3d F) {
    return lambda;
}

} // namespace energy