// matrixUtils.hpp
#pragma once

#include "EigenTypes.hpp"
#include <Eigen/Core>
#include <Eigen/SVD>
#include <vector>

namespace Utils {

// Get lambda and mu values based on the Young's Modulus (E) and Poisson Ratio (V)
double solveMu(double E, double V);

double solveLambda(double E, double V);

// Cross product matrix (little hat)
Utils::Matrix3d crossProdMatrix(Utils::Vector3d x);

// Vectorize a matrix
Utils::Vector9d vectorizeMatrix(Utils::Matrix3d& mat);

// Matricize a 9d vector
Utils::Matrix3d matricizeVector9d(Utils::Vector9d& vec);

// Compute the Rotation Variant SVD
// Input an empty U, Sigma, V
void rotationVariantSVD(Utils::Matrix3d& mat, Utils::Matrix3d& U, Utils::Vector3d& Sigma, Utils::Matrix3d& V);

// Compute Polar Decomposition given rotation variant SVD
// Returns a vector containing R and then S
void polarDecomposition(Utils::Matrix3d& mat, Utils::Matrix3d& R, Utils::Matrix3d& S);

// Compute the deformation gradient given original positions p and new positions q
Utils::Matrix3d computeF(Utils::Matrix3d Dm_Inv, std::vector<Utils::Vector3d> tet_new);

// Compute F w.r.t new positions
// Gave up and just followed HOBAK
Eigen::MatrixXd computedFdx(const Utils::Matrix3d& Dm_inv);

// Compute the oh so painful rotation gradient that ARAP wants
// Equal to second invariant's first derivative (g2)!
Utils::Matrix9d computeRotationGradient(Utils::Matrix3d& F);

// Stole this function from HOBAK (MATRIX_UTIL.CPP)
// This is the pre-vectorized form of gj
Utils::Matrix3d partialJpartialF(Utils::Matrix3d& F);

// Twist and Flip eigenmatrices for getting the projected Hessian
// Taken in part from HOBAK
void buildTwistAndFlipEigenvectors(const Utils::Matrix3d& U, const Utils::Matrix3d& V, Utils::Matrix9d& Q);

// Scaling eigenmatrices for getting the projected Hessian
// Taken in part from HOBAK
void buildScalingEigenvectors(const Utils::Matrix3d& U, const Utils::Matrix3d& eigenvecs, const Utils::Matrix3d& V, Utils::Matrix9d& Q);

// A jackpot version where we don't need to build submatrix A
void buildScalingEigenvectors(const Utils::Matrix3d& U, const Utils::Matrix3d& V, Utils::Matrix9d& Q);

// COMPUTE INVARIANTS, GRADIENTS, AND THEIR HESSIANS
double computeI1(Utils::Matrix3d& F);

double computeI2(Utils::Matrix3d& F);

double computeI3(Utils::Matrix3d& F);

Utils::Vector9d computeg1(Utils::Matrix3d F);

Utils::Vector9d computeg2(Utils::Matrix3d& F);

Utils::Vector9d computegj(Utils::Matrix3d& F);

Utils::Matrix9d computeH1(Utils::Matrix3d& U, Utils::Vector3d& Sigma, Utils::Matrix3d& V);
Utils::Matrix9d computeH1(Utils::Matrix3d& F);

Utils::Matrix9d computeH2(Utils::Matrix3d& F);

Utils::Matrix9d computeHj(Utils::Matrix3d& F);

} // namespace Utils