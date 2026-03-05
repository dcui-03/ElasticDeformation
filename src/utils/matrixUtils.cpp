#include "matrixUtils.hpp"

#include "EigenTypes.hpp"
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <iostream>
#include <vector>

namespace Utils {

// Get lambda and mu values based on the Young's Modulus (E) and Poisson Ratio (V)
double solveMu(double E, double V) {
    return (E/(2*(1+V)));
}

double solveLambda(double E, double V) {
    double eps = 1e-5;
    // Extra careful about not exploding at V = 0.5
    if (V <= 0.5+eps && V >= 0.5) {
        V = 0.5+eps;
    } else if (V >= 0.5-eps && V <= 0.5) {
        V = 0.5 - eps;
    }
    return (E*V/((1+V)*(1-2*V)));
}

// Converts a 3n x 1 vector to an nx3 vector
Eigen::MatrixXd foldVector3d(Eigen::VectorXd& p) {
    int numvals = p.size() / 3;
    Eigen::MatrixXd p_folded(numvals, 3);

    for (int i = 0; i < numvals; i++) {
        p_folded(i, 0) = p(3*i);
        p_folded(i, 1) = p(3*i + 1);
        p_folded(i, 2) = p(3*i + 2);
    }
    return p_folded;
}

double computeTetVolume(std::vector<Eigen::Vector3d> tet_verts) {
    Eigen::Vector3d diff1 = tet_verts[1] - tet_verts[0];
    Eigen::Vector3d diff2 = tet_verts[2] - tet_verts[0];
    Eigen::Vector3d diff3 = tet_verts[3] - tet_verts[0];
    double vol = ((diff3.dot((diff1).cross(diff2)) / 6.0));
    return vol;
}

// Cross product matrix (little hat)
Utils::Matrix3d crossProdMatrix(Utils::Vector3d x) {
    Utils::Matrix3d x_hat;
    x_hat << 0, -x(2), x(1),
             x(2), 0, -x(0),
             -x(1), x(0), 0;
    return x_hat;
}

// Vectorize a matrix
Utils::Vector9d vectorizeMatrix(Utils::Matrix3d& mat) {
    Utils::Vector9d vec;
    int counter = 0;
    for (int c = 0; c < 3; c++) {
        for (int r = 0; r < 3; r++) {
            vec(counter) = mat(r, c);
            counter++;
        }
    }
    return vec;
}

// Matricize a 9d vector
Utils::Matrix3d matricizeVector9d(Utils::Vector9d& vec) {
    Utils::Matrix3d mat; mat.setZero();
    int counter = 0;
    for (int c = 0; c < 3; c++) {
        for (int r = 0; r < 3; r++) {
            mat(r, c) = vec(counter);
            counter++;
        }
    }
    return mat;
}

// Compute the Rotation Variant SVD
// Input an empty U, Sigma, V
void rotationVariantSVD(Utils::Matrix3d& mat, Utils::Matrix3d& U, Utils::Vector3d& Sigma, Utils::Matrix3d& V) {
    // Compute SVD to get Sigma, U and V
    Eigen::JacobiSVD<Utils::Matrix3d> svd(mat, Eigen::ComputeFullU | Eigen::ComputeFullV);
    U = svd.matrixU();
    Sigma = svd.singularValues();
    V = svd.matrixV();

    // Compute L and remove reflections from U and V
    Utils::Matrix3d L;
    L.setIdentity();
    L(2,2) = (U*V.transpose()).determinant();
    Sigma(2) = Sigma(2) * L(2,2); // To keep it a vector (taken from HOBAK)

    double u_det = U.determinant();
    double v_det = V.determinant();
    if (u_det < 0 && v_det > 0) {
        U = U * L;
    } else if (u_det > 0 && v_det < 0) {
        V = V * L;
    }

    return;
}

// Compute Polar Decomposition given rotation variant SVD
// Returns a vector containing R and then S
void polarDecomposition(Utils::Matrix3d& mat, Utils::Matrix3d& R, Utils::Matrix3d& S) {
    Utils::Matrix3d U, V; Utils::Vector3d Sigma;
    rotationVariantSVD(mat, U, Sigma, V);

    // Put together R and S from the inputs
    R = U * V.transpose();
    S = V * Sigma.asDiagonal() * V.transpose();

    return;
}

// Compute the deformation gradient given original positions p and new positions q
Utils::Matrix3d computeF(Utils::Matrix3d Dm_Inv, std::vector<Utils::Vector3d> tet_new) {
    Utils::Matrix3d Ds;
    Ds.setZero();
    Ds.col(0) = tet_new[1] - tet_new[0];
    Ds.col(1) = tet_new[2] - tet_new[0];
    Ds.col(2) = tet_new[3] - tet_new[0];

    return (Ds * Dm_Inv);
}

// Compute F w.r.t new positions
// Gave up and just followed HOBAK
Eigen::MatrixXd computedFdx(const Utils::Matrix3d& Dm_inv) {
    Eigen::MatrixXd dFdx(9, 12);
    dFdx.setZero();

    double m = Dm_inv (0, 0);
    double n = Dm_inv (0, 1);
    double o = Dm_inv (0, 2);
    double p = Dm_inv (1, 0);
    double q = Dm_inv (1, 1);
    double r = Dm_inv (1, 2);
    double s = Dm_inv (2, 0);
    double t = Dm_inv (2, 1);
    double u = Dm_inv (2, 2);

    double t1 = - m - p - s;
    double t2 = - n - q - t;
    double t3 = - o - r - u;

    dFdx(0, 0) = t1;
    dFdx(0, 3) = m;
    dFdx(0, 6) = p;
    dFdx(0, 9) = s;
    dFdx(1, 1) = t1;
    dFdx(1, 4) = m;
    dFdx(1, 7) = p;
    dFdx(1, 10) = s;
    dFdx(2, 2) = t1;
    dFdx(2, 5) = m;
    dFdx(2, 8) = p;
    dFdx(2, 11) = s;
    dFdx(3, 0) = t2;
    dFdx(3, 3) = n;
    dFdx(3, 6) = q;
    dFdx(3, 9) = t;
    dFdx(4, 1) = t2;
    dFdx(4, 4) = n;
    dFdx(4, 7) = q;
    dFdx(4, 10) = t;
    dFdx(5, 2) = t2;
    dFdx(5, 5) = n;
    dFdx(5, 8) = q;
    dFdx(5, 11) = t;
    dFdx(6, 0) = t3;
    dFdx(6, 3) = o;
    dFdx(6, 6) = r;
    dFdx(6, 9) = u;
    dFdx(7, 1) = t3;
    dFdx(7, 4) = o;
    dFdx(7, 7) = r;
    dFdx(7, 10) = u;
    dFdx(8, 2) = t3;
    dFdx(8, 5) = o;
    dFdx(8, 8) = r;
    dFdx(8, 11) = u;

    return dFdx;
}

// Compute the oh so painful rotation gradient that ARAP wants
// Equal to second invariant's first derivative (g2)!
Utils::Matrix9d computeRotationGradient(Utils::Matrix3d& F) {
    // pg. 69
    Eigen::JacobiSVD<Utils::Matrix3d> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Utils::Matrix3d U = svd.matrixU();
    Utils::Vector3d Sigma = svd.singularValues();
    Utils::Matrix3d V = svd.matrixV();

    return computeH1(U, Sigma, V);
}

// Stole this function from HOBAK (MATRIX_UTIL.CPP)
// This is the pre-vectorized form of gj
Utils::Matrix3d partialJpartialF(Utils::Matrix3d& F) {
    Utils::Matrix3d pJpF;
    pJpF.col(0) = F.col(1).cross(F.col(2));
    pJpF.col(1) = F.col(2).cross(F.col(0));
    pJpF.col(2) = F.col(0).cross(F.col(1));
    return pJpF;
}

// Twist and Flip eigenmatrices for getting the projected Hessian
// Taken in part from HOBAK
void buildTwistAndFlipEigenvectors(const Utils::Matrix3d& U, const Utils::Matrix3d& V, Utils::Matrix9d& Q) {
    // There's no pretty way of doing this, I'm afraid...
    // Twist Eigenvectors
    Utils::Matrix3d T3; Utils::Matrix3d T4; Utils::Matrix3d T5;
    T3 << 0, -1, 0,
          1, 0, 0,
          0, 0, 0;
    T4 << 0, 0, 0,
          0, 0, 1,
          0, -1, 0;
    T5 << 0, 0, 1,
          0, 0, 0,
          -1, 0, 0;

    Utils::Matrix3d Q3 = (1/sqrt(2))* (U * T3 * V.transpose());
    Utils::Matrix3d Q4 = (1/sqrt(2))* (U * T4 * V.transpose());
    Utils::Matrix3d Q5 = (1/sqrt(2))* (U * T5 * V.transpose());
    // Flip Eigenvectors
    Utils::Matrix3d T6; Utils::Matrix3d T7; Utils::Matrix3d T8;
    T6 << 0, 1, 0,
          1, 0, 0,
          0, 0, 0;
    T7 << 0, 0, 0,
          0, 0, 1,
          0, 1, 0;
    T8 << 0, 0, 1,
          0, 0, 0,
          1, 0, 0;

    Utils::Matrix3d Q6 = (1/sqrt(2))* (U * T6 * V.transpose());
    Utils::Matrix3d Q7 = (1/sqrt(2))* (U * T7 * V.transpose());
    Utils::Matrix3d Q8 = (1/sqrt(2))* (U * T8 * V.transpose());

    Q.col(3) = vectorizeMatrix(Q3);
    Q.col(4) = vectorizeMatrix(Q4);
    Q.col(5) = vectorizeMatrix(Q5);
    Q.col(6) = vectorizeMatrix(Q6);
    Q.col(7) = vectorizeMatrix(Q7);
    Q.col(8) = vectorizeMatrix(Q8);
    return;
}

// Scaling eigenmatrices for getting the projected Hessian
// Taken in part from HOBAK
void buildScalingEigenvectors(const Utils::Matrix3d& U, const Utils::Matrix3d& eigenvecs, const Utils::Matrix3d& V, Utils::Matrix9d& Q) {
    // Instead of using the D matrices as in the text,
    // we know the eigenvals are related to the eigenvecs of the submatrix A,
    // so just use those directly... way easier!
    Utils::Matrix3d Q0 = U * (eigenvecs.col(0)).asDiagonal() * V.transpose();
    Utils::Matrix3d Q1 = U * (eigenvecs.col(1)).asDiagonal() * V.transpose();
    Utils::Matrix3d Q2 = U * (eigenvecs.col(2)).asDiagonal() * V.transpose();

    Q.col(0) = vectorizeMatrix(Q0);
    Q.col(1) = vectorizeMatrix(Q1);
    Q.col(2) = vectorizeMatrix(Q2);
    return;
}

// A jackpot version where we don't need to build submatrix A
void buildScalingEigenvectors(const Utils::Matrix3d& U, const Utils::Matrix3d& V, Utils::Matrix9d& Q) {
    Utils::Vector3d x = {1,0,0};
    Utils::Vector3d y = {0,1,0};
    Utils::Vector3d z = {0,0,1};
    // Since off diagonals are 0, we can just use this directly
    Utils::Matrix3d Q0 = U * x.asDiagonal() * V.transpose();
    Utils::Matrix3d Q1 = U * y.asDiagonal() * V.transpose();
    Utils::Matrix3d Q2 = U * z.asDiagonal() * V.transpose();

    Q.col(0) = vectorizeMatrix(Q0);
    Q.col(1) = vectorizeMatrix(Q1);
    Q.col(2) = vectorizeMatrix(Q2);
    return;
}

// COMPUTE INVARIANTS, GRADIENTS, AND THEIR HESSIANS
double computeI1(Utils::Matrix3d& F) {
    // Polar Decomposition
    Utils::Matrix3d R; Utils::Matrix3d S;
    polarDecomposition(F, R, S);
    return S.trace();
}

double computeI2(Utils::Matrix3d& F) {
    Utils::Matrix3d I_C = F.transpose() * F;
    return I_C.trace();
}

double computeI3(Utils::Matrix3d& F) {  // i.e., J
    return F.determinant();
}

Utils::Vector9d computeg1(Utils::Matrix3d F) {
    // Polar Decomposition
    Utils::Matrix3d R; Utils::Matrix3d S;
    polarDecomposition(F, R, S);
    return vectorizeMatrix(R);
}

Utils::Vector9d computeg2(Utils::Matrix3d& F) {
    Utils::Matrix3d temp = 2*F;
    return vectorizeMatrix(temp);
}

Utils::Vector9d computegj(Utils::Matrix3d& F) {
    Utils::Matrix3d pJpF = partialJpartialF(F);
    return vectorizeMatrix(pJpF);
}

Utils::Matrix9d computeH1(Utils::Matrix3d& U, Utils::Vector3d& Sigma, Utils::Matrix3d& V) {
    Utils::Matrix3d temp0; Utils::Matrix3d temp1; Utils::Matrix3d temp2;
    temp0 << 0, -1, 0,     1, 0, 0,    0, 0, 0;
    temp1 << 0, 0, 0,     0, 0, 1,    0, -1, 0;
    temp2 << 0, 0, 1,     0, 0, 0,    -1, 0, 0;
    Utils::Matrix3d Q0 = U * temp0 * V.transpose();
    Utils::Matrix3d Q1 = U * temp1 * V.transpose();
    Utils::Matrix3d Q2 = U * temp2 * V.transpose();

    double eigen0 = 2/(Sigma(0) + Sigma(1));
    double eigen1 = 2/(Sigma(1) + Sigma(2));
    double eigen2 = 2/(Sigma(2) + Sigma(0));

    // Assemble matrix
    Utils::Matrix9d H1;
    H1.setZero();
    H1 += eigen0 * (vectorizeMatrix(Q0) * vectorizeMatrix(Q0).transpose());
    H1 += eigen1 * (vectorizeMatrix(Q1) * vectorizeMatrix(Q1).transpose());
    H1 += eigen2 * (vectorizeMatrix(Q2) * vectorizeMatrix(Q2).transpose());
    return H1;
}

Utils::Matrix9d computeH1(Utils::Matrix3d& F) {
    Utils::Matrix3d U; Utils::Vector3d Sigma; Utils::Matrix3d V;
    rotationVariantSVD(F, U, Sigma, V);
    return computeH1(U, Sigma, V);
}

Utils::Matrix9d computeH2(Utils::Matrix3d& F) {
    return 2 * Utils::Matrix9d::Identity(9, 9);
}

Utils::Matrix9d computeHj(Utils::Matrix3d& F) {
        Utils::Matrix9d Hj;
        Hj.setZero();

        Utils::Matrix3d f0 = crossProdMatrix(F.col(0));
        Utils::Matrix3d f1 = crossProdMatrix(F.col(1));
        Utils::Matrix3d f2 = crossProdMatrix(F.col(2));
        // Insert into Hj (start index, matrix size)
        Hj.block(0, 3, 3, 3) = -f2;
        Hj.block(0, 6, 3, 3) = f1;
        Hj.block(3, 0, 3, 3) = f2;
        Hj.block(3, 6, 3, 3) = -f0;
        Hj.block(6, 0, 3, 3) = -f1;
        Hj.block(6, 3, 3, 3) = f0;

        return Hj;
    }

} // namespace Utils