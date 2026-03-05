#include "Static.hpp"

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <iostream>
#include <string>
#include "utils/matrixUtils.hpp"
#include "utils/EigenTypes.hpp"
#include "mesh/tetmesh.hpp"
#include "energy/energy.hpp"
#include "energy/arap.hpp"
#include "energy/snh.hpp"

// Tools for Static Solving position-based static equilibrium
namespace Solver {

Static::Static(Geom::mesh& m, Geom::mesh& m_start, energy::Energy& e, std::vector<std::vector<Utils::Vector3d>>& c, Utils::Vector3d& a_ext) :
    integrator(m, m_start, e, c, a_ext) {
}

// Line search using a starting alpha 
// returns -1.0 when an error has occurred.
double Static::LineSearch(double start_alpha) {
    double min_alpha = 1e-6;
    bool inversion = true;
    double alpha = start_alpha;

    Eigen::VectorXd dx;
    double residual;
    do {
        residual = TimeStep(dx);

        // Check for inversion
        if (residual < 0.0) {
            alpha = alpha/2.0;
            std::cout << "    - Step size resulted in inverted elements. Reducing step size by half: " << alpha << std::endl;
        } else {   // acceptable step
            std::cout << "    - Step size of " << alpha << " is successful." << std::endl;
            inversion = false;
        }

        // Reduce alpha
        if (alpha <= min_alpha && inversion) {
            std::cout << "    - Step size reduced to minimum (<=" << min_alpha << ") without convergence. Terminating." << std::endl;
            break;
        }
    } while(inversion);
    // Completed one iteration successfully, so we can add to positions
    if (!inversion) {
        // Update positions using the constraints
        std::cout << "    Updating positions" << std::endl;
        // Update positions
        positions_t += alpha * (dx);
        return residual;
    }
    return -1.0;
}

// Take a single Newton step for our solver
double Static::TimeStep(Eigen::VectorXd& dx) {
    Eigen::SparseMatrix<double> _A;
    Eigen::VectorXd _b;
    Eigen::VectorXd _z;
    // Build system matrices
    std::cout << "    Building System" << std::endl;
    double success = buildSystem(_A, _b, _z);
    // If building system failed, then that means there was an element inversion
    if (success < 0.0) {
        return -1.0;
    }

    // Solve system using Eigen's Conjugate Gradient
    std::cout << "    Solving System" << std::endl;
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper> solver;
    solver.compute(_A);
    Eigen::VectorXd _y = solver.solve(_b);

    dx = _y + _z;
    // If we have a bad value, terminate the simulator
    if (!dx.allFinite()) {
        return -1.0;
    }
    // Compute residual
    double res_size = _b.norm();
    return res_size;
}

// Assemble LHS
Eigen::SparseMatrix<double> Static::assembleLHS(Eigen::SparseMatrix<double>& K) {
    return -1 * K;
}

// Compute the force vector (RHS)
Eigen::VectorXd Static::assembleRHS(Eigen::VectorXd& f_ext, Eigen::VectorXd& f_int) {
    return (f_int - f_ext);
}

// Assemble system
// Note we're doing this and returning just the matrices so that I'm not hogging extra memory during the system solve
double Static::buildSystem(Eigen::SparseMatrix<double>& A_global, Eigen::VectorXd& b_global, Eigen::VectorXd& z_global) {
    Eigen::SparseMatrix<double> K(_m.dim*_m.n, _m.dim*_m.n);  // Stiffness Matrix
    Eigen::SparseMatrix<double> C(_m.dim*_m.n, _m.dim*_m.n);  // Damping Matrix
    std::vector<T> tripletK;
    tripletK.reserve((_m.dim*_m.verts_per_el)*(_m.dim*_m.verts_per_el)*_m.num_t);   // (3x4)x(3x4) x tets

    // Eigen::VectorXd f_ext(_m.dim * n);  // External forces are already pre-computed
    Eigen::VectorXd f_int(_m.dim * _m.n);    // Internal elastic force
    f_int.setZero();

    double J_min = 0.05;   // Minimum element size

    // Loop over elements to generate stiffness matrix
    for (int t = 0; t < _m.num_t; t++) {
        // Get information about this element
        std::vector<int> v_idxs = _m.t[t];
        double vol = _m.t_vols[t];
        std::vector<Eigen::Vector3d> verts_new(v_idxs.size());
        // Get current positions of vertices
        for (int v = 0; v < _m.verts_per_el; v++) {
            verts_new[v] = Eigen::Vector3d({positions_t[_m.dim * v_idxs[v]],
                                            positions_t[_m.dim * v_idxs[v]+1],
                                            positions_t[_m.dim * v_idxs[v]+2]});
        }

        // Get the Dm of the element (3x3)
        Eigen::MatrixXd D_mInv = _m.D_mInv[t];
        // Compute Deformation gradient (3x3)
        Utils::Matrix3d F = Utils::computeF(D_mInv, verts_new);
        // If determinant inverts the element
        double J = Utils::computeI3(F);
        if (J <= J_min) {
            return -1.0;
        }
        // Compute dF/dx (3x3)
        Eigen::MatrixXd dFdx = Utils::computedFdx(D_mInv);
        // Compute PK1 for force computation (3x3)
        Utils::Matrix3d PK1 = _e.PK1(F);
        // Compute elastic force itself (12x1)
        Eigen::VectorXd localF = elementElasticForce(vol, dFdx, PK1);
        // Compute Hessian
        // Eigen::MatrixXd Hessian = _e.computeHessian(F);
        // Compute definite version of the Hessian (clamped/projected)
        Utils::Matrix9d clampedHessian = _e.computePSDHessian(F);
        // Multiply Hessian to get force derivative (12x12)
        Eigen::MatrixXd localK = elementElasticGradient(vol, clampedHessian, dFdx);
        // Symmetrize K (should already be symmetrical!)
        // localK = 0.5 * (localK + localK.transpose());
        // Insert into K
        for (int rV = 0; rV < 4; rV++) {
            for (int cV = 0; cV < 4; cV++) {
                for (int rd = 0; rd < 3; rd++) {
                    for (int cd = 0; cd < 3; cd++) {
                        int gr = 3*v_idxs[rV] + rd;
                        int gc = 3*v_idxs[cV] + cd;
                        tripletK.push_back(T(gr, gc, localK(3*rV + rd, 3*cV + cd)));
                    }
                }
            }
        }
        // Insert into elastic forces
        for (int v = 0; v < _m.verts_per_el; v++) {
            for (int d = 0; d < _m.dim; d++) {
                f_int(_m.dim*(v_idxs[v])+d) += (localF(_m.dim*v+d));
            }
        }
    }
    K.setFromTriplets(tripletK.begin(), tripletK.end());

    // Construct A using K and M
    std::cout << "    internal force norm: " << f_int.norm() << std::endl;

    Eigen::SparseMatrix<double> A = assembleLHS(K);
    Eigen::VectorXd b = assembleRHS(f_ext, f_int);

    // Apply Constraints via Prefiltering
    preFilteringA(A, A_global);
    //assembleZ(z_global);
    z_global.resize(b.size());
    z_global.setZero();
    preFilteringb(A, b, z_global, b_global);

    return 1.0;
}

} // namespace Solver
