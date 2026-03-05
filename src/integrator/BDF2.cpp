#include "BDF2.hpp"

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

// Tools for BDF2 Velocity-based implicit integration
namespace Solver {

BDF2::BDF2(Geom::mesh& m, Geom::mesh& m_p, Geom::mesh& m_start, energy::Energy& e, double t, bool R, std::vector<std::vector<Utils::Vector3d>>& c, Utils::Vector3d& a_ext):
    integrator(m, m_start, e, t, R, c, a_ext) {
    initializePV_p(m_p);
}

// Initialize previous states
void BDF2::initializePV_p(Geom::mesh& m_p) {
    positions_p.resize(_m.dim * _m.n);
    velocities_p.resize(_m.dim * _m.n);

    for (int v = 0; v < _m.n; v++) {    // Initialize current position to original
        //std::cout << "Vertex " << v << ": " << _m.v[v](0) << ", " << _m.v[v](1) << ", " << _m.v[v](2) << std::endl;
        positions_p(_m.dim*v) = m_p.v[v](0);
        positions_p(_m.dim*v+1) = m_p.v[v](1);
        positions_p(_m.dim*v+2) = m_p.v[v](2);
    }
    velocities_p.setZero();   // Set initial velocities to 0
    return;
}

// Take a single timestep
double BDF2::TimeStep(Eigen::VectorXd& dv) {
    Eigen::SparseMatrix<double> _A;
    Eigen::VectorXd _b;
    Eigen::VectorXd _z;
    // Build system matrices
    std::cout << "    Building System" << std::endl;
    buildSystem(_A, _b, _z);

    // Solve system using Eigen's Conjugate Gradient
    std::cout << "    Solving System" << std::endl;
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper> solver;
    solver.compute(_A);
    Eigen::VectorXd _y = solver.solve(_b);

    // Update velocities using the constraints
    std::cout << "    Updating positions and velocities" << std::endl;
    // Store temporary variables for position and velocity update
    Eigen::VectorXd temp_p(_m.dim*_m.n);
    Eigen::VectorXd temp_v(_m.dim*_m.n);
    for (int vert = 0; vert < _m.dim*_m.n; vert++) {
        temp_p[vert] = positions_t[vert];
        temp_v[vert] = velocities_t[vert];
    }
    // Update current positions and velocities
    velocities_t += (_y + _z);
    // Update positions
    positions_t = (4.0/3.0) * positions_t - (1.0/3.0) * positions_p + (2.0/3.0)*h*velocities_t;
    // Update old positions and velocities
    positions_p = temp_p;
    velocities_p = temp_v;

    // If we have a bad value, terminate the simulator
    if (!velocities_t.allFinite() || !positions_t.allFinite()) {
        return 0.0;
    }
    return 1.0;
}

// Assemble the velocity coffecient (LHS)
Eigen::SparseMatrix<double> BDF2::assembleLHS(Eigen::SparseMatrix<double>& M, Eigen::SparseMatrix<double>& C, Eigen::SparseMatrix<double>& K) {
    return (M - (2.0/3.0)*h*C - (4.0/9.0)*h*h*K);
}

// Compute the force vector (RHS)
Eigen::VectorXd BDF2::assembleRHS(Eigen::VectorXd& f_elast, Eigen::VectorXd& f_damp, Eigen::VectorXd& f_ext, Eigen::SparseMatrix<double>& K) {
    return ((1.0/3.0)*(M*(velocities_t - velocities_p)) + (2.0/3.0)*h*((f_elast + f_damp + f_ext) + (1.0/3.0)*K*(positions_t - positions_p + 2 * h * velocities_t)));
}

// Assemble system
// Note we're doing this and returning just the matrices so that I'm not hogging extra memory during the system solve
double BDF2::buildSystem(Eigen::SparseMatrix<double>& A_global, Eigen::VectorXd& b_global, Eigen::VectorXd& z_global) {
    Eigen::SparseMatrix<double> K(_m.dim*_m.n, _m.dim*_m.n);  // Stiffness Matrix
    Eigen::SparseMatrix<double> C(_m.dim*_m.n, _m.dim*_m.n);  // Damping Matrix
    std::vector<T> tripletK;
    tripletK.reserve((_m.dim*_m.verts_per_el)*(_m.dim*_m.verts_per_el)*_m.num_t);   // (3x4)x(3x4) x tets

    // Eigen::VectorXd f_ext(_m.dim * n);  // External forces are already pre-computed
    Eigen::VectorXd f_damp(_m.dim * _m.n); // Damping force
    Eigen::VectorXd f_elast(_m.dim * _m.n);    // Elastic force
    f_damp.setZero();
    f_elast.setZero();

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
                f_elast(_m.dim*(v_idxs[v])+d) += (localF(_m.dim*v+d));
            }
        }
    }
    K.setFromTriplets(tripletK.begin(), tripletK.end());

    // Construct A using K and M
    if (Rayleigh) { // With Rayleigh Damping
        C = (alpha)*M + (beta)*K;
        f_damp = -1 * (C * velocities_t);
    }   // Else leave it as 0

    Eigen::SparseMatrix<double> A = assembleLHS(M, C, K);
    Eigen::VectorXd b = assembleRHS(f_elast, f_damp, f_ext, K);

    // Apply Constraints via Prefiltering
    preFilteringA(A, A_global);
    assembleZ(z_global);
    preFilteringb(A, b, z_global, b_global);

    return 1.0;
}

} // namespace Solver
