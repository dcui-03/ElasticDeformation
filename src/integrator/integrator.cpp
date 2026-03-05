#include "integrator.hpp"

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <iostream>
#include <string>
#include "utils/matrixUtils.hpp"
#include "utils/EigenTypes.hpp"
#include "mesh/mesh.hpp"
#include "mesh/tetmesh.hpp"
#include "energy/energy.hpp"
#include "energy/arap.hpp"
#include "energy/snh.hpp"

// A generic implicit integration template

namespace Solver {

// Constructor
integrator::integrator(Geom::mesh& m, Geom::mesh& m_start, energy::Energy& e, double t, bool R, std::vector<std::vector<Utils::Vector3d>>& c, Utils::Vector3d& a_ext):
        _m(m), _e(e), h(t), Rayleigh(R), constraints(c) {
    initializePV(m_start);
    initializeDamping();
    initializeMass();
    initializeExtForce(a_ext);
    assemblePrefiltering();
};

// Constructor for Solves without timestepping or Rayleigh Damping
integrator::integrator(Geom::mesh& m, Geom::mesh& m_start, energy::Energy& e, std::vector<std::vector<Utils::Vector3d>>& c, Utils::Vector3d& a_ext):
        _m(m), _e(e), constraints(c) {
    initializePV(m_start);
    initializeMass();
    initializeExtForce(a_ext);
    assemblePrefiltering();
};

// Set initial velocities to 0 and positions to the original
void integrator::initializePV(Geom::mesh& m_start) {
    positions_t.resize(_m.dim * _m.n);
    velocities_t.resize(_m.dim * _m.n);

    for (int v = 0; v < _m.n; v++) {    // Initialize current position to original
        //std::cout << "Vertex " << v << ": " << _m.v[v](0) << ", " << _m.v[v](1) << ", " << _m.v[v](2) << std::endl;
        positions_t(_m.dim*v) = m_start.v[v](0);
        positions_t(_m.dim*v+1) = m_start.v[v](1);
        positions_t(_m.dim*v+2) = m_start.v[v](2);
    }
    velocities_t.setZero();   // Set initial velocities to 0
    return;
}

// Set damping, if so desired
void integrator::initializeDamping() {
    // Set alpha and beta here
    return;
}

// Precompute the mass matrix
void integrator::initializeMass() {
    M.resize(_m.dim*_m.n, _m.dim*_m.n);
    std::vector<T> tripletM;
    tripletM.reserve(_m.dim*_m.n);
    // Iterate over elements
    for (int t = 0; t < _m.num_t; t++) {
        double vert_vol = _m.t_vols[t]/_m.verts_per_el;

        for (int i = 0; i < _m.verts_per_el; i++) {
            for (int d = 0; d < _m.dim; d++) { // Each component of each vertex gets the same value
                tripletM.push_back(T(_m.dim * _m.t[t][i] + d, _m.dim * _m.t[t][i] + d, vert_vol));
            }
        }
    }
    M.setFromTriplets(tripletM.begin(), tripletM.end());
    return;
}

// Precompute the external force vector
void integrator::initializeExtForce(Utils::Vector3d& a_ext) {
    f_ext.resize(_m.dim*_m.n);
    f_ext.setZero();
    if (a_ext.size() != _m.dim) {
        std::cout << "External force acceleration is incompatible with problem dimension." << std::endl;
        return;
    }
    // Iterate over elements
    for (int t = 0; t < _m.num_t; t++) {
        double vert_vol = _m.t_vols[t];
        // Compute force as volume * acceleration
        Eigen::VectorXd local_fext = vert_vol * a_ext;
        
        // Insert each dim into external force vector
        for (int i = 0; i < _m.t[t].size(); i++) {
            for (int d = 0; d < _m.dim; d++) {
                f_ext(_m.dim * _m.t[t][i] + d) += local_fext(d)/_m.verts_per_el;
            }
        }
    }
    return;
}

// Compute local S matrix
Eigen::Matrix3d integrator::localS3D(const std::vector<Utils::Vector3d> c) {
    Eigen::Matrix3d I;
    I.setIdentity();

    if (c.size() == 0) {
        return I;
    } else if (c.size() >= _m.dim) {
        return I.setZero();
    }

    for (int i = 0; i < c.size(); i++) {
        I -= (c[i]) * (c[i]).transpose();
    }

    return I;
}

// Assemble prefiltering matrix S for PPCG
void integrator::assemblePrefiltering() {
    S.resize(_m.dim*_m.n, _m.dim*_m.n);
    std::vector<T> tripletS;
    tripletS.reserve(_m.dim*_m.dim*_m.n);

    IminusS.resize(_m.dim*_m.n, _m.dim*_m.n);
    // Iterate over vertices
    for (int v = 0; v < _m.n; v++) {
        Eigen::MatrixXd S_local = localS3D(constraints[v]);

        for (int r = 0; r < _m.dim; r++) {
            for (int c = 0; c < _m.dim; c++) {
                tripletS.push_back(T(_m.dim*v+r, _m.dim*v+c, S_local(r, c)));
            }
        }
    }
    S.setFromTriplets(tripletS.begin(), tripletS.end());
    // Also compute I - S here
    Eigen::SparseMatrix<double> I(_m.dim*_m.n, _m.dim*_m.n);
    I.setIdentity();
    IminusS = I - S;
    return;
}

// Line search using a starting alpha 
// returns -1.0 when an error has occurred.
double integrator::LineSearch(double start_alpha) {
    double min_alpha = 0.001;
    bool inversion = true;
    double alpha = start_alpha;

    Eigen::VectorXd dx;
    double residual;
    do {
        residual = TimeStep(dx);

        // Check for inversion
        if (residual == -2.0) {
            alpha = alpha/2.0;
            std::cout << "    - Step size resulted in inverted elements. Reducing step size by half: " << alpha << std::endl;
        }

        // Reduce alpha
        if (alpha <= min_alpha && inversion) {
            std::cout << "    - Step size reduced to minimum (<=" << min_alpha << ") without convergence. Terminating." << std::endl;
            break;
            return -1.0;
        }
    } while(inversion);
    // Completed one iteration, so we can add to positions
    if (!inversion) {
        return residual;
    }
    return -1.0;
}

// Main timestepping function
double integrator::TimeStep(Eigen::VectorXd& dv) {
    Eigen::SparseMatrix<double> _A;
    Eigen::VectorXd _b;
    Eigen::VectorXd _z;
    // Build system matrices
    buildSystem(_A, _b, _z);

    // Solve system using Eigen's Conjugate Gradient
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper> solver;
    solver.compute(_A);
    Eigen::VectorXd _y = solver.solve(_b);

    // Update velocities using the constraints
    dv = _y + _z;
    velocities_t += dv;
    // Update positions
    positions_t += h*velocities_t;

    // If we have a bad value, terminate the simulator
    if (!velocities_t.allFinite() || !positions_t.allFinite()) {
        return 0.0;
    }
    return 1.0;
}

// Assemble per-iteration z value
// Here we just do a simple velocity cancellation
void integrator::assembleZ(Eigen::VectorXd& z) {
    z.resize(_m.dim*_m.n);
    z.setZero();

    for (int v = 0; v < _m.n; v++) {
        Eigen::Vector3d vv;
        vv << velocities_t(_m.dim*v+0),
            velocities_t(_m.dim*v+1),
            velocities_t(_m.dim*v+2);

        Eigen::Matrix3d Sv = localS3D(constraints[v]);
        Eigen::Matrix3d Cv = Eigen::Matrix3d::Identity() - Sv;

        Eigen::Vector3d zv = -Cv * vv;

        z(_m.dim*v+0) = zv(0);
        z(_m.dim*v+1) = zv(1);
        z(_m.dim*v+2) = zv(2);
    }
    //z = IminusS * z;
    return;
}

// Check if a tet is inverted. If so, add it to a list
std::vector<int> integrator::checkElementInversion() {
    double threshold = 1e-2;
    Eigen::MatrixXd folded_pos = Utils::foldVector3d(positions_t);
    std::vector<int> inverted_tets;
    for (int e = 0; e < _m.num_t; e++) {
        std::vector<Utils::Vector3d> tet_pos;
        tet_pos.push_back(folded_pos.row(_m.t[e][0]));
        tet_pos.push_back(folded_pos.row(_m.t[e][1]));
        tet_pos.push_back(folded_pos.row(_m.t[e][2]));
        tet_pos.push_back(folded_pos.row(_m.t[e][3]));

        double volume = Utils::computeTetVolume(tet_pos);
        if (volume < threshold) {
            inverted_tets.push_back(e);
        }
    }
    return inverted_tets;
}

// Compute elastic force for an element (minus mass)
Eigen::MatrixXd integrator::elementElasticForce(double vol, Eigen::MatrixXd& dFdx, Utils::Matrix3d& PK1) {
    Eigen::MatrixXd localF = -1 * vol * ((dFdx).transpose() * Utils::vectorizeMatrix(PK1));   // dFdx is 9x12, vec(PK1) is 9x1
    return localF;
}

// Compute elastic force gradient (K) for an element
Eigen::MatrixXd integrator::elementElasticGradient(double vol, Utils::Matrix9d& Hessian, Eigen::MatrixXd& dFdx) {
    Eigen::MatrixXd localK = -1 * vol * (dFdx.transpose() * Hessian * dFdx); // dFdx is 9x12, Hessian is 9x9
    return localK;
}

// Assemble the velocity coffecient (LHS)
Eigen::SparseMatrix<double> integrator::assembleLHS(Eigen::SparseMatrix<double>& M, Eigen::SparseMatrix<double>& C, Eigen::SparseMatrix<double>& K) {
    return (M + h*C - h*h*K);
}

// Compute the force vector (RHS)
Eigen::VectorXd integrator::assembleRHS(Eigen::VectorXd& f_elast, Eigen::VectorXd& f_damp, Eigen::VectorXd& f_ext, Eigen::SparseMatrix<double>& K) {
    return (h*(f_elast + f_damp + f_ext) + (h*h * K)*velocities_t);
}

// Prefiltering step gets its own methods in case we ever need to change it
void integrator::preFilteringA(const Eigen::SparseMatrix<double>& A, Eigen::SparseMatrix<double>& A_global) {
    A_global = (S*A*S + IminusS);
    return;
}

void integrator::preFilteringb(const Eigen::SparseMatrix<double>& A, const Eigen::VectorXd& b, const Eigen::VectorXd& z, Eigen::VectorXd& b_global) {
    b_global = S * (b - A*z);
    return;
}

// Assemble system
// Note we're doing this and returning just the matrices so that I'm not hogging extra memory during the system solve
double integrator::buildSystem(Eigen::SparseMatrix<double>& A_global, Eigen::VectorXd& b_global, Eigen::VectorXd& z_global) {
    Eigen::SparseMatrix<double> K(_m.dim*_m.n, _m.dim*_m.n);  // Stiffness Matrix
    Eigen::SparseMatrix<double> C(_m.dim*_m.n, _m.dim*_m.n);  // Damping Matrix
    std::vector<T> tripletK;
    tripletK.reserve((_m.dim*_m.n));

    // Eigen::VectorXd f_ext(_m.dim * n);  // External forces are already pre-computed
    Eigen::VectorXd f_damp(_m.dim * _m.n); // Damping force
    Eigen::VectorXd f_elast(_m.dim * _m.n);    // Elastic force
    f_damp.setZero();
    f_elast.setZero();

    // Template Identity
    for (int v = 0; v < _m.n; v++) {
        for (int d = 0; d < _m.dim; d++) {
            tripletK.push_back(T(_m.dim*v + d, _m.dim*v + d, 1));
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

}   // namespace Solver