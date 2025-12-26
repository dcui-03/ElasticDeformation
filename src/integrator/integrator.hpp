// integrator.hpp
#pragma once

#include "utils/EigenTypes.hpp"
#include "mesh/mesh.hpp"
#include "mesh/tetmesh.hpp"
#include "energy/energy.hpp"
#include "energy/arap.hpp"
#include "energy/snh.hpp"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <string>

namespace Solver {

class integrator {
    public:
        typedef Eigen::Triplet<double> T;

        integrator(Geom::mesh& m, energy::Energy& e, double t, bool R, std::vector<std::vector<Utils::Vector3d>>& c, Utils::Vector3d& a_ext);

        // Main function for timestepping procedure
        virtual int TimeStep();

        // Current State
        Eigen::VectorXd positions_t;
        Eigen::VectorXd velocities_t;
    protected:
        // Function for initializing and pre-computing
        virtual void initializePV();

        virtual void initializeDamping();

        virtual void initializeMass();

        virtual void initializeExtForce(Utils::Vector3d& a_ext);

        virtual Eigen::Matrix3d localS3D(const std::vector<Utils::Vector3d> c);

        virtual void assemblePrefiltering();

        virtual void assembleZ(Eigen::VectorXd& z);

        // Force Computations
        virtual Eigen::MatrixXd elementElasticForce(double vol, Eigen::MatrixXd& dFdx, Utils::Matrix3d& PK1);

        virtual Eigen::MatrixXd elementElasticGradient(double vol, Utils::Matrix9d& Hessian, Eigen::MatrixXd& dFdx);

        // Prefiltering functions
        virtual void preFilteringA(const Eigen::SparseMatrix<double>& A, Eigen::SparseMatrix<double>& A_global);
        virtual void preFilteringb(const Eigen::SparseMatrix<double>& A, const Eigen::VectorXd& b, const Eigen::VectorXd& z, Eigen::VectorXd& b_global);

        // Assembling system parts
        virtual Eigen::SparseMatrix<double> assembleLHS(Eigen::SparseMatrix<double>& M, Eigen::SparseMatrix<double>& C, Eigen::SparseMatrix<double>& K);

        virtual Eigen::VectorXd assembleRHS(Eigen::VectorXd& f_elast, Eigen::VectorXd& f_damp, Eigen::VectorXd& f_ext, Eigen::SparseMatrix<double>& K);

        // Assemble all system parts
        virtual void buildSystem(Eigen::SparseMatrix<double>& A_global, Eigen::VectorXd& b_global, Eigen::VectorXd& z_global);


        // VARIABLES
        Geom::mesh& _m;
        energy::Energy& _e;

        double h;   // Timestep size

        // Damping Force
        bool Rayleigh = true;
        double alpha = 0.01;   // For Rayleigh damping
        double beta = 0.01;    // For Rayleigh damping

        // External Force
        bool ExternalForce = true;
        //Eigen::VectorXd extAcceleration;    // Direction of force
        Eigen::VectorXd f_ext;  // Actual vector containing forces on vertices

        // A mapping from vertices to matrix rows (if we need it)
        std::map<int, int> VertToRow();

        // Vector for pinned components
        // i.e., initialize to 0. If we pin the x, y, or z component, 
        std::vector<std::vector<Utils::Vector3d>>& constraints;
        Eigen::SparseMatrix<double> M;
        // Matrix for prefiltering
        Eigen::SparseMatrix<double> S;
        Eigen::SparseMatrix<double> IminusS;
    
    private:
        std::string integrator_type = "Template";
};

} // namespace Solver