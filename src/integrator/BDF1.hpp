// BDF1.hpp
#pragma once

#include "integrator.hpp"
#include "energy/energy.hpp"
#include "energy/arap.hpp"
#include "energy/snh.hpp"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <string>

// A class for BDF1
namespace Solver {

class BDF1 : public integrator {
    public:
        using integrator::integrator;
        int TimeStep() override;
    protected:
        // The integrator class already defines most of our functions for us, so BDF1 just needs to
        // call the appropriate functions and assemble the system.

        // Since prefiltering is standardized, we also let the parent function do that
        
        // Assembling system parts
        Eigen::SparseMatrix<double> assembleLHS(Eigen::SparseMatrix<double>& M, Eigen::SparseMatrix<double>& C, Eigen::SparseMatrix<double>& K) override;

        Eigen::VectorXd assembleRHS(Eigen::VectorXd& f_elast, Eigen::VectorXd& f_damp, Eigen::VectorXd& f_ext, Eigen::SparseMatrix<double>& K) override;

        // Assemble all system parts
        void buildSystem(Eigen::SparseMatrix<double>& A_global, Eigen::VectorXd& b_global, Eigen::VectorXd& z_global) override;
    
    private:
        std::string integrator_type = "BDF-1";
};

} // namespace Solver