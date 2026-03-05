// Static2.hpp
#pragma once

#include "integrator.hpp"
#include "energy/energy.hpp"
#include "energy/arap.hpp"
#include "energy/snh.hpp"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <string>

// A class for Static solver (i.e., static equilibrium)
namespace Solver {

class Static : public integrator {
    public:
        Static(Geom::mesh& m, Geom::mesh& m_start, energy::Energy& e, std::vector<std::vector<Utils::Vector3d>>& c, Utils::Vector3d& a_ext);
        
        double LineSearch(double start_alpha) override;
        // Take a single Newton step in solver and return residual
        double TimeStep(Eigen::VectorXd& dx) override;
    protected:
        Eigen::SparseMatrix<double> assembleLHS(Eigen::SparseMatrix<double>& K);

        Eigen::VectorXd assembleRHS(Eigen::VectorXd& f_ext, Eigen::VectorXd& f_int);

        // Assemble all system parts and return Stiffness
        double buildSystem(Eigen::SparseMatrix<double>& A_global, Eigen::VectorXd& b_global, Eigen::VectorXd& z_global) override;
    private:
        //Eigen::VectorXd positions_0;
        std::string integrator_type = "Static";
};

} // namespace Solver