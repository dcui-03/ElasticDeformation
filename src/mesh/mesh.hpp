// mesh.hpp
#pragma once

#include "utils/EigenTypes.hpp"
#include <Eigen/Core>
#include <string>
#include <vector>

namespace Geom {

class mesh {
    public:
        mesh(std::vector<Eigen::Vector3d>& verts, std::vector<std::vector<int>>& connectivity);
        mesh(std::vector<Eigen::Vector3d>& verts, std::vector<std::vector<int>>& connectivity, double _rho);

        std::vector<Eigen::Vector3d>& v;     // ordered list of vertex positions
        std::vector<std::vector<int>>& t;    // List of vertex indices per element (tri/tet)
        std::vector<double> t_vols;         // volume of each element
        std::vector<Eigen::Matrix3d> D_mInv;   // List of rest deformation matrices
        double rho = 1.0;     // mass density
        int num_t;       // number of tets
        int n;   // number of vertices

        int verts_per_el = 4;   // Number of vertices in 1 element
        int dim; // Dimension of problem/number of ex. xyz components (2 for triangle mesh in 2D, 3 for thin shells in 3D, 3 for tets)
    protected:
        void initializeVariables();

        void setElementDim();

        virtual double computeElementVolume(std::vector<Eigen::Vector3d> el_verts);

        virtual Eigen::Matrix3d computeDmInv(std::vector<Eigen::Vector3d> el_verts);
    private:
        int el_dim; // Dimension of elements (2 for triangle mesh in 2D, 2 for thin shells in 3D, 3 for tets)
        std::string mesh_type = "mesh";
};

}   // namespace Geom