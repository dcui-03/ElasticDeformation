#include "mesh.hpp"

#include "utils/EigenTypes.hpp"
#include "polyscope/polyscope.h"
#include "glm/vec3.hpp"
#include <Eigen/Dense>
#include <iostream>
#include <array>
#include <vector>

namespace Geom {
// Constructors
mesh::mesh(std::vector<Eigen::Vector3d>& verts, std::vector<std::vector<int>>& connectivity): v(verts), t(connectivity) {
    initializeVariables();
    setElementDim();
}

mesh::mesh(std::vector<Eigen::Vector3d>& verts, std::vector<std::vector<int>>& connectivity, double _rho): v(verts), t(connectivity), rho(_rho) {
    initializeVariables();
    setElementDim();
}

// Initialization
void mesh::initializeVariables() {
    num_t = t.size();
    n = v.size();
    verts_per_el = t[0].size();
    dim = v[0].size();
    // Compute element volumes and D_m
    t_vols.resize(num_t);
    D_mInv.resize(num_t);
    for (int e = 0; e < num_t; e++) {
        std::vector<Eigen::Vector3d> vert_pos(verts_per_el);
        for (int vert = 0; vert < verts_per_el; vert++) {
            vert_pos[vert] = v[t[e][vert]];
        }
        double temp_vol = computeElementVolume(vert_pos);
        // If any volume is negative swap any two vertices so that volume is positive
        if (temp_vol < 0.0) {
            int p0 = t[e][0]; int p1 = t[e][1];
            t[e][0] = p1; t[e][1] = p0;
            Eigen::Vector3d vert_temp = vert_pos[0];
            vert_pos[0] = vert_pos[1]; vert_pos[1] = vert_temp;
            temp_vol = -1 * temp_vol;
        }
        t_vols[e] = temp_vol;
        D_mInv[e] = computeDmInv(vert_pos);
    }
    return;
}

void mesh::setElementDim() {
    el_dim = 3; // Change this if needed
    return;
}

// Placeholder
double mesh::computeElementVolume(std::vector<Eigen::Vector3d> el_verts) {
    return 1;
}

// Placeholder
Eigen::Matrix3d mesh::computeDmInv(std::vector<Eigen::Vector3d> el_verts) {
    Eigen::MatrixXd temp(3, 3);
    temp.setIdentity();
    return temp;
}

} // namespace Geom