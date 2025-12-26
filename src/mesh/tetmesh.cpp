#include "tetmesh.hpp"

#include "utils/EigenTypes.hpp"
#include <Eigen/Dense>
#include <iostream>
#include <cmath>

namespace Geom {

tetmesh::tetmesh(std::vector<Eigen::Vector3d>& verts, std::vector<std::vector<int>>& connectivity): mesh(verts, connectivity) {
    initializeVariables();
    setElementDim();
}
    
/*void tetmesh::initializeVariables() {
    num_t = t.size();
    n = v.size();
    verts_per_el = t[0].size();
    // Compute element volumes and D_m
    t_vols.resize(num_t);
    D_m.resize(num_t);
    std::cout << "Looping over tets and computing D_m and volumes" << std::endl;
    for (int e = 0; e < num_t; e++) {
        std::vector<Utils::Vector3d> vert_pos(verts_per_el);
        for (int vert = 0; vert < verts_per_el; vert++) {
            vert_pos[vert] = v[t[e][vert]];
        }
        std::cout << "Element volume" << std::endl;
        t_vols[e] = computeElementVolume(vert_pos);
        std::cout << "Element D_m" << std::endl;
        D_m[e] = computeDm(vert_pos);
    }
    return;
}*/

// Stole a lot of this function's structure from HOBAK (TET_MESH.cpp)
double tetmesh::computeElementVolume(std::vector<Utils::Vector3d> tet_verts) {
    Utils::Vector3d diff1 = tet_verts[1] - tet_verts[0];
    Utils::Vector3d diff2 = tet_verts[2] - tet_verts[0];
    Utils::Vector3d diff3 = tet_verts[3] - tet_verts[0];
    double vol = ((diff3.dot((diff1).cross(diff2)) / 6.0));
    return vol;
}

// Compute D_m, which we use to compute the deformation gradient and dF/dx
Utils::Matrix3d tetmesh::computeDmInv(std::vector<Utils::Vector3d> tet_old) {
    Utils::Matrix3d Dm;
    Dm.setZero();
    Dm.col(0) = tet_old[1] - tet_old[0];
    Dm.col(1) = tet_old[2] - tet_old[0];
    Dm.col(2) = tet_old[3] - tet_old[0];
    return Dm.inverse();
}

}   // namespace Geom