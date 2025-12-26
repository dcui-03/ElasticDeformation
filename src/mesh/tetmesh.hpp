// tetmesh.hpp
#pragma once

#include "mesh.hpp"
#include "utils/EigenTypes.hpp"
#include <Eigen/Core>
#include <string>
#include <vector>

namespace Geom {

class tetmesh : public mesh {
    public:
        tetmesh(std::vector<Eigen::Vector3d>& verts, std::vector<std::vector<int>>& connectivity);
    protected:
        double computeElementVolume(std::vector<Utils::Vector3d> el_verts) override;

        Utils::Matrix3d computeDmInv(std::vector<Utils::Vector3d> el_verts) override;
    private:
        int el_dim = 3;
        std::string mesh_type = "tetmesh";
};

}   // namespace Geom