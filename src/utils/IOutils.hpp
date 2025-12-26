// Some functions to help with loading files for playback
#include <Eigen/Dense>

#include "EigenTypes.hpp"
#include <glm/vec3.hpp>
#include <array>
#include <string>
#include <vector>

namespace IO {

// Function to load TOBJ files (as specified in DD and HOBAK)
bool loadTOBJ(
    const std::string& path,
    std::vector<Utils::Vector3d>& outVertices,
    std::vector<std::vector<int>>& outTets,
    bool indicesAreOneBased = false
);

// Function to convert our tobj loader output to a format polyscope can read
void polyscopeTetConverter(
    const std::vector<Utils::Vector3d>& V_eigen,
    const std::vector<std::vector<int>>& T_in,
    std::vector<glm::vec3>& V,
    std::vector<std::array<size_t, 4>>& T
);


// writes current state to a txt file
void writeStateToTxt(
    const std::string& filePath,
    const Eigen::VectorXd& p,
    const Eigen::VectorXd& v,
    const double dim,
    std::vector<int> progress
);

// Reads the txt file and extracts the vertex positions in a format we can input to polyscope
void readStateTxtToPS(
    const std::string& filePath,
    std::vector<glm::vec3>& V_glm,
    int expectedVertexCount
);

} // namespace IO