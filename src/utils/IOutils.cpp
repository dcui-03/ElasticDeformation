// Some functions to help with loading files for playback
#include "IOutils.hpp"

#include "EigenTypes.hpp"
#include <Eigen/Dense>
#include <glm/vec3.hpp>
#include <array>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <stdexcept>

namespace IO {

// GLM::vec3 to Eigen::Vector3d converter
Eigen::Vector3d glmToEigen(glm::vec3 input) {
    Eigen::Vector3d output;
    output(0) = static_cast<double>(input.x);
    output(1) = static_cast<double>(input.y);
    output(2) = static_cast<double>(input.z);
    return output;
}

// Eigen::Vector3d to GLM::vec3 converter
glm::vec3 eigenToGLM(Eigen::Vector3d input) {
    glm::vec3 output;
    output.x = static_cast<float>(input(0));
    output.y = static_cast<float>(input(1));
    output.z = static_cast<float>(input(2));
    return output;
}

// Function to load TOBJ files (as specified in DD and HOBAK)
bool loadTOBJ(
    const std::string& path,
    std::vector<Utils::Vector3d>& outVertices,
    std::vector<std::vector<int>>& outTets,
    bool indicesAreOneBased
) {
    std::ifstream in(path);
    if (!in.is_open()) {
        throw std::runtime_error("loadTOBJ: failed to open file: " + path);
    }
    outVertices.clear();
    outTets.clear();

    std::string line;
    size_t lineNo = 0;

    while (std::getline(in, line)) {
        ++lineNo;

        // Strip comments (anything after '#')
        if (auto pos = line.find('#'); pos != std::string::npos) {
            line = line.substr(0, pos);
        }

        // Skip empty / whitespace-only lines
        auto firstNonWs = line.find_first_not_of(" \t\r\n");
        if (firstNonWs == std::string::npos) continue;

        std::istringstream ss(line);
        std::string tag;
        ss >> tag;

        if (tag == "v") {
            double x, y, z;
            if (!(ss >> x >> y >> z)) {
                throw std::runtime_error(
                    "loadTOBJ: malformed vertex line at " + std::to_string(lineNo) + ": " + line
                );
            }
            outVertices.emplace_back(x, y, z);
        }
        else if (tag == "t") {
            int a, b, c, d;
            if (!(ss >> a >> b >> c >> d)) {
                throw std::runtime_error(
                    "loadTOBJ: malformed tet line at " + std::to_string(lineNo) + ": " + line
                );
            }

            if (indicesAreOneBased) { a--; b--; c--; d--; }

            auto inRange = [&](int idx) {
                return idx >= 0 && idx < static_cast<int>(outVertices.size());
            };
            // Optional: you can defer this check until after reading the whole file,
            // but doing it here catches errors early if tets reference future vertices.
            if (!inRange(a) || !inRange(b) || !inRange(c) || !inRange(d)) {
                throw std::runtime_error(
                    "loadTOBJ: tet index out of range at line " + std::to_string(lineNo) +
                    " (vertices read so far: " + std::to_string(outVertices.size()) + ")"
                );
            }

            outTets.push_back({a, b, c, d});
        }
        else {
            // Ignore unknown tags (or change to throw if you want strict parsing)
            throw std::runtime_error("loadTOBJ: unknown tag at line " + std::to_string(lineNo));
            continue;
        }

        // Extra tokens on a line are ignored; change if you want strictness
    }
    
    return true;
}

// Function to convert our tobj loader output to a format polyscope can read
void polyscopeTetConverter(
    const std::vector<Utils::Vector3d>& V_eigen,
    const std::vector<std::vector<int>>& T_in,
    std::vector<glm::vec3>& V,
    std::vector<std::array<size_t, 4>>& T
) {
    V.clear();
    T.clear();
    // 1) Convert vertex positions to glm::vec3
    V.reserve(V_eigen.size());
    for (const auto& p : V_eigen) {
        V.emplace_back(
            static_cast<float>(p.x()),
            static_cast<float>(p.y()),
            static_cast<float>(p.z())
        );
    }

    // 2) Convert tets to array<4> of size_t (recommended)
    T.reserve(T_in.size());

    const size_t nV = V_eigen.size();
    for (const auto& tet : T_in) {
        if (tet.size() != 4) {
            throw std::runtime_error("polyscopeTetConverter: tet does not have 4 indices.");
        }
        // range check + cast
        for (int k = 0; k < 4; ++k) {
            if (tet[k] < 0 || static_cast<size_t>(tet[k]) >= nV) {
                throw std::runtime_error("polyscopeTetConverter: tet index out of range.");
            }
        }
        T.push_back({{
            static_cast<size_t>(tet[0]),
            static_cast<size_t>(tet[1]),
            static_cast<size_t>(tet[2]),
            static_cast<size_t>(tet[3])
        }});
    }

    return;
}

void convertPolyscopeToEigen(
    const std::vector<glm::vec3>& ps_verts,
    const std::vector<std::array<size_t, 4>>& ps_connectivity,
    std::vector<Eigen::Vector3d>& verts,
    std::vector<std::vector<int>>& tets)
{
    // ---------- Copy vertex positions ----------
    verts.resize(ps_verts.size());

    for (size_t i = 0; i < ps_verts.size(); ++i) {
        const glm::vec3& p = ps_verts[i];
        verts[i] = Eigen::Vector3d(
            static_cast<double>(p.x),
            static_cast<double>(p.y),
            static_cast<double>(p.z)
        );
    }

    // ---------- Copy connectivity ----------
    tets.resize(ps_connectivity.size());

    for (size_t i = 0; i < ps_connectivity.size(); ++i) {

        const std::array<size_t, 4>& tet = ps_connectivity[i];

        tets[i].resize(4);

        for (int j = 0; j < 4; ++j) {

            // Optional safety check
            if (tet[j] >= ps_verts.size()) {
                throw std::runtime_error(
                    "convertPolyscopeTetToEigen(): connectivity index out of bounds."
                );
            }

            tets[i][j] = static_cast<int>(tet[j]);
        }
    }
}

// writes current state to a txt file
void writeStateToTxt(
    const std::string& filePath,
    const Eigen::VectorXd& p,
    const Eigen::VectorXd& v,
    const double dim,
    std::vector<int> progress
) {
    double eps = 1e-5;
    std::ofstream out(filePath);
    if (!out.is_open()) {
        throw std::runtime_error("writeStateToTxt: failed to open file: " + filePath);
    }
    // Comment a line referencing the frame number
    out << "# Frame " << progress[0] << " of " << progress[1] << "\n";
    const int n = static_cast<int>(p.size());
    for (int i = 0; i < n; i += dim) {
        for (int d = 0; d < dim; d++) {
            out << std::to_string(p[i+d])     << " ";
        }
        for (int d = 0; d < dim; d++) {
            double temp_v = v[i+d];
            if (temp_v <= eps && temp_v >= -eps) {
                temp_v = 0.0;
            }
            out << std::to_string(v[i+d])     << " ";
        }
        out << "\n";
    }
    return;
}

// Reads the txt file and extracts the vertex positions in a format we can input to polyscope
void readStateTxtToPS(
    const std::string& filePath,
    std::vector<glm::vec3>& V_glm,
    int expectedVertexCount
) {
    std::ifstream in(filePath);
    if (!in.is_open()) {
        throw std::runtime_error("readPVTxtToGLM: failed to open file: " + filePath);
    }

    V_glm.clear();
    V_glm.reserve(expectedVertexCount);

    std::string line;
    size_t lineNo = 0;

    while (std::getline(in, line)) {
        ++lineNo;

        // Skip blank or comment lines
        auto firstNonWs = line.find_first_not_of(" \t\r\n");
        if (firstNonWs == std::string::npos) continue;
        if (line[firstNonWs] == '#') continue;

        if (V_glm.size() >= expectedVertexCount) {
            throw std::runtime_error(
                "readPVTxtToGLM: file has more rows than expected vertices "
                "(extra data at line " + std::to_string(lineNo) + ")."
            );
        }

        std::istringstream ss(line);

        double px, py, pz;
        double vx, vy, vz; // read but ignored

        if (!(ss >> px >> py >> pz >> vx >> vy >> vz)) {
            throw std::runtime_error(
                "readPVTxtToGLM: malformed line " + std::to_string(lineNo) +
                " (expected 6 numbers): " + line
            );
        }

        V_glm.emplace_back(
            static_cast<float>(px),
            static_cast<float>(py),
            static_cast<float>(pz)
        );
    }

    return;
}

} // namespace IO