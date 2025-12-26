#include "polyscope/polyscope.h"
//#include "polyscope/surface_mesh.h"
#include "polyscope/volume_mesh.h"
#include "polyscope/point_cloud.h"

#include <Eigen/Core>
#include <chrono>
#include <iostream>
#include <string>
#include <cmath>
#include <tuple>
#include <array>
#include <vector>
#include <stdexcept>
#include <glm/glm.hpp>
#include <glm/vec3.hpp>

//#include "args/args.hxx"
#include "imgui.h"

// My files
#include "energy/arap.hpp"
#include "energy/snh.hpp"
#include "energy/energy.hpp"

#include "integrator/BDF1.hpp"
#include "integrator/BDF2.hpp"
#include "integrator/integrator.hpp"

#include "mesh/tetmesh.hpp"
#include "mesh/mesh.hpp"

#include "utils/IOutils.hpp"
#include "utils/matrixUtils.hpp"

// Main file for Deformation and visualization with Polyscope

// For timing simulation computation time
using Clock = std::chrono::high_resolution_clock;

bool loadingMode = false;

// Handles for the tet mesh and point cloud
polyscope::VolumeMesh* psMesh;
polyscope::PointCloud* psConstraintPC;

// VARIABLES FOR PARSING AND WRITING FILES
std::string InputPath;
std::string PlaybackPath;
int playbackStartIdx = 0;
int playbackEndIdx = 0;

// autoplaying
bool autoPlaying = false;
int currentSimFrame = 0;
int currentPSFrame = 0;
int framesPerUpdate = 6;

// scrubbing
bool scrubMode = false;
int scrubFrame = 0;

// TETMESH COPY
std::unique_ptr<Geom::mesh> m;
// TODO find a mu and lambda to initialize to
float PR = 0.4;   // PR is usually between 0 and 0.5 (0.5 may break solve)
double YM = 5e4;      // YM is usually between 1e4 and 1e6
float logE = std::log10(YM);

// Timestepping variables
int seconds = 10;
int fps = 5;
double spf = 1.0/fps;

// Simulator type
int sIdx = 0;   // 0 for BDF-1, 1 for BDF-2

// Energy and Forces
int eIdx = 0;   // 0 for SNH, 1 for ARAP
bool R = true;
bool G = true;
float Gx = 0;
float Gy = 0;
float Gz = -9.8;

// Constrained vertices
bool constraintMode = false;
std::vector<glm::vec3> psConstraints;
std::vector<std::vector<Eigen::Vector3d>> constraints;


// Compute Deformation given the current inputs
void solveDeformation() {
    std::cout << "====== HYPERELASTIC DEFORMATION ======\n" << std::endl;
    // Get time variables
    double t = 1.0/(fps);
    playbackEndIdx = seconds * fps;
    std::cout << "Total Frame Count: " << playbackEndIdx << std::endl;
    std::cout << "    - seconds: " << seconds << std::endl;
    std::cout << "    - t: " << t << "\n" << std::endl;

    // Declare energy
    std::unique_ptr<energy::Energy> e;
    double mu = Utils::solveMu(static_cast<double>(YM), static_cast<double>(PR));
    double lambda = Utils::solveLambda(static_cast<double>(YM), static_cast<double>(PR));
    
    if (eIdx == 0) {    // SNH
        std::cout << "ENERGY: Stable Neo-Hookean\n" << std::endl;
        e = std::make_unique<energy::SNH>(mu, lambda);
    } else {            // ARAP
        std::cout << "ENERGY: As-Rigid-As-Possible\n" << std::endl;
        e = std::make_unique<energy::ARAP>(mu, lambda);
    }
    std::cout << "Lamé Parameters (mu, lambda): (" << mu << ", " << lambda << ")" << std::endl;

    // Damping
    if (R) {
        std::cout << "DAMPING: Rayleigh Damping" << std::endl;
    } else {
        std::cout << "DAMPING: None" << std::endl;
    }
    // Decide on external force
    Eigen::Vector3d gravity = {static_cast<double>(Gx), static_cast<double>(Gy), static_cast<double>(Gz)}; // Change this for 2D
    if (G) {
        std::cout << "EXT. FORCE: USER-DEFINED" << std::endl;
        std::cout << "    (" << gravity(0) << ", " << gravity(1) << ", " << gravity(2) << ")\n" << std::endl;
    } else {
        std::cout << "EXT. FORCE: None\n" << std::endl;
        gravity(0) = 0.0; gravity(1) = 0.0; gravity(2) = 0.0;
    }

    // Constraints
    std::cout << "CONSTRAINTS: " << std::endl;
    for (int v = 0; v < constraints.size(); v++) {
        if (constraints[v].size() != 0) {
            std::cout << "    - Vertex " << v << ": " << "[x, y, z]" << std::endl;  // Locking all components for simplicity
        }
    }
    std::cout << std::endl;

    // Declare solver
    std::unique_ptr<Solver::integrator> timestepper;

    if (sIdx == 0) {
        // BDF1
        std::cout << "SOLVER: " << "BDF-1\n" << std::endl; // Change if using a different solver
        timestepper = std::make_unique<Solver::BDF1>(*m, *e, t, R, constraints, gravity);
    } else {
        // BDF2
        std::cout << "SOLVER: " << "BDF-2\n" << std::endl; // Change if using a different solver
        timestepper = std::make_unique<Solver::BDF2>(*m, *e, t, R, constraints, gravity);
    }

    std::vector<int> progress(2);
    progress[0] = playbackStartIdx;
    progress[1] = playbackEndIdx;

    // Save the starting state
    std::string framePath = PlaybackPath.substr(0, PlaybackPath.size()-4) + "_f0.txt";   // Update the saved file name
    IO::writeStateToTxt(framePath, timestepper->positions_t, timestepper->velocities_t, m->dim, progress);
    double totalTime = 0.0; // Timekeeper!

    // Iteratively solve each frame
    std::cout << "====== NOW SOLVING DEFORMATION ======" << std::endl; // Change if using a different solver
    for (int frame = 0; frame < playbackEndIdx; frame++) {
        auto t_start = Clock::now(); // Time starts NOW!
        progress[0]++;  // Update the frame count
        std::cout << "COMPUTING FRAME " << progress[0] << " OF " << progress[1] << std::endl;
        framePath = PlaybackPath.substr(0, PlaybackPath.size()-4) + "_f" + std::to_string(progress[0]) + ".txt";   // Update the saved file name
        // Solve the system
        int success = timestepper->TimeStep();
        // save positions and velocities
        IO::writeStateToTxt(framePath, timestepper->positions_t, timestepper->velocities_t, m->dim, progress);
        std::cout << "    Frame saved to <" << framePath << ">" << std::endl;
        auto t_end = Clock::now();
        std::chrono::duration<double> elapsed = t_end - t_start;
        std::cout << "    Frame computation time: " << elapsed.count() << " seconds\n" << std::endl;
        totalTime += elapsed.count();

        // If simulator failed or caught a bad value like a nan then terminate here
        if (success == 0) {
            std::cout << "SIMULATION TERMINATED at frame " << frame << " due to unexpected output (ex. nan)" << std::endl;
            std::cout << "See dumped frame file for more information" << std::endl;
            playbackEndIdx = frame;
            break;
        }
    }

    std::cout << "\n====== ALL FRAMES COMPUTED ======" << std::endl;
    // TODO maybe print some quantities like timing, etc.?
    std::cout << "Total Simulation Time: " << totalTime << " seconds" << std::endl;
    std::cout << "Avg. Simulation Time Per Frame: " << totalTime / (playbackEndIdx - playbackStartIdx) << " seconds\n" << std::endl;
    std::cout << "Press the PLAYBACK button to visualize deformation." << std::endl;
    return;
}

void playback(std::string& framePath) {
    // Note there are playbackEndIdx+1 frames since we include the rest frame
    // read this frame's positions and velocities
    std::vector<glm::vec3> V_glm;
    IO::readStateTxtToPS(framePath, V_glm, m->n);
    // Update positions in mesh
    psMesh->updateVertexPositions(V_glm);
}

// Add constraints of a vertex to the constraint vector
// So far, can only lock all dims in the major axis directions
// TODO Future: Enable locking in arbitrary direction individually (as long as 3 direcs are orthoNORMAL)
void addConstraintToVector(int v_idx) {    // Need to generalize this to all dims.
    // For now, just lock all components
    Eigen::Vector3d i0 = {1, 0, 0};
    Eigen::Vector3d i1 = {0, 1, 0};
    Eigen::Vector3d i2 = {0, 0, 1};

    constraints[v_idx].push_back(i0);
    constraints[v_idx].push_back(i1);
    constraints[v_idx].push_back(i2);
    return;
}

// A user-defined callback, for creating control panels (etc)
// Use ImGUI commands to build whatever you want here, see
// https://github.com/ocornut/imgui/blob/master/imgui.h
void myCallback() {
    ImGuiIO& io = ImGui::GetIO();
    
    bool mouseDown = ImGui::IsMouseDown(0);
    bool mouseClicked = ImGui::IsMouseClicked(0);
    bool mouseReleased = ImGui::IsMouseReleased(0);
    int selectedVertex = -1;
    glm::vec2 screen{io.MousePos.x, io.MousePos.y};

    polyscope::PickResult pick = polyscope::pickAtScreenCoords(screen);

    // Deformation stuff
    if (ImGui::Button("Compute Deformation")) {
        if (loadingMode) {
            std::cout << "Currently in loading mode. Run without flags to solve for a deformation." << std::endl;
        } else {
            constraintMode = false;
            scrubMode = false;
            autoPlaying = false;
            solveDeformation();
        }
    }

    if (ImGui::Button(autoPlaying ? "Stop Playback" : "Playback")) {
        scrubMode = false;  // Enforce this
        if (playbackEndIdx == 0) {
            std::cout << "No simulation initialized. Cannot perform playback" << std::endl;
        } else {
            autoPlaying = !autoPlaying;
            if (autoPlaying) {
                std::cout << "Playing Back Simulation." << std::endl;
                currentSimFrame = playbackStartIdx;
                // Get the nearest number of frames at 60 fps
                // I'm paranoid...
                framesPerUpdate = std::max(1, static_cast<int>(60/fps));
            } else {
                std::cout << "Stopping Playback." << std::endl;
            }
        }
    }

    // Updates the frame every specified number of frames
    if (autoPlaying) {
        constraintMode = false; // Enforce this to avoid confusion
        scrubMode = false;
        if (currentPSFrame % framesPerUpdate == 0) {    // One update every num frames
            std::string framePath = PlaybackPath.substr(0, PlaybackPath.size()-4) + "_f" + std::to_string(currentSimFrame) + ".txt";   // Update the saved file name
            playback(framePath);
            // Update simulation frame
            // If we run over the last index, then restart the simulation
            if (currentSimFrame == playbackEndIdx) {
                currentSimFrame = playbackStartIdx;
            } else {
                currentSimFrame++;
            }
        }
        // update polyscope frame
        currentPSFrame = (currentPSFrame+1)%framesPerUpdate;
    }

    ImGui::SameLine();
    if (ImGui::Button(scrubMode ? "Stop Scrubbing" : "Scrubbing")) {
        // Enforce no other modes allowed
        autoPlaying = false;
        constraintMode = false;
        if (playbackEndIdx == 0) {
            std::cout << "No simulation initialized. Cannot perform scrubbing" << std::endl;
        } else {
            scrubMode = !scrubMode;
            if (scrubMode) {
                std::cout << "Starting Scrub Mode." << std::endl;
                currentSimFrame = playbackStartIdx;
                // Get the nearest number of frames at 60 fps
                // I'm paranoid...
                framesPerUpdate = std::max(1, static_cast<int>(60/fps));
            } else {
                std::cout << "Stopping Scrub Mode" << std::endl;
            }
        }
    }

    if (ImGui::SliderInt("Scrub Frame", &scrubFrame, playbackStartIdx, playbackEndIdx) && scrubMode) {   // We can only scrub when not autoplaying
        // Allow frame selection
        std::string framePath = PlaybackPath.substr(0, PlaybackPath.size()-4) + "_f" + std::to_string(scrubFrame) + ".txt";   // Update the saved file name
        playback(framePath);
    }


    if (ImGui::Button(constraintMode ? "Stop Applying Constraints" : "Apply Constraint")) {
        if (loadingMode) {
            std::cout << "Currently in loading mode. Run without flags to solve for a deformation." << std::endl;
        } else if (autoPlaying) {
            std::cout << "Currently in Autoplaying mode. Stop Playback to add constraints." << std::endl;
        } else {
            constraintMode = !constraintMode;
        }
    }
    
    // Check if selected a mesh vertex
    if (constraintMode && mouseClicked && pick.isHit && pick.structure == psMesh) {
        polyscope::VolumeMeshPickResult meshPick = psMesh->interpretPickResult(pick);

        if (meshPick.elementType == polyscope::VolumeMeshElement::VERTEX) {
            selectedVertex = static_cast<int>(meshPick.index);
            // Add to constraint vector
            addConstraintToVector(selectedVertex);
            std::cout << "Added Vertex " << selectedVertex << " to constraint set." << std::endl;
            // Create point and add to point cloud
            Eigen::Vector3d v_xyz = m->v[selectedVertex];
            glm::vec3 temp_pnt;
            psConstraints.emplace_back(static_cast<float>(v_xyz(0)),
                                       static_cast<float>(v_xyz(1)),
                                       static_cast<float>(v_xyz(2)));
            psConstraintPC = polyscope::registerPointCloud("Constraints", psConstraints);
            psConstraintPC->setPointRadius(0.008);
        } else {
            std::cout << "Did not click on handle." << std::endl;
        }
    }

    ImGui::SameLine();
    // Reset the constraint set
    if (ImGui::Button("Reset Constraints")) {
        constraintMode = false;
        // Clear constraint vector
        for (int i = 0; i < constraints.size(); i++) {
            constraints[i].clear();
        }
        // Clear point cloud
        psConstraints.clear();
        std::cout << "Cleared all constraints." << std::endl;
        psConstraintPC = polyscope::registerPointCloud("Constraints", psConstraints);
    }

    // Simulation Time
    ImGui::SliderInt("Simulation Time (s)", &seconds, 1, 30);
    if (ImGui::SliderInt("Steps Per Second", &fps, 1, 100) && !autoPlaying) {   // We can only scrub when not autoplaying
        spf = 1.0/fps;
    }
    ImGui::Text("Timestep = %.3e", spf);

    // Simulator
    ImGui::RadioButton("BDF-1", &sIdx, 0);
    ImGui::SameLine();
    ImGui::RadioButton("BDF-2", &sIdx, 1);

    // Energies
    ImGui::RadioButton("SNH", &eIdx, 0);
    ImGui::SameLine();
    ImGui::RadioButton("ARAP", &eIdx, 1);

    // Material Parameters
    ImGui::SliderFloat("Poisson's Ratio", &PR, 0., 0.5, "%.2f");
    if (ImGui::SliderFloat("Young's Modulus (log10)", &logE, 4.0f, 6.0f, "%.2f")) {
        YM = std::pow(10.0, logE);
    }
    ImGui::Text("E = %.3e", YM);

    // Force Parameters
    ImGui::Checkbox("Rayleigh Damping", &R);
    ImGui::SameLine();
    ImGui::Checkbox("External Force", &G);

    ImGui::SliderFloat("Ext. Force (x)", &Gx, -10.0, 10.0);
    ImGui::SliderFloat("Ext. Force (y)", &Gy, -10.0, 10.0);
    ImGui::SliderFloat("Ext. Force (z)", &Gz, -10.0, 10.0);
}

int main(int argc, char **argv) {
    if (argc < 3) {
        std::cout << "Too few arguments. Usage: ./deformation <mesh file path> <simulation file path>" << std::endl;
        return 1;
    }
    InputPath = argv[1];
    PlaybackPath = argv[2]; // Sort of unsafe may need to add this to parse inputs
    // Parse inputs (add flags)
    for (int i = 3; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--p" && i + 2 < argc) {          // Playback  --p <start frame idx>  <number of frames to play>
            try {
                loadingMode = true;
                playbackStartIdx = std::stoi(argv[i+1]);   // This will get the file ex. savedfile_f3.txt
                scrubFrame = playbackStartIdx;
                playbackEndIdx = playbackStartIdx + std::stoi(argv[i+2]);  // This will go up to the file ex. savedfile_f7.txt
            } catch (const std::invalid_argument& e) {
                std::cerr << "Invalid input for --p.\n";
                return 1;
            }
        }
    }


    // Initialize polyscope
    polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::None; // Disable ground plane

    // Set the callback function
    polyscope::state::userCallback = myCallback;

    polyscope::init();

    // Set camera view
    polyscope::view::setUpDir(polyscope::UpDir::ZUp);      // Z up
    polyscope::view::setFrontDir(polyscope::FrontDir::NegYFront); // -Y forward

    // Set projection to orthographic
    polyscope::view::setProjectionMode(polyscope::ProjectionMode::Orthographic);

    // Load our tet mesh object
    std::vector<Eigen::Vector3d> V;
    std::vector<std::vector<int>> T;
    std::cout << "\nLoading TOBJ file" << std::endl;
    IO::loadTOBJ(InputPath, V, T, false);
    std::cout << "Loading TOBJ data to internal tetmesh object" << std::endl;
    m = std::make_unique<Geom::tetmesh>(V, T);
    // Initialize constraints
    constraints.resize(m->n);

    std::vector<glm::vec3> psV;
    std::vector<std::array<size_t, 4>> psT;
    std::cout << "\nConverting to Polyscope TetMesh" << std::endl;
    IO::polyscopeTetConverter(V, T, psV, psT);

    // Register tet mesh with PS
    std::cout << "Registering Tet Mesh to Polyscope" << std::endl;
    psMesh = polyscope::registerTetMesh("Simulation Mesh", psV, psT);

    // Empty point cloud for constraints
    psConstraintPC = polyscope::registerPointCloud("Constraints", psConstraints);

    std::cout << "\nTotal Vertex Count: " << m->n << std::endl;
    std::cout << "Total Tet Count: " << m->num_t << "\n" << std::endl;

    // Give control to the polyscope gui
    polyscope::show();

    return EXIT_SUCCESS;
}