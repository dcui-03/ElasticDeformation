#include "polyscope/polyscope.h"
//#include "polyscope/surface_mesh.h"
#include "polyscope/volume_mesh.h"
#include "polyscope/point_cloud.h"

#include <Eigen/Core>
#include <chrono>
#include <iostream>
#include <string>
#include <cmath>
#include <map>
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

#include "integrator/Static.hpp"
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
std::vector<glm::vec3> psV;
// Eigen::MatrixXd psV
std::vector<std::array<size_t, 4>> psT;
// Eigen::MatrixXi TT
polyscope::VolumeMesh* psMesh;
polyscope::PointCloud* psConstraintPC;
bool showPC;

// VARIABLES FOR PARSING AND WRITING FILES
std::string InputPath;
std::string PlaybackPath;
int playbackStartIdx = 0;
int playbackEndIdx = 0;

// move vertices
int selectedVertex = -1;
bool gizmoMode = false; // This allows the user to create a gizmo
bool activeGizmo = false; // This tells us if there is an active gizmo
static polyscope::TransformationGizmo* vertexGizmo = nullptr;

// autoplaying
bool autoPlaying = false;
int currentSimFrame = 0;
int currentPSFrame = 0;
int framesPerUpdate = 6;

// scrubbing
bool scrubMode = false;
int scrubFrame = 0;

// TETMESH COPY
std::vector<Eigen::Vector3d> V;
std::vector<std::vector<int>> T;
std::vector<Eigen::Vector3d> V_current;
std::vector<std::vector<int>> T_current;
std::unique_ptr<Geom::mesh> m;
// TODO find a mu and lambda to initialize to
float PR = 0.35;   // PR is usually between 0 and 0.5 (0.5 may break solve)
double YM = 5e4;      // YM is usually between 1e4 and 1e6
float logE = std::log10(YM);

// Timestepping variables
int seconds = 10;
int fps = 5;
double spf = 1.0/fps;

// Simulator type
int sIdx = 2;   // 0 for Static, 1 for BDF-1, 2 for BDF-2

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
std::map<int, int> VtoConstraint;   // Maps polyscope mesh vertex to constraint vertex
// Constraint directions
bool customConstraint = false;
float Cx = 1.0;
float Cy = 1.0;
float Cz = 1.0;
// Canonical Frame axes (remove these when implementing arbitrary constraints)
Eigen::Vector3d i0 = {1, 0, 0};
Eigen::Vector3d i1 = {0, 1, 0};
Eigen::Vector3d i2 = {0, 0, 1};
std::vector<Eigen::Vector3d> constraintDirecs = {i0, i1, i2};
bool xAxis = true;
bool yAxis = true;
bool zAxis = true;


// Compute Deformation given the current inputs
void solveDeformation() {
    std::cout << "====== HYPERELASTIC DEFORMATION ======\n" << std::endl;
    // Copy current mesh state from Polyscope
    std::unique_ptr<Geom::mesh> m_start = std::make_unique<Geom::tetmesh>(V_current, T_current);

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
    bool abs = (PR >= 0.49);

    
    if (eIdx == 0) {    // SNH
        std::cout << "ENERGY: Stable Neo-Hookean\n" << std::endl;
        e = std::make_unique<energy::SNH>(mu, lambda, abs);
    } else {            // ARAP
        std::cout << "ENERGY: As-Rigid-As-Possible\n" << std::endl;
        e = std::make_unique<energy::ARAP>(mu, lambda, abs);
    }
    std::cout << "Lamé Parameters (mu, lambda): (" << mu << ", " << lambda << ")" << std::endl;
    if (abs) {
        std::cout << "Material is highly incompressible (Poisson's Ratio -> 0.5). For stability, PSD Hessian computed using ABS rule." << std::endl;
    } else {
        std::cout << "Material is compressible. PSD Hessian computed using CLAMP rule." << std::endl;
    }

    // Damping
    if (R && sIdx != 0) {
        std::cout << "DAMPING: Rayleigh Damping" << std::endl;
    } else {
        std::cout << "DAMPING: None" << std::endl;
        R = false;
    }
    // Decide on external force
    Eigen::Vector3d gravity = {static_cast<double>(Gx), static_cast<double>(Gy), static_cast<double>(Gz)}; // Change this for 2D
    if (G && sIdx != 0) {
        std::cout << "EXT. FORCE: USER-DEFINED" << std::endl;
        std::cout << "    (" << gravity(0) << ", " << gravity(1) << ", " << gravity(2) << ")\n" << std::endl;
    } else {
        std::cout << "EXT. FORCE: None\n" << std::endl;
        gravity(0) = 0.0; gravity(1) = 0.0; gravity(2) = 0.0;
        G = false;
    }

    // Constraints
    std::cout << "CONSTRAINTS: " << std::endl;
    for (int v = 0; v < constraints.size(); v++) {
        if (constraints[v].size() != 0) {
            std::cout << "    - Vertex " << v << ": ";
            for (int c = 0; c < constraints[v].size(); c++) {
                std::cout << "[" << constraints[v][c][0] << ", " <<
                                    constraints[v][c][1] << ", " << 
                                    constraints[v][c][2] << "]";
                if (c != constraints.size()-1) {
                    std::cout << ";  ";
                }
            }
            std::cout << std::endl;
        }
    }
    std::cout << std::endl;

    // Number of frames
    std::vector<int> progress(2);
    progress[0] = playbackStartIdx;
    progress[1] = playbackEndIdx;
    int num_frames = progress[1] - progress[0];
    
    // Declare solver
    std::unique_ptr<Solver::integrator> timestepper;

    if (sIdx == 0) {
        num_frames = 0;
        // Static solver
        std::cout << "SOLVER: " << "Static\n" << std::endl; // Change if using a different solver
        timestepper = std::make_unique<Solver::Static>(*m, *m_start, *e, constraints, gravity);
    } else if (sIdx == 1) {
        // BDF1
        std::cout << "SOLVER: " << "BDF-1\n" << std::endl; // Change if using a different solver
        timestepper = std::make_unique<Solver::BDF1>(*m, *m_start, *e, t, R, constraints, gravity);
    } else if (sIdx == 2) {
        // BDF2
        std::cout << "SOLVER: " << "BDF-2\n" << std::endl; // Change if using a different solver
        // NOTE: The second input should be the previous state for BDF2. Here we assume it was the same as the current start
        timestepper = std::make_unique<Solver::BDF2>(*m, *m_start, *m_start, *e, t, R, constraints, gravity);
    }

    // Save the starting state
    std::string framePath = PlaybackPath.substr(0, PlaybackPath.size()-4) + "_f0.txt";   // Update the saved file name
    IO::writeStateToTxt(framePath, timestepper->positions_t, timestepper->velocities_t, m->dim, progress);
    double totalTime = 0.0; // Timekeeper!

    Eigen::VectorXd delta;   // Dx or dv

    // Iteratively solve each frame
    std::cout << "====== NOW SOLVING DEFORMATION ======" << std::endl;
    if (sIdx != 0) {    // Time integration
        for (int frame = 0; frame < playbackEndIdx; frame++) {
            auto t_start = Clock::now(); // Time starts NOW!
            progress[0]++;  // Update the frame count
            std::cout << "COMPUTING FRAME " << progress[0] << " OF " << progress[1] << std::endl;
            framePath = PlaybackPath.substr(0, PlaybackPath.size()-4) + "_f" + std::to_string(progress[0]) + ".txt";   // Update the saved file name
            // Solve the system
            int success = timestepper->TimeStep(delta);
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
    } else {    // Quasi-static
        double residual;
        int max_frames = 5e2;
        double min_residual = 1e-5;
        double alpha = 1.0;
        std::cout << "Solving static equilibrium via Newton iterations." << std::endl;
        do {
            auto t_start = Clock::now(); // Time starts NOW!
            std::cout << "COMPUTING ITERATION " << num_frames + 1 << std::endl;
            residual = timestepper->LineSearch(alpha);
            if (residual < 0.0) {
                std::cout << "\nSOLVER TERMINATED at iteration " << num_frames + 1 << " due to unexpected output (ex. nan, inverted element)" << std::endl;
                break;
            }
            std::cout << "    Residual: " << residual << "\n" << std::endl;

            auto t_end = Clock::now();
            std::chrono::duration<double> elapsed = t_end - t_start;
            std::cout << "    Iteration computation time: " << elapsed.count() << " seconds\n" << std::endl;
            totalTime += elapsed.count();
            num_frames++;

            // save positions and velocities
            framePath = PlaybackPath.substr(0, PlaybackPath.size()-4) + "_f" + std::to_string(num_frames) + ".txt";   // Update the saved file name
            IO::writeStateToTxt(framePath, timestepper->positions_t, timestepper->velocities_t, m->dim, progress);
            std::cout << "    Frame saved to <" << framePath << ">" << std::endl;
        } while (residual > min_residual && num_frames < max_frames);
        if (residual <= min_residual) {
            std::cout << "\nStatic Solver Converged in " << num_frames << " iterations." << std::endl;
        } else if (num_frames >= max_frames) {
            std::cout << "\nExceeded iteration count (" << max_frames << "). Stopping." << std::endl;
        }
        // Change UI parameter so we can scrub all frames
        playbackEndIdx = std::min(num_frames, max_frames);
    }

    std::cout << "\n====== ALL FRAMES COMPUTED ======" << std::endl;
    // TODO maybe print some quantities like timing, etc.?
    std::cout << "Total Simulation Time: " << totalTime << " seconds" << std::endl;
    std::cout << "Avg. Simulation Time Per Frame: " << totalTime / (std::max(1, num_frames)) << " seconds\n" << std::endl;
    std::cout << "Press the PLAYBACK button to visualize deformation." << std::endl;
    return;
}

// Display the constraint vectors that the user applied to each constrained vertex
void displayConstraintDirecs() {
    // Iterate over map and build 3 vector lists
    std::vector<Eigen::Vector3d> constraints0(psConstraints.size());
    std::vector<Eigen::Vector3d> constraints1(psConstraints.size());
    std::vector<Eigen::Vector3d> constraints2(psConstraints.size());
    for (const auto& cVert : VtoConstraint) {
        constraints0[VtoConstraint[cVert.first]] = constraints[cVert.first][0];
        if (constraints[cVert.first].size() >= 2) {
            constraints1[VtoConstraint[cVert.first]] = constraints[cVert.first][1];
        } else {
            constraints1[VtoConstraint[cVert.first]] = Eigen::Vector3d::Zero();
        }
        if (constraints[cVert.first].size() >= 3) {
            constraints2[VtoConstraint[cVert.first]] = constraints[cVert.first][2];
        } else {
            constraints2[VtoConstraint[cVert.first]] = Eigen::Vector3d::Zero();
        }
    }
    auto* c0_handle = psConstraintPC->addVectorQuantity("c0", constraints0);
    c0_handle->setEnabled(true);
    c0_handle->setVectorLengthScale(0.04);
    // Also show second and third (filled with 0's)
    auto* c1_handle = psConstraintPC->addVectorQuantity("c1", constraints1);
    c1_handle->setEnabled(true);
    c1_handle->setVectorLengthScale(0.04);
    auto* c2_handle = psConstraintPC->addVectorQuantity("c2", constraints2);
    c2_handle->setEnabled(true);
    c2_handle->setVectorLengthScale(0.04);
}

// Playback a mesh state by loading a file
void playback(std::string& framePath) {
    // Note there are playbackEndIdx+1 frames since we include the rest frame
    // read this frame's positions and velocities
    IO::readStateTxtToPS(framePath, psV, m->n);
    // Update positions in mesh
    psMesh->updateVertexPositions(psV);
    // Update positions in point cloud
    for (const auto& pair : VtoConstraint) {
        psConstraints[pair.second] = psV[pair.first];
    }
    psConstraintPC = polyscope::registerPointCloud("Constraints", psConstraints);
    // Update constraint vectors
    displayConstraintDirecs();
}

// Add a gizmo at the specified vertex idx
void beginVertexEdit(int vert_idx) {
    glm::vec3 startpos = psV[vert_idx];

    activeGizmo = true;
    if (!vertexGizmo) {
    vertexGizmo = polyscope::addTransformationGizmo("vertex_editor");
    vertexGizmo->setAllowTranslation(true);
    vertexGizmo->setAllowRotation(false);
    vertexGizmo->setAllowScaling(false);
    vertexGizmo->setInteractInLocalSpace(false);
    //vertexGizmo->setGizmoSize(0.5f);
    }

    // place gizmo at vertex position
    vertexGizmo->setPosition(startpos);
}

// Remove gizmo
void endVertexEdit() {
    if (!activeGizmo) return;

    if (vertexGizmo) {
        vertexGizmo->remove();
        vertexGizmo = nullptr;
    }
    activeGizmo = false;
    return;
}

// Add specified constraint of a vertex to the constraint vector
bool addConstraintToVector(int v_idx, Eigen::Vector3d constraint_direc) {    // Need to generalize this to all dims.
    // If our vertex is already fully constrained, or the constraint direction is essentially 0, ignore.
    if (constraints[v_idx].size() >= 3) {
        std::cout << "Vertex is already fully constrained. Did not add to constraint list." << std::endl;
        return false;
    } else if (constraint_direc.norm() <= 1e-5) {
        std::cout << "Constraint is essentially 0. Did not add to constraint list." << std::endl;
        return false;
    }
    // Normalize in place
    constraint_direc.normalize();
    Eigen::Vector3d proj = constraint_direc;
    // Check which vectors are not already represented
    for (int j = 0; j < constraints[v_idx].size(); j++) {
        // Project out each existing constraint
        proj -= (constraint_direc.dot(constraints[v_idx][j])) * constraints[v_idx][j];
        if (proj.norm() == 1e-5) {  // If projection yielded 0, then we have a dependent vector
            std::cout << "Constraint is linearly dependent on existing constraints. Did not add to constraint list." << std::endl;
            return false;
        }
    }
    proj.normalize();
    // Add constraint vectors in, so long as we don't exceed 3
    int num_total = constraints[v_idx].size();
    // If we would have 3 constraints and they are all orthogonal, then
    // just replace with canonical axes for computational simplicity/stability
    if (num_total == 2) {
        std::cout << "Constraint added. Vertex " << v_idx << " is now fully constrained." << std::endl;
        constraints[v_idx][0] = i0;
        constraints[v_idx][1] = i1;
        constraints[v_idx].push_back(i2);
    } else {
        constraints[v_idx].push_back(proj);
        std::cout << "Constraint added. Vertex " << v_idx << " now has " << num_total+1 << " constraints." << std::endl;
    }
    return true;
}
// Wipe all constraints
void resetConstraintVector(bool resetPC) {
    // Clear Constraints
    for (int i = 0; i < constraints.size(); i++) {
        constraints[i].clear();
    }
    if (resetPC) {
        // Clear point cloud
        psConstraints.clear();
        VtoConstraint.clear();
        std::cout << "Cleared all constraints." << std::endl;
        psConstraintPC = polyscope::registerPointCloud("Constraints", psConstraints);
    }
    return;
}

// For the direct solver, test whether there are at least three constrained vertices
// If so, constrain in all direcs and return true
bool testDirectSolverConstraints() {
    std::vector<int> fullyConstrain;
    // Get a list of which vertices are constrained
    for (int v = 0; v < constraints.size(); v++) {
        if (constraints[v].size() != 0) {
            fullyConstrain.push_back(v);
        }
    }
    // Terminate here, if not enough constrained
    if (fullyConstrain.size() < 3) {
        return false;
    }
    // Apply rigid constraint to all directions for constrained vertices
    resetConstraintVector(false);
    std::vector<int> allDirecs = {0, 1, 2};
    for (int v = 0; v < fullyConstrain.size(); v++) {
        addConstraintToVector(fullyConstrain[v], i0);
        addConstraintToVector(fullyConstrain[v], i1);
        addConstraintToVector(fullyConstrain[v], i2);
    }
    return true;
}

// Modify a vertex position in polyscope
void modifyVertexPositions(int vert_idx, glm::vec3 new_pos) {
    psV[vert_idx] = new_pos;
    V_current[vert_idx] = IO::glmToEigen(new_pos);
    return;
}

// Reset all vertex positions in Polyscope, as well as the constraints
void resetVertexPositions() {
    // Reset Mesh
    IO::polyscopeTetConverter(V, T, psV, psT);
    IO::loadTOBJ(InputPath, V_current, T_current, false);
    psMesh = polyscope::registerTetMesh("Simulation Mesh", psV, psT);
    // Reset constraint PC
    for (const auto& pair : VtoConstraint) {
        psConstraints[pair.second] = psV[pair.first];
    }
    psConstraintPC = polyscope::registerPointCloud("Constraints", psConstraints);

    //resetConstraintVector(true);
    return;
}

// A user-defined callback, for creating control panels (etc)
void myCallback() {
    ImGuiIO& io = ImGui::GetIO();
    float totalWidth = ImGui::GetContentRegionAvail().x;
    float spacing = ImGui::GetStyle().ItemSpacing.x;
    float itemWidthThird1 = (totalWidth - spacing) * 0.3f;
    float itemWidthThird2 = (totalWidth - spacing) * 0.27f;
    
    bool mouseDown = ImGui::IsMouseDown(0);
    bool mouseClicked = ImGui::IsMouseClicked(0);
    bool mouseReleased = ImGui::IsMouseReleased(0);
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
            gizmoMode = false;
            endVertexEdit();

            // If using direct solver, we need to make sure there is at least one constraint
            if (sIdx == 0) {
                bool success = testDirectSolverConstraints();
                if (!success) {
                    std::cout << "Direct Solver failed. At least three constrained vertices must be specified to produce a non-singular system." << std::endl;
                    return;
                }
                // Direct solve has no time steps, so just compute 1 frame
                playbackStartIdx = 0;
                playbackEndIdx = 1;
                seconds = 1;
                fps = 1;
            } 
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

    if (ImGui::Button("Reset All Settings")) {
        // Turn off all modes first
        endVertexEdit();
        gizmoMode = false;
        scrubMode = false;
        constraintMode = false;
        loadingMode = false;
        // Reset variables
        PR = 0.35;
        YM = 5e4;
        seconds = 10;
        fps = 5;
        spf = 1.0/fps;
        int sIdx = 2;
        eIdx = 0;
        R = true;
        G = true;
        Gx = 0;
        Gy = 0;
        Gz = -9.8;
        customConstraint = false;
        Cx = 1.0;
        Cy = 1.0;
        Cz = 1.0;
        std::cout << "All settings reset. NOTE: Constraints and Mesh Vertices are NOT reset." << std::endl;
    }

    if (ImGui::Button(gizmoMode ? "Stop Moving Vertex" : "Move Vertex")) {
        if (loadingMode) {
            std::cout << "Currently in loading mode. Run without flags to solve for a deformation." << std::endl;
        } else if (autoPlaying) {
            std::cout << "Currently in Autoplaying mode. Stop Playback to move vertices." << std::endl;
        } else if (constraintMode) {
            std::cout << "Currently in Constraint mode. Stop applying constraints to move vertices." << std::endl;
        }
        else {
            gizmoMode = !gizmoMode;
            selectedVertex = -1;
            endVertexEdit();
            if (gizmoMode) {
                std::cout << "Hiding Constraints during editing." << std::endl;
                showPC = false;
                psConstraintPC->setEnabled(showPC);
            } else {
                std::cout << "Unhiding Constraints." << std::endl;
                showPC = true;
                psConstraintPC->setEnabled(showPC);
            }
        }
    }

    ImGui::SameLine();
    if (ImGui::Button("Reset Vertices")) {
        if (loadingMode) {
            std::cout << "Currently in loading mode. Run without flags to solve for a deformation." << std::endl;
        } else if (autoPlaying) {
            std::cout << "Currently in Autoplaying mode. Stop Playback to reset vertices." << std::endl;
        } else if (constraintMode) {
            std::cout << "Currently in Constraint mode. Stop applying constraints to reset vertices." << std::endl;
        }
        std::cout << "Resetting vertex locations. Note that constraints will NOT be reset." << std::endl;
        resetVertexPositions();
        gizmoMode = false;
        showPC = true;
        psConstraintPC->setEnabled(showPC);
    }


    if (ImGui::Button(constraintMode ? "Stop Applying Constraints" : "Apply Constraint")) {
        if (loadingMode) {
            std::cout << "Currently in loading mode. Run without flags to solve for a deformation." << std::endl;
        } else if (autoPlaying) {
            std::cout << "Currently in Autoplaying mode. Stop Playback to add constraints." << std::endl;
        } else if (gizmoMode) {
            std::cout << "Currently in Vertex Moving mode. Stop this mode to add constraints." << std::endl;
            std::cout << "NOTE: It is recommended that you move vertices before applying constraints." << std::endl;
        } else {
            constraintMode = !constraintMode;
            if (constraintMode) {
                std::cout << "Note: Use Axis checkboxes to select which axes to constrain." << std::endl;
            }
            selectedVertex = -1;
            showPC = true;
            psConstraintPC->setEnabled(showPC);
        }
    }

    ImGui::SameLine();
    // Reset the constraint set
    if (ImGui::Button("Reset Constraints")) {
        constraintMode = false;
        resetConstraintVector(true);
    }

    // Which axes to constrain
    ImGui::Text("Constrain Axes:");
    ImGui::SameLine();
    ImGui::Checkbox("X", &xAxis);
    ImGui::SameLine();
    ImGui::Checkbox("Y", &yAxis);
    ImGui::SameLine();
    ImGui::Checkbox("Z", &zAxis);
    ImGui::SameLine();
    ImGui::Checkbox("Custom", &customConstraint);

    if (customConstraint) { // Custom constraint should exist on its own
        xAxis = false;
        yAxis = false;
        zAxis = false;
    }
    ImGui::PushItemWidth(itemWidthThird1);
    ImGui::SliderFloat("x", &Cx, -1.0, 1.0, "%.2f");
    ImGui::PopItemWidth();
    ImGui::SameLine();
    ImGui::PushItemWidth(itemWidthThird1);
    ImGui::SliderFloat("y", &Cy, -1.0, 1.0, "%.2f");
    ImGui::PopItemWidth();
    ImGui::SameLine();
    ImGui::PushItemWidth(itemWidthThird1);
    ImGui::SliderFloat("z", &Cz, -1.0, 1.0, "%.2f");
    ImGui::PopItemWidth();

    // Gizmo Mode activation
    if (gizmoMode) {
        // Gizmo moved
        if (activeGizmo) {
            glm::vec3 newPos = vertexGizmo->getPosition();
            modifyVertexPositions(selectedVertex, newPos);
            psMesh->updateVertexPositions(psV);

            // Also update constraint vert if is constrained
            auto it = VtoConstraint.find(selectedVertex);
            if (it != VtoConstraint.end()) {
                psConstraints[it->second] = newPos;
                psConstraintPC->updatePointPositions(psConstraints);
            }
        }
        // New vertex selected
        if (mouseClicked && pick.isHit && pick.structure == psMesh) {
            polyscope::VolumeMeshPickResult meshPick = psMesh->interpretPickResult(pick);

            if (meshPick.elementType == polyscope::VolumeMeshElement::VERTEX) {
                int tempSelectedVertex = static_cast<int>(meshPick.index);
                // If we are switching to a different vertex, then remove the current gizmo and create a new one
                if ((tempSelectedVertex != selectedVertex)) {
                    endVertexEdit();
                    selectedVertex = tempSelectedVertex;
                    std::cout << "Moving Vertex " << selectedVertex << std::endl;
                    beginVertexEdit(selectedVertex);
                }
            } else if (io.WantCaptureMouse || ImGui::IsAnyItemActive() || ImGui::IsMouseDragging(0)) {
                // If user is dragging gizmo (or any UI item is active), do NOT dismiss
                return; // skip click-away logic this frame
            } else {
                std::cout << "Did not click on handle or gizmo." << std::endl;
                selectedVertex = -1;
                endVertexEdit();    // Remove gizmo
            }
        }
    }

    // Constraint mode activation
    if (constraintMode && mouseClicked && pick.isHit) {
        bool hit_object = false;
        // If we hit a mesh, get the mesh vertex we hit
        if (pick.structure == psMesh) {
            polyscope::VolumeMeshPickResult meshPick = psMesh->interpretPickResult(pick);
            if (meshPick.elementType == polyscope::VolumeMeshElement::VERTEX) {
                selectedVertex = static_cast<int>(meshPick.index);
                hit_object = true;
            } else {
                std::cout << "Did not click on mesh vertex." << std::endl;
            }
        } else if (pick.structure == psConstraintPC) {  // If we hit the point cloud, get the underlying mesh vertex
            polyscope::PointCloudPickResult pcPick = psConstraintPC->interpretPickResult(pick);
            selectedVertex = static_cast<int>(pcPick.index);
            // Figure out which point cloud point we hit
            // TODO: This is quite slow if we have tons of constraints...
            for (const auto& cVert : VtoConstraint) {
                if (cVert.second == selectedVertex) {
                    selectedVertex = cVert.first;
                    break;
                }
            }
            hit_object = true;
        } else {    // If we did not hit either object, quit
            std::cout << "Did not click on mesh vertex." << std::endl;
        }
        // Add constraint if we hit a mesh vertex
        if (hit_object) {
            bool success;
            int num_constraints_before = constraints[selectedVertex].size();
            // Add the specified constraint(s)
            if (customConstraint) {
                Eigen::Vector3d new_constraint = {static_cast<double>(Cx),
                                                  static_cast<double>(Cy),
                                                  static_cast<double>(Cz)};
                success = addConstraintToVector(selectedVertex, new_constraint);
            } else {
                bool successX = false;
                bool successY = false;
                bool successZ = false;
                if (xAxis) {
                    successX = addConstraintToVector(selectedVertex, i0);
                } if (yAxis) {
                    successY = addConstraintToVector(selectedVertex, i1);
                } if (zAxis) {
                    successZ = addConstraintToVector(selectedVertex, i2);
                }
                success = (successX || successY || successZ);
            }
            // Add to constraint vector
            // Revisualize our constrained vertices if any new ones appeared
            if ((num_constraints_before == 0) && success) {
                std::cout << "Added Vertex " << selectedVertex << " to constraint set." << std::endl;
                // Create point and add to point cloud
                glm::vec3 v_xyz = psV[selectedVertex];
                psConstraints.emplace_back(v_xyz);
                psConstraintPC = polyscope::registerPointCloud("Constraints", psConstraints);
                psConstraintPC->setPointRadius(0.008);

                VtoConstraint[selectedVertex] = psConstraints.size()-1;
            }
            // Update the vector visualization
            if (success) {
                displayConstraintDirecs();
            }
        }
    }

    // Simulation Time
    ImGui::SliderInt("Simulation Time (s)", &seconds, 1, 120);
    if (ImGui::SliderInt("Steps Per Second", &fps, 1, 100) && !autoPlaying) {   // We can only scrub when not autoplaying
        spf = 1.0/fps;
    }
    ImGui::Text("Timestep = %.3e", spf);

    // Simulator
    ImGui::RadioButton("Direct", &sIdx, 0);
    ImGui::SameLine();
    ImGui::RadioButton("BDF-1", &sIdx, 1);
    ImGui::SameLine();
    ImGui::RadioButton("BDF-2", &sIdx, 2);

    // If clicked Direct Solver, let's just go ahead and nix the other forces
    if (sIdx == 0) {
        R = false;
        G = false;
    }

    // Energies
    ImGui::RadioButton("SNH", &eIdx, 0);
    ImGui::SameLine();
    ImGui::RadioButton("ARAP", &eIdx, 1);

    // Material Parameters
    ImGui::SliderFloat("Poisson's Ratio", &PR, 0., 0.499, "%.3f");
    if (ImGui::SliderFloat("Young's Mod. (log10)", &logE, 4.0f, 6.0f, "%.2f")) {
        YM = std::pow(10.0, logE);
    }
    ImGui::Text("E = %.3e", YM);

    // Force Parameters
    ImGui::Checkbox("Rayleigh Damping", &R);
    ImGui::SameLine();
    ImGui::Checkbox("External Force", &G);

    ImGui::PushItemWidth(itemWidthThird2);
    ImGui::SliderFloat("F(x)", &Gx, -10.0, 10.0, "%.2f");
    ImGui::PopItemWidth();
    ImGui::SameLine();
    ImGui::PushItemWidth(itemWidthThird2);
    ImGui::SliderFloat("F(y)", &Gy, -10.0, 10.0, "%.2f");
    ImGui::PopItemWidth();
    ImGui::SameLine();
    ImGui::PushItemWidth(itemWidthThird2);
    ImGui::SliderFloat("F(z)", &Gz, -10.0, 10.0, "%.2f");
    ImGui::PopItemWidth();
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
    std::cout << "\nLoading TOBJ file" << std::endl;
    IO::loadTOBJ(InputPath, V, T, false);
    IO::loadTOBJ(InputPath, V_current, T_current, false);
    std::cout << "Loading TOBJ data to internal tetmesh object" << std::endl;
    m = std::make_unique<Geom::tetmesh>(V, T);
    // Initialize constraints
    constraints.resize(m->n);

    // Tetrahedralize OBJ inputs
    // TODO
    /*
    if (false) {
        Eigen::MatrixXd V;
        Eigen::MatrixXi F;
        Eigen::MatrixXi TF;
        // Read file using libigl
        // igl::copyleft::tetgen::tetrahedralize(psV, psF,"pq1.414Y", TV,TT,TF);
    }
    */

    std::cout << "\nConverting to Polyscope TetMesh" << std::endl;
    IO::polyscopeTetConverter(V, T, psV, psT);

    // Register tet mesh with PS
    std::cout << "Registering Tet Mesh to Polyscope" << std::endl;
    psMesh = polyscope::registerTetMesh("Simulation Mesh", psV, psT);

    // Empty point cloud for constraints
    psConstraintPC = polyscope::registerPointCloud("Constraints", psConstraints);
    psConstraintPC->setEnabled(showPC);

    std::cout << "\nTotal Vertex Count: " << m->n << std::endl;
    std::cout << "Total Tet Count: " << m->num_t << "\n" << std::endl;

    // Give control to the polyscope gui
    polyscope::show();

    return EXIT_SUCCESS;
}