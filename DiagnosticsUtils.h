#pragma once
#include <vector>
#include <array>
#include <fstream>
#include <string>
#include <cmath>

struct DiagnosticsState {
    int target_node = -1;
    std::vector<int> diagnostic_nodes;

    std::ofstream diag_stress;
    std::ofstream diag_stretch;
    std::ofstream diag_pforce;
    std::ofstream diag_disp;

    std::ofstream diag_stress_allsteps;
    std::ofstream diag_pforce_allsteps;
    std::ofstream diag_stretch_allsteps;
    std::ofstream diag_disp_allsteps;
};

inline int find_nearest_x_node(const std::vector<std::array<double, 2>>& coord_bar, double target_x) {
    int best = 0;
    double best_dist = 1e100;
    for (int i = 0; i < static_cast<int>(coord_bar.size()); ++i) {
        const double d = std::abs(coord_bar[i][0] - target_x);
        if (d < best_dist) {
            best_dist = d;
            best = i;
        }
    }
    return best;
}

inline void init_diagnostics_nodes_if_enabled(
    DiagnosticsState& ds,
    bool track_numerical_noise,
    const std::vector<std::array<double, 2>>& coord_bar,
    int ndivx,
    double target_x_coord_bar
) {
    ds.target_node = find_nearest_x_node(coord_bar, target_x_coord_bar);
    ds.diagnostic_nodes.clear();

    if (!track_numerical_noise) return;

    // Keep your same logic here. If your original code used center + neighbors, replicate it:
    // Example: center node and +/- some offsets. Adjust to match your existing behavior.
    const int center = ds.target_node;
    ds.diagnostic_nodes.push_back(center);

    // simple symmetric neighbor selection (change to match original if needed)
    const int left = std::max(0, center - 1);
    const int right = std::min(static_cast<int>(coord_bar.size()) - 1, center + 1);
    if (left != center) ds.diagnostic_nodes.push_back(left);
    if (right != center && right != left) ds.diagnostic_nodes.push_back(right);
}

inline void open_allsteps_files_if_needed(
    DiagnosticsState& ds,
    bool track_numerical_noise,
    bool do_snap,
    const std::string& tag,
    const std::vector<std::array<double, 2>>& coord_bar
) {
    if (!track_numerical_noise || !do_snap) return;

    const std::string suffix = "_at_loop=" + tag + ".txt";
    ds.diag_stress_allsteps.open("diagnostic_stress_all_ADR_steps" + suffix);
    ds.diag_pforce_allsteps.open("diagnostic_pforce_all_ADR_steps" + suffix);
    ds.diag_stretch_allsteps.open("diagnostic_max_stretch_all_ADR_steps" + suffix);
    ds.diag_disp_allsteps.open("diagnostic_disp_bar_all_ADR_steps" + suffix);

    ds.diag_stress_allsteps << "# ADR_step";
    ds.diag_pforce_allsteps << "# ADR_step";
    ds.diag_stretch_allsteps << "# ADR_step";
    ds.diag_disp_allsteps << "# ADR_step";

    for (int ni : ds.diagnostic_nodes) {
        ds.diag_stress_allsteps << " node" << ni << "(x=" << coord_bar[ni][0] << ")";
        ds.diag_pforce_allsteps << " node" << ni;
        ds.diag_stretch_allsteps << " node" << ni;
        ds.diag_disp_allsteps << " node" << ni;
    }
    ds.diag_stress_allsteps << '\n';
    ds.diag_pforce_allsteps << '\n';
    ds.diag_stretch_allsteps << '\n';
    ds.diag_disp_allsteps << '\n';
}

inline void close_allsteps_files(DiagnosticsState& ds) {
    if (ds.diag_stress_allsteps.is_open()) ds.diag_stress_allsteps.close();
    if (ds.diag_pforce_allsteps.is_open()) ds.diag_pforce_allsteps.close();
    if (ds.diag_stretch_allsteps.is_open()) ds.diag_stretch_allsteps.close();
    if (ds.diag_disp_allsteps.is_open()) ds.diag_disp_allsteps.close();
}