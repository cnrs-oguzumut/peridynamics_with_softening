#pragma once
#include <vector>
#include <array>
#include <cmath>
#include <iostream>

inline std::vector<bool> build_node_in_defect(
    const std::vector<std::array<double, 2>>& coord_bar,
    bool introduce_defect,
    double defect_x_center,
    double defect_y_center,
    double defect_radius,
    double s_initial,
    double defect_strength_factor,
    double& defect_s_initial_out
) {
    const int totnode = static_cast<int>(coord_bar.size());
    std::vector<bool> node_in_defect(totnode, false);

    defect_s_initial_out = s_initial * defect_strength_factor;

    if (!introduce_defect) return node_in_defect;

    std::cout << "\n=== DEFECT INTRODUCTION ===\n";
    std::cout << "Defect center: (" << defect_x_center << ", " << defect_y_center << ")\n";
    std::cout << "Defect radius: " << defect_radius << " m\n";
    std::cout << "Strength factor: " << defect_strength_factor
              << " (bonds fail at " << (defect_strength_factor * 100) << "% of normal)\n";

    int defect_node_count = 0;
    for (int i = 0; i < totnode; ++i) {
        const double dx = coord_bar[i][0] - defect_x_center;
        const double dy = coord_bar[i][1] - defect_y_center;
        const double dist = std::sqrt(dx * dx + dy * dy);
        if (dist <= defect_radius) {
            node_in_defect[i] = true;
            ++defect_node_count;
        }
    }

    std::cout << "Nodes in defect zone: " << defect_node_count << "\n";
    std::cout << "Normal critical stretch: " << s_initial << "\n";
    std::cout << "Defect critical stretch: " << defect_s_initial_out << "\n";
    std::cout << "===========================\n\n";

    return node_in_defect;
}
