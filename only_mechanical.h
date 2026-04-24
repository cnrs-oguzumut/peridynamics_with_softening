#include <vector>
#include <array>
#include <fstream>
#include <string>
#include <filesystem>
#include "DefectUtils.h"
#include "DiagnosticsUtils.h"
inline void mechanical(
    // geometry / sizes
    double dx_bar,
    int ndivx,
    int totnode,
    int totint,


    // time/load control
    int num_loop,
    double& epsilonn, // updated each outer increment
    double del_epsilonn,
    double dt_mechanical_bar,
    int nt_mechanical,
    double tolerance,
    int max_iteration,


    // material / constants
    double stress_coeff,
    double s_initial,
    double s_c,
    double s_k,
    double beta_softening,
    double Rd,
    int degradation_model,
    double vol_bar,
    double bc_mechanical_film_bar, //!it should be updated for film_substrate system
    double& cn,
     double c_int_ratio, double h_interface_bar, double delta_bar,  // NEW

    // fields (in/out)
    const std::vector<std::array<double,2>>& coord_bar,
    std::vector<std::array<double,2>>& disp_bar,
    std::vector<std::array<double,2>>& vel_bar,
    std::vector<std::array<double,2>>& vel_barhalf,
    std::vector<std::array<double,2>>& vel_barhalfold,
    std::vector<std::array<double,2>>& pforce_mechanical,
    std::vector<std::array<double,2>>&  massvec_mechanical,
    std::vector<std::array<double,2>>& bforce_mechanical,
    std::vector<std::array<double,2>>& pforceold,
    std::vector<double>& strain,
    std::vector<double>& dmg,

    // topology
    std::vector<int>& numfam,
    std::vector<int>& pointfam,
    int* nodefam,
    int** fail,
    double** fac,
    double** scr_mechanical,
    double** bond_s_history,
    double** bond_dix,
    double** bond_diy,
    double** bond_idist,
    double** bond_inv_idist,
    double** bond_weight,
    int**    bond_is_sub,


    // energies / outputs
    std::vector<double>& strain_energy,
    std::vector<double>& kinetic_energy,
    std::vector<double>& fracture_energy,
    std::vector<double>& damping_dissipation_energy,
    std::vector<double>& substrate_energy,
    std::vector<double>& total_energy_density,

    std::vector<double>& sum_ac_of_stress,
    std::vector<std::array<double,4>>& stress,
    const int* write_snaps, int write_snaps_count,


    /* ovito start gate */
    int ovito_start_increment,
    double target_x_coord_bar,  // NEW: x-coord_barinate to track

    /* noise tracking flag */
    bool track_numerical_noise,

    /* stress convergence tracking flag */
    bool track_stress_convergence,

    /* fixed nt_mechanical flag */
    bool fixed_nt_mechanical,

    /* fixed damping_coefficient flag */
    bool fixed_damping_coefficient,

    /* defect parameters */
    bool introduce_defect,
    double defect_x_center,
    double defect_y_center,
    double defect_radius,
    double defect_strength_factor,
    bool stop_after_first_crack,
    int& first_crack_increment_out
)

{
//helper
auto want_snapshot = [&](int k) -> bool {
    for (int i = 0; i < write_snaps_count; ++i)
        if (write_snaps[i] == k) return true;
    return false;
};

// Find the node index in the film closest to target_x_coord_bar: we want to track ADR equilibirium
int target_node = -1;
double min_dist = 1e10;
for (int i = 0; i < totnode; ++i) {
    double dist = std::abs(coord_bar[i][0] - target_x_coord_bar);
    if (dist < min_dist) {
        min_dist = dist;
        target_node = i;
    }
}

//! ==================================================================
//! DEFECT: NODE-BASED CRITICAL STRETCH
//! ==================================================================
double defect_s_initial = s_initial;
std::vector<bool> node_in_defect = build_node_in_defect(
    coord_bar,
    introduce_defect,
    defect_x_center,
    defect_y_center,
    defect_radius,
    s_initial,
    defect_strength_factor,
    defect_s_initial
);
//! ==================================================================
//! To check numerical noise--------------
DiagnosticsState diag;
init_diagnostics_nodes_if_enabled(
    diag,
    track_numerical_noise,
    coord_bar,
    ndivx,
    target_x_coord_bar
);
if (track_numerical_noise) {
    diag.diag_stress.open("diagnostic_stress.txt", std::ios::trunc);
    diag.diag_pforce.open("diagnostic_pforce.txt", std::ios::trunc);
    diag.diag_disp.open("diagnostic_disp_bar.txt", std::ios::trunc);
    diag.diag_stretch.open("diagnostic_max_stretch.txt", std::ios::trunc);
}
//!End:To check numerical noise------------------------------------------

// Open convergence tracking file
std::ofstream convergence_file;
if (track_stress_convergence) {
    convergence_file.open("convergence_of_stress.txt");
    //convergence_file << "# Each row = one outer loop\n";
    //convergence_file << "# Tracking stress[" << target_node << "][0] at x=" << coord_bar[target_node][0] << "\n";
}
std::ofstream outfile16;
outfile16.open("total_energy_density_matrix.txt");

std::ofstream outfile18;
outfile18.open("total_kinetic_energy_density_matrix.txt");

std::ofstream outfile19;
outfile19.open("total_damping_dissipation_energy_density_matrix.txt");

std::ofstream outfile20;
outfile20.open("total_fracture_energy_density_matrix.txt");

std::ofstream outfile21;
outfile21.open("total_strain_energy_density_matrix.txt");

std::ofstream outfile22;
outfile22.open("total_substrate_energy_density_matrix.txt");


//std::ofstream outfile14;
//outfile14.open("cn.txt");//ploting damping coefficient for adaptive dynamic relaxation

std::ofstream outfile1("Average_stress_of_interface.txt");
outfile1 << "epsilon" << ' ' << "xx-stress" << '\n';
std::ofstream outfile2("Final_ADR_stress.txt");
std::ofstream outfile3("Final_ADR_strain.txt");
std::ofstream outfile4("Final_ADR_strain_energy.txt");
std::ofstream outfile5("Final_ADR_kinetic_energy.txt");
std::ofstream outfile6("Final_ADR_fracture_energy.txt");
std::ofstream outfile7("Final_ADR_damping_dissipation_energy.txt");
std::ofstream outfile17("Final_ADR_substrate_energy.txt");

// --- Early-stop tracking variables ---
int first_crack_increment = -1;   // will store the increment when first dmg > 0
int early_stop_increment  = -1;   // will store first_crack_increment + 5

// Buffers to accumulate ADR-step data for the cracking increment
// (written retroactively if this turns out to be the cracking increment)
std::vector<std::vector<double>> buf_stress, buf_strain, buf_strain_energy;
std::vector<std::vector<double>> buf_kinetic, buf_fracture, buf_damp, buf_sub, buf_total;
const std::string ovito_output_dir = "ovito_files";
const std::string result_output_dir = "result";
std::filesystem::create_directories(ovito_output_dir);
std::filesystem::create_directories(result_output_dir);

//double cn = 0.01;//0.01;// cn unit:1/s
const int max_extended_loop = 10000;  // upper bound when waiting for crack
const int effective_num_loop = (stop_after_first_crack) ? max_extended_loop : num_loop;
for (int loop_counter=0;loop_counter<effective_num_loop;loop_counter++)
{
    std::cout << "increment " << loop_counter << '\n';


    //damping_dissipation_energy[i] = 0.0;
    epsilonn += del_epsilonn;
    int counterADR = 0;
    double maxval = 0.0;

    buf_stress.clear(); buf_strain.clear(); buf_strain_energy.clear();
    buf_kinetic.clear(); buf_fracture.clear(); buf_damp.clear();
    buf_sub.clear(); buf_total.clear();

    double sum_interface_stress=0.0;
    double sum_interface_strain=0.0;

    // per-increment ovito file
    /*std::ofstream outfile15("for_ovito_" + std::to_string(loop_counter + 1) + ".xyz");
    outfile15 << totint << '\n' << '\n';*/



    bool do_snap = want_snapshot(loop_counter);
    if (do_snap) {
            const std::string tag = std::to_string(loop_counter);
            std::ofstream("stress"            + tag + ".txt", std::ios::trunc).close();
            std::ofstream("strain"            + tag + ".txt", std::ios::trunc).close();
            std::ofstream("strain_energy"     + tag + ".txt", std::ios::trunc).close();
            std::ofstream("kinetic_energy"    + tag + ".txt", std::ios::trunc).close();
            std::ofstream("fracture_energy"   + tag + ".txt", std::ios::trunc).close();
            std::ofstream("damping_dissipation_energy"+ tag + ".txt", std::ios::trunc).close();
            std::ofstream("substrate_energy"+ tag + ".txt", std::ios::trunc).close();
            std::ofstream("total_energy_density"+ tag + ".txt", std::ios::trunc).close();
             // Open diagnostic files for all ADR steps at this increment
       // Open diagnostic files for all ADR steps at this increment
   open_allsteps_files_if_needed(
    diag,
    track_numerical_noise,
    do_snap,
    tag,
    coord_bar
);
    }

    // Continue each load increment from the converged displacement field of the
    // previous increment. Only dynamic relaxation state is reset, matching
    // ../pery/pery.py: the solution path is continuous, but each increment
    // starts relaxation from zero velocity.
    for (int i = 0;i<totnode;i++){
        pforceold[i][0] = 0.0;
        vel_bar[i][0] = 0.0;
        vel_bar[i][1] = 0.0;
       // velhalf[i][0] = 0.0;
       // velhalf[i][1] = 0.0;
//        velhalfold[i][0] = 0.0;
      //  velhalfold[i][1] = 0.0;
    }

//!ADR inner loop------------------------------------------------------
for (int tt_mech = 1; (fixed_nt_mechanical ? (tt_mech <= nt_mechanical) : (tt_mech <= max_iteration)); tt_mech++) { //to control the loop without tolerance
   // while (true){
        double total_energy_of_system=0.0;
        double total_kinetic_energy_of_system =0.0;
        double total_fracture_energy_of_system =0.0;
        double total_strain_energy_of_system =0.0;
        double total_substrate_energy_of_system =0.0;
        double total_damping_dissipation_energy_of_system=0.0;

        maxval=0.0;
        double cn1 = 0.0;
        double cn2 = 0.0;
    //    std::cout << "counterADR " << counterADR << '\n';
    //maxval=1.0;
    //cout<<"counterADR "<<counterADR<<endl;
   // const double ctime_mech=(counterADR+1)*dt_mechanical_bar;

    //!BCs: Left & Right load
    ///bforce_mechanical[0][0]=-10.0e9*epsilonn/dx;//-appres/dx;//
    ///bforce_mechanical[ndivx-1][0]=10.0e9*epsilonn/dx;//appres/dx;//


    const int prescribed_boundary_nodes = 3;
    for (int k = 0; k < prescribed_boundary_nodes && k < totnode; ++k) {
        const int left = k;
        const int right = totnode - 1 - k;

        disp_bar[left][0] = epsilonn * coord_bar[left][0];
        vel_bar[left][0] = 0.0;

        if (right != left) {
            disp_bar[right][0] = epsilonn * coord_bar[right][0];
            vel_bar[right][0] = 0.0;
        }
    }


    //double vel_barocity=20.0;
    /*for (int i = totint;i<totint+nbnd_mechanical;i++){
        //bforce_mechanical[i][0]=-appres/dx;

        //vel_bar[i][0] = -vel_barocity;
        disp_bar[i][0] = epsilonn*coord_bar[i][0];//-vel_barocity* tt_mech * dt_mechanical_bar;
    }
    for (int i = totleft;i<totleft+nbnd_mechanical;i++){
        //bforce_mechanical[i][0]=appres/dx;
        //mark[i]=10;
      // vel_bar[i][0] = vel_barocity;
        disp_bar[i][0] = epsilonn*coord_bar[i][0];//vel_barocity* tt_mech * dt_mechanical_bar;
    }*/

    #pragma omp parallel
    {
    //!BCs: strain applied on the substrate
   /* #pragma omp for schedule(static)
    for (int i=0;i<n_substrate;i++){
        disp_bar[i][0]=epsilonn*coord_bar[i][0];
    }*/

    //!loop over the film material points
    #pragma omp for schedule(static) reduction(+:sum_interface_stress,sum_interface_strain)
    for (int i = 0;i<totnode;i++){
        double dmgpar1=0.0;
        double dmgpar2=0.0;
        pforce_mechanical[i][0] = 0.0;
        pforce_mechanical[i][1] = 0.0;
        //
        sum_ac_of_stress[i]=0.0;

        const bool prescribed_force_node =
            (i < prescribed_boundary_nodes) ||
            (i >= totnode - prescribed_boundary_nodes);
        if (prescribed_force_node) {
            dmg[i] = 0.0;
            strain_energy[i] = 0.0;
            substrate_energy[i] = 0.0;
            stress[i][0] = 0.0;
            continue;
        }
        //sum_ad_of_stress[i]=0.0;
        //sum_bc_of_stress[i]=0.0;
        //sum_bd_of_stress[i]=0.0;
        strain_energy[i] = 0.0;
       // fracture_energy[i]=0.0;
        //sum_strain[i]=0.0;

        // cache row pointers once (only for the sake of optimization!)
        const std::array<double,2>&        coord_bar_i = coord_bar[i];
        const std::array<double,2>&        disp_bar_i = disp_bar[i];

        double u_substrate_at_i = epsilonn * coord_bar[i][0];

     //!loop over neighboring material points (j) inside horizon
        for (int j =0 ;j<numfam[i];j++){
                const int cnode = nodefam[pointfam[i]+j-1]-1;//neighbor (-1 is for the index)


                const std::array<double,2>&  coord_bar_j = coord_bar[cnode];
                const std::array<double,2>& disp_bar_j  = disp_bar[cnode];

                 // geometry: initial & current bond vectors
                const double dix       = bond_dix[i][j];
                const double diy       = bond_diy[i][j];
                const double idist     = bond_idist[i][j];
                const double inv_idist = bond_inv_idist[i][j];
                const double w_const   = bond_weight[i][j];     // vol_bar * scr * fac
                const bool   sub_nbr   = (bond_is_sub[i][j] != 0);
                //const double dix = coord_bar_j[0] - coord_bar_i[0];
                //const double diy = coord_bar_j[1] - coord_bar_i[1];

        // ==============================
        // 1. FILM-FILM BOND FORCE
        // ==============================

                const double nx  = (coord_bar_j[0] + disp_bar_j[0]) - (coord_bar_i[0] + disp_bar_i[0]);
                const double ny  = (coord_bar_j[1] + disp_bar_j[1]) - (coord_bar_i[1] + disp_bar_i[1]);

                // distances (avoid pow; square by multiply)
                //const double idist2 = dix*dix + diy*diy;
                const double nlen2  = nx*nx;// + ny*ny;
                //const double idist  = std::sqrt(idist2);
                const double nlength   = std::sqrt(nlen2);

                //const double inv_idist = 1.0 / idist;
                const double inv_nlen  = 1.0 / nlength;

                const double stretch = (nlength - idist) * inv_idist;

                double& s_history = bond_s_history[i][j];
                if (stretch > s_history) {
                    s_history = stretch;
                }
                const double s_eff = s_history;

                double mu_history = 1.0;     // envelope force ratio relative to c * s_initial
                double f_env_history = 0.0;  // envelope force magnitude at s_history
                double k_unload = 0.0;       // unloading/reloading secant stiffness
                bool on_envelope = (std::abs(stretch - s_history) < 1.0e-14);

                double dforce1_mechanical, d_strain_energy, d_fracture_energy;


                /*idist = std::sqrt(
                std::pow(coord_bar[cnode][0] - coord_bar[i][0], 2) +
                std::pow(coord_bar[cnode][1] - coord_bar[i][1], 2));
                nlength = std::sqrt(
                std::pow((coord_bar[cnode][0] + disp_bar[cnode][0]) - (coord_bar[i][0] + disp_bar[i][0]), 2) +
                std::pow((coord_bar[cnode][1] + disp_bar[cnode][1]) - (coord_bar[i][1] + disp_bar[i][1]), 2));
                stretch = (nlength - idist) / idist;*/

        if (fail[i][j] == 1) { // bond not fully failed yet

        const double dirx = nx * inv_nlen;

        if (degradation_model < 0) {
            // Pure elastic mode: matches the Python reference solver.
            dforce1_mechanical =
                bc_mechanical_film_bar * stretch * w_const * dirx;

            const double s2 = stretch * stretch;
            d_strain_energy =
                0.25 * bc_mechanical_film_bar * s2 * idist * w_const;
            d_fracture_energy = 0.0;

            sum_ac_of_stress[i] += dforce1_mechanical * dix;

            dmgpar1 += vol_bar * fac[i][j];
            dmgpar2 += vol_bar * fac[i][j];
        }
        else if (degradation_model == 0) {
            // ==================================================
            // Original brittle model
            // ==================================================
            if (stretch <= s_initial) {
                dforce1_mechanical =
                    bc_mechanical_film_bar * stretch * w_const * dirx;

                const double s2 = stretch * stretch;
                d_strain_energy =
                    0.25 * bc_mechanical_film_bar * s2 * idist * w_const;
                d_fracture_energy = 0.0;

                sum_ac_of_stress[i] += dforce1_mechanical * dix;

                dmgpar1 += vol_bar * fac[i][j];
                dmgpar2 += vol_bar * fac[i][j];
            }
            else {
                fail[i][j] = 0;
                dforce1_mechanical = 0.0;
                d_strain_energy = 0.0;

                const double cs2 = s_initial * s_initial;
                d_fracture_energy =
                    0.25 * bc_mechanical_film_bar * cs2 * idist * w_const;
                fracture_energy[i] += d_fracture_energy;

                dmgpar2 += vol_bar * fac[i][j];
            }
        }
        else if (degradation_model == 1) {
            // ==================================================
            // Linear softening model with unloading/reloading
            // ==================================================

            if (s_history <= s_initial) {
                // purely elastic history
                dforce1_mechanical =
                    bc_mechanical_film_bar * stretch * w_const * dirx;

                const double s2 = stretch * stretch;
                d_strain_energy =
                    0.25 * bc_mechanical_film_bar * s2 * idist * w_const;
                d_fracture_energy = 0.0;

                sum_ac_of_stress[i] += dforce1_mechanical * dix;

                dmgpar1 += vol_bar * fac[i][j];
                dmgpar2 += vol_bar * fac[i][j];
            }
            else if (s_history < s_c && std::abs(coord_bar[i][0] - 0.5) < 0.45) {
                // envelope value at the history stretch
                mu_history = (s_c - s_history) / (s_c - s_initial);
                f_env_history = bc_mechanical_film_bar * s_initial * mu_history;
                k_unload = f_env_history / s_history;

                if (on_envelope) {
                    // active loading on the softening envelope
                    dforce1_mechanical =
                        f_env_history * w_const * dirx;
                } else {
                    // unloading/reloading with degraded stiffness
                    dforce1_mechanical =
                        k_unload * stretch * w_const * dirx;
                }

                d_strain_energy =
                    0.25 * k_unload * stretch * stretch * idist * w_const;
                d_fracture_energy = 0.0;

                sum_ac_of_stress[i] += dforce1_mechanical * dix;

                dmgpar1 += mu_history * vol_bar * fac[i][j];
                dmgpar2 += vol_bar * fac[i][j];
            }
            else if (std::abs(coord_bar[i][0] - 0.5) < 0.45) {
                fail[i][j] = 0;
                dforce1_mechanical = 0.0;
                d_strain_energy = 0.0;

                const double cs2 = s_initial * s_initial;
                d_fracture_energy =
                    0.25 * bc_mechanical_film_bar * cs2 * idist * w_const;
                fracture_energy[i] += d_fracture_energy;

                dmgpar2 += vol_bar * fac[i][j];
            }
        }
           else if (degradation_model == 2) {
            // ==================================================
            // Bilinear softening model (trilinear) with unloading/reloading
            // ==================================================

            if (s_history <= s_initial) {
                dforce1_mechanical =
                    bc_mechanical_film_bar * stretch * w_const * dirx;

                const double s2 = stretch * stretch;
                d_strain_energy =
                    0.25 * bc_mechanical_film_bar * s2 * idist * w_const;
                d_fracture_energy = 0.0;

                sum_ac_of_stress[i] += dforce1_mechanical * dix;

                dmgpar1 += vol_bar * fac[i][j];
                dmgpar2 += vol_bar * fac[i][j];
            }
            else if (s_history < s_c) {
                if (s_history < s_k) {
                    mu_history =
                        1.0 - (1.0 - beta_softening) *
                        (s_history - s_initial) / (s_k - s_initial);
                } else {
                    mu_history =
                        beta_softening * (s_c - s_history) / (s_c - s_k);
                }

                f_env_history = bc_mechanical_film_bar * s_initial * mu_history;
                k_unload = f_env_history / s_history;

                if (on_envelope) {
                    dforce1_mechanical =
                        f_env_history * w_const * dirx;
                } else {
                    dforce1_mechanical =
                        k_unload * stretch * w_const * dirx;
                }

                d_strain_energy =
                    0.25 * k_unload * stretch * stretch * idist * w_const;
                d_fracture_energy = 0.0;

                sum_ac_of_stress[i] += dforce1_mechanical * dix;

                dmgpar1 += mu_history * vol_bar * fac[i][j];
                dmgpar2 += vol_bar * fac[i][j];
            }
            else {
                fail[i][j] = 0;
                dforce1_mechanical = 0.0;
                d_strain_energy = 0.0;

                const double cs2 = s_initial * s_initial;
                d_fracture_energy =
                    0.25 * bc_mechanical_film_bar * cs2 * idist * w_const;
                fracture_energy[i] += d_fracture_energy;

                dmgpar2 += vol_bar * fac[i][j];
            }
        }
           else if (degradation_model == 3) {
            // ==================================================
            // Cornelissen nonlinear softening model with unloading/reloading
            // ==================================================

            if (s_history <= s_initial) {
                dforce1_mechanical =
                    bc_mechanical_film_bar * stretch * w_const * dirx;

                const double s2 = stretch * stretch;
                d_strain_energy =
                    0.25 * bc_mechanical_film_bar * s2 * idist * w_const;
                d_fracture_energy = 0.0;

                sum_ac_of_stress[i] += dforce1_mechanical * dix;

                dmgpar1 += vol_bar * fac[i][j];
                dmgpar2 += vol_bar * fac[i][j];
            }
            else if (s_history < s_c) {
                const double eta =
                    (s_history - s_initial) / (s_c - s_initial);

                mu_history =
                    (1.0 + std::pow(3.0 * eta, 3.0)) * std::exp(-6.93 * eta)
                    - 28.0 * eta * std::exp(-6.93);

                mu_history = std::max(0.0, mu_history);

                f_env_history = bc_mechanical_film_bar * s_initial * mu_history;
                k_unload = f_env_history / s_history;

                if (on_envelope) {
                    dforce1_mechanical =
                        f_env_history * w_const * dirx;
                } else {
                    dforce1_mechanical =
                        k_unload * stretch * w_const * dirx;
                }

                d_strain_energy =
                    0.25 * k_unload * stretch * stretch * idist * w_const;
                d_fracture_energy = 0.0;

                sum_ac_of_stress[i] += dforce1_mechanical * dix;

                dmgpar1 += mu_history * vol_bar * fac[i][j];
                dmgpar2 += vol_bar * fac[i][j];
            }
            else {
                fail[i][j] = 0;
                dforce1_mechanical = 0.0;
                d_strain_energy = 0.0;

                const double cs2 = s_initial * s_initial;
                d_fracture_energy =
                    0.25 * bc_mechanical_film_bar * cs2 * idist * w_const;
                fracture_energy[i] += d_fracture_energy;

                dmgpar2 += vol_bar * fac[i][j];
            }
        }
          else if (degradation_model == 4) {
            // ==================================================
            // Exponential softening model with unloading/reloading
            // ==================================================

            if (s_history <= s_initial) {
                dforce1_mechanical =
                    bc_mechanical_film_bar * stretch * w_const * dirx;

                const double s2 = stretch * stretch;
                d_strain_energy =
                    0.25 * bc_mechanical_film_bar * s2 * idist * w_const;
                d_fracture_energy = 0.0;

                sum_ac_of_stress[i] += dforce1_mechanical * dix;

                dmgpar1 += vol_bar * fac[i][j];
                dmgpar2 += vol_bar * fac[i][j];
            }
            else {
                const double alpha =
                    (bc_mechanical_film_bar * s_initial) / Rd;

                mu_history =
                    std::exp(-alpha * (s_history - s_initial));

                f_env_history = bc_mechanical_film_bar * s_initial * mu_history;
                k_unload = f_env_history / s_history;

                if (on_envelope) {
                    dforce1_mechanical =
                        f_env_history * w_const * dirx;
                } else {
                    dforce1_mechanical =
                        k_unload * stretch * w_const * dirx;
                }

                d_strain_energy =
                    0.25 * k_unload * stretch * stretch * idist * w_const;
                d_fracture_energy = 0.0;

                sum_ac_of_stress[i] += dforce1_mechanical * dix;

                dmgpar1 += mu_history * vol_bar * fac[i][j];
                dmgpar2 += vol_bar * fac[i][j];
            }
        }
        else {
            // ==================================================
            // Unknown degradation model -> fall back to brittle
            // ==================================================
            if (stretch <= s_initial) {
                dforce1_mechanical =
                    bc_mechanical_film_bar * stretch * w_const * dirx;

                const double s2 = stretch * stretch;
                d_strain_energy =
                    0.25 * bc_mechanical_film_bar * s2 * idist * w_const;
                d_fracture_energy = 0.0;

                sum_ac_of_stress[i] += dforce1_mechanical * dix;

                dmgpar1 += vol_bar * fac[i][j];
                dmgpar2 += vol_bar * fac[i][j];
            }
            else {
                fail[i][j] = 0;
                dforce1_mechanical = 0.0;
                d_strain_energy = 0.0;

                const double cs2 = s_initial * s_initial;
                d_fracture_energy =
                    0.25 * bc_mechanical_film_bar * cs2 * idist * w_const;
                fracture_energy[i] += d_fracture_energy;

                dmgpar2 += vol_bar * fac[i][j];
            }
        }
    }
    else { // already failed bond
        dforce1_mechanical = 0.0;
        d_strain_energy = 0.0;
        d_fracture_energy = 0.0;
        dmgpar2 += vol_bar * fac[i][j];
    }

        // ==============================
        // 2. FILM-SUBSTRATE INTERFACE FORCE (Virtual Substrate)
        // ==============================
            // Virtual substrate position at x[cnode]
                /*double u_substrate_at_j = epsilonn * coord_bar[cnode][0];

                // Initial 2D distance (includes vertical offset h_interface)
                double dx_init = coord_bar_j[0] - coord_bar_i[0];
                double dist_init_substrate = std::sqrt(dx_init * dx_init + h_interface_bar * h_interface_bar);

                // Only compute interface force if within horizon
                if (dist_init_substrate <= delta_bar) {
                    // Deformed distance (substrate moves, film moves)
                    double dx_def = (coord_bar_j[0] + u_substrate_at_j) - (coord_bar_i[0] + disp_bar_i[0]);
                    double dist_def_substrate = std::sqrt(dx_def * dx_def + h_interface_bar * h_interface_bar);

                    // Interface bond stretch
                    double stretch_interface = (dist_def_substrate - dist_init_substrate) / dist_init_substrate;

                    // Force direction (only x-component contributes)
                    double dir_x_interface = dx_def / dist_def_substrate;

                    // Interface stiffness
                    double c_interface = bc_mechanical_film_bar * c_int_ratio;

                    // Interface force (add to total force)
                    dforce1_mechanical += c_interface * stretch_interface * dir_x_interface * w_const;

                    // Interface strain energy (add to total energy)
                    d_strain_energy += 0.25 * c_interface * stretch_interface * stretch_interface
                                     * dist_init_substrate * vol_bar;
        }*/


            //! Apply a different film-substrate stiffness by multiplying a coefficient
           /* if (sub_nbr){ //(coord_bar[cnode][1]-y_film)<=0.0
                dforce1_mechanical*=force_coeff;
            }*/

            //sum of the PD bond forces for material point i exerted by neighbors inside its horizon
            pforce_mechanical[i][0] += dforce1_mechanical;//N per voulme
            //pforce_mechanical[i][1] += dforce2_mechanical;
            strain_energy[i]+=d_strain_energy;//Strain Energy Density: energy per vol_barume
           // fracture_energy[i]+=d_fracture_energy;//energy per voulme


            //!Definition of a no-fail zone: No damage is allowed in the substrate
            // Use defect critical stretch if node is in defect zone
           /* double crit_stretch_local = s_initial;
            if (node_in_defect[i]) {
                crit_stretch_local = defect_s_initial;
            }*/

           // if (stretch > s_initial /*& abs(coord_bar[i][0])<0.25*/)
           //     fail[i][j] = 0;

               // dmgpar1 += fail[i][j] * vol_bar * fac[i][j];
              //  dmgpar2 += vol_bar * fac[i][j];

          //ac_of_stress=dforce1_mechanical*(coord_bar[cnode][0] - coord_bar[i][0]);//xx component of stress
          //ad_of_stress=dforce1_mechanical*(coord_bar[cnode][1] - coord_bar[i][1]);
          //bc_of_stress=dforce2_mechanical*(coord_bar[cnode][0] - coord_bar[i][0]);
          //bd_of_stress=dforce2_mechanical*(coord_bar[cnode][1] - coord_bar[i][1]);
         //   sum_ac_of_stress[i] += dforce1_mechanical * dix;
          //sum_ad_of_stress[i]+=ad_of_stress;
          //sum_bc_of_stress[i]+=bc_of_stress;
          //sum_bd_of_stress[i]+=bd_of_stress;
        }


        double substrate_force = 0.0;
        double substrate_energy_density = 0.0;
        const double c_interface = bc_mechanical_film_bar * c_int_ratio;
        for (int js = 0; js < totnode; ++js) {
            const double dx_init_substrate = coord_bar[js][0] - coord_bar_i[0];
            const double dist_init_substrate =
                std::sqrt(dx_init_substrate * dx_init_substrate + h_interface_bar * h_interface_bar);

            if (dist_init_substrate > delta_bar + 1.0e-12) continue;

            const double u_substrate_j = epsilonn * coord_bar[js][0];
            const double dx_def_substrate =
                (coord_bar[js][0] + u_substrate_j) - (coord_bar_i[0] + disp_bar_i[0]);
            const double dist_def_substrate =
                std::sqrt(dx_def_substrate * dx_def_substrate + h_interface_bar * h_interface_bar);
            const double stretch_substrate =
                (dist_def_substrate - dist_init_substrate) / dist_init_substrate;
            const double dir_x_substrate = dx_def_substrate / dist_def_substrate;

            substrate_force += c_interface * stretch_substrate * dir_x_substrate * vol_bar;
            substrate_energy_density +=
                0.25 * c_interface * stretch_substrate * stretch_substrate
                * dist_init_substrate * vol_bar;
        }
        pforce_mechanical[i][0] += substrate_force;

        dmg[i] = 1.0 - dmgpar1 / dmgpar2;
        //kinetic_energy[i]=0.5 * massvec_mechanical [i][0] * (vel_bar[i][0] * vel_bar[i][0] /*+ vel_bar[i][1] * vel_bar[i][1]*/);//energy per voulme
        substrate_energy[i] = substrate_energy_density;
        //damping_dissipation_energy[i]= 0.5 * massvec_mechanical[i][0] * cn * (pow(vel_bar[i][0],2) + pow(vel_bar[i][1],2)) * dt_mechanical_bar;//energy per voulme
        stress[i][0]=stress_coeff*sum_ac_of_stress[i]; //sigma_xx
      //stress[i][1]=stress_coeff*sum_ad_of_stress[i]; //sigma_xy
      //stress[i][2]=stress_coeff*sum_bc_of_stress[i]; //sigma_yx
      //stress[i][3]=stress_coeff*sum_bd_of_stress[i]; //sigma_yy
      //strain[i]=(disp_bar[i+1][0]-disp_bar[i][0])/dx;//disp_bar[i][0]/coord_bar[i][0];//sum_strain[i]/numfam[i];
    }


 //!Adaptive dynamic relaxation--------------------------------------------------------------------------
  /*  #pragma omp for reduction(+:cn1, cn2) schedule(static)
    for (int i = 0;i<totnode;i++){
        if (vel_barhalfold[i][0]!=0.0) {
            cn1 -= disp_bar[i][0] * disp_bar[i][0] * (pforce_mechanical[i][0] / massvec_mechanical[i][0] - pforceold[i][0] / massvec_mechanical[i][0]) / (dt_mechanical_bar * vel_barhalfold[i][0]);
        }
        //for 2D
      /* if (vel_barhalfold[i][1]!=0) {
            cn1 -= disp_bar[i][1] * disp_bar[i][1] * (pforce_mechanical[i][1] / massvec_mechanical[i][1] - pforceold[i][1] / massvec_mechanical[i][1]) / (dt_mechanical_bar * vel_barhalfold[i][1]);
        }*/
     /*   cn2 +=  disp_bar[i][0] * disp_bar[i][0];
        //for 2D
        //cn2 += disp_bar[i][1] * disp_bar[i][1];
    }

    // Compute damping on a single thread, then synchronize
    #pragma omp single
    {
    if (cn2 != 0.0) {
        if ((cn1 / cn2) > 0.0) cn = 2.0 * std::sqrt(cn1 / cn2);
        else cn = 0.0;
    } else {
        cn = 0.0;
    }
    if (cn > 2.0) cn = 1.9;
    }
    #pragma omp barrier*/
    //-------------------------------------------------
    // compute once per step (exclude Dirichlet nodes)
if (!fixed_damping_coefficient) {
double num = 0.0, den = 0.0;
#pragma omp parallel for reduction(+:num,den)
for (int i = 0; i < totnode; ++i) {
    const double m  = massvec_mechanical[i][0];

    // centered vel_barocity at time n (you have vel_bar^n already)
    const double v  = vel_bar[i][0];

    // acceleration-like change from forces (using current & previous internal force)
    const double a_cur =  pforce_mechanical[i][0] / m;
    const double a_prv =  pforceold[i][0]         / m;
    const double da    = (a_cur - a_prv) / dt_mechanical_bar;

    // use disp_barlacement as weight (like ADR), but you can also use 1.0
    const double u = disp_bar[i][0];

    // Only accumulate where v != 0 to avoid blowups
    if (std::abs(v) > 1e-18) {
        num -= u*u * da / v;   // “stiffness–like” quotient
        den += u*u;
    }
}

// propose cn
double cn_opt = 0.0;
if (den > 0.0 && num > 0.0) cn_opt = 2.0 * std::sqrt(num/den);  // same spirit as ADR

// smoothing & clamps
const double max_cndt = 0.8;                        // keep cn*dt below ~O(1)
const double cn_cap   = max_cndt / dt_mechanical_bar;   // convert to 1/s cap
cn_opt = std::min(cn_opt, cn_cap);
cn = 0.7*cn + 0.3*cn_opt;                 // simple low-pass
}
//!---------------------------------------------------------------------------------------------------

    //maxval=0.0;
    //cn=0.01;//constant damping coefficient
    #pragma omp for schedule(static) \
    reduction(max:maxval) \
    reduction(+:total_energy_of_system,total_kinetic_energy_of_system,total_damping_dissipation_energy_of_system,total_fracture_energy_of_system,total_strain_energy_of_system,total_substrate_energy_of_system)
    for (int i = 0;i<totnode;i++){
            const bool prescribed_boundary_node =
                (i < prescribed_boundary_nodes) ||
                (i >= totnode - prescribed_boundary_nodes);
            const double prescribed_disp_x = epsilonn * coord_bar[i][0];
            const double old_disp_bar = disp_bar[i][0];
            const double m  = massvec_mechanical[i][0];
            // v^n (centered) from half-step vel_barocities
           // vel_bar[i][0] = 0.5 * (vel_barhalfold[i][0] + vel_barhalf[i][0]);
            const double ax = (pforce_mechanical[i][0] + bforce_mechanical[i][0]- cn * vel_bar[i][0]) / m;

        //! Integrate acceleration over time
       /* if (counterADR==0) {
            vel_barhalf[i][0] = (dt_mechanical_bar / massvec_mechanical[i][0]) * (pforce_mechanical[i][0] + bforce_mechanical[i][0]) / 2.0;
            //vel_barhalf[i][1] = (dt_mechanical_bar / massvec_mechanical[i][1]) * (pforce_mechanical[i][1] + bforce_mechanical[i][1]) / 2.0;
        }
        else
            {
            vel_barhalf[i][0] = ((2.0 - cn * dt_mechanical_bar) * vel_barhalfold[i][0] + 2.0 * dt_mechanical_bar / massvec_mechanical[i][0] * (pforce_mechanical[i][0] + bforce_mechanical[i][0])) / (2.0 + cn * dt_mechanical_bar);
            //vel_barhalf[i][1] = ((2.0 - cn * dt_mechanical_bar) * vel_barhalfold[i][1] + 2.0 * dt_mechanical_bar / massvec_mechanical[i][1] * (pforce_mechanical[i][1] + bforce_mechanical[i][1])) / (2.0 + cn * dt_mechanical_bar);
            }*/

            /*if ((loop_counter==0)&& (counterADR==0 || counterADR==1)){//only for the first time increment and first ADR loop
               //   if (abs(vel_barhalf[i][0])>= maxval)
               // {
                    maxval=10.0;//abs(vel_barhalf[i][0])*1.3e6;
               // }
            }*/
            /*if ((loop_counter==0)&& (counterADR<1000)){
                //double relative_change =abs((disp_bar[i][0]-disp_barold[i][0]));
                //if (relative_change>= maxval)
                //{
                    maxval=10.0;//relative_change;
                    //cout<<maxval<<endl;
               // }
            }
            else {
                    maxval=abs((disp_bar[i][0]-disp_barold[i][0])/disp_barold[i][0]);
           // cout <<"stops at "<<counterADR<<endl;
                double relative_change =abs((disp_bar[i][0]-disp_barold[i][0])/disp_barold[i][0]);
                if (relative_change>= maxval)
                {
                    maxval=relative_change;
                }
            }*/

          /*  vel_bar[i][0] = 0.5 * (vel_barhalfold[i][0] + vel_barhalf[i][0]);
            //vel_bar[i][1] = 0.5 * (vel_barhalfold[i][1] + vel_barhalf[i][1]);

            //disp_barold[i][0]=disp_bar[i][0];
            disp_bar[i][0] += vel_barhalf[i][0] * dt_mechanical_bar;
            //disp_bar[i][1] = 0.0;//disp_bar[i][1] + vel_barhalf[i][1] * dt_mechanical_bar;

            vel_barhalfold[i][0] = vel_barhalf[i][0];

            //vel_barhalfold[i][1] = vel_barhalf[i][1];*/
            pforceold[i][0] = pforce_mechanical[i][0];
            //pforceold[i][1] = pforce_mechanical[i][1];


            //!Calculation of vel_barocity of material point i
            //!by integrating the acceleration of material point i
            if (prescribed_boundary_node) {
                vel_bar[i][0] = 0.0;
                disp_bar[i][0] = prescribed_disp_x;
            } else {
                vel_bar[i][0] +=  ax * dt_mechanical_bar;
                //vel_bar[i-1][1] = vel_bar[i-1][1] + acc[i-1][1]*dt_mechanical_bar;

                disp_bar[i][0] +=  vel_bar[i][0] * dt_mechanical_bar;
            }
           // disp_bar[i-1][1] = 0.0;//disp_bar[i-1][1] + vel_bar[i-1][1] * dt_mechanical_bar;

            const double v0 = vel_bar[i][0]/*, v1 = vel_bar[i][1]*/;
            const double speed = std::abs(v0);
            if (std::isfinite(speed) && speed > maxval) maxval = speed;

            kinetic_energy[i] = 0.5 * m * v0*v0;
            damping_dissipation_energy[i] += /*0.5 **/ m * cn * (v0*v0 /*+ v1*v1*/) * dt_mechanical_bar;//energy per voulme
            total_energy_density[i] = strain_energy[i]+damping_dissipation_energy[i]+kinetic_energy[i]+fracture_energy[i]+substrate_energy[i]; //just remember damping_dissipation_energy is history-dependent and increasing (not conservative): it is only energy per node
            total_energy_of_system += total_energy_density[i]*vol_bar;
            total_kinetic_energy_of_system += kinetic_energy[i] *vol_bar;
            total_damping_dissipation_energy_of_system += damping_dissipation_energy[i] *vol_bar;
            total_fracture_energy_of_system +=fracture_energy[i]*vol_bar;
            total_strain_energy_of_system +=strain_energy[i]*vol_bar;
            total_substrate_energy_of_system +=substrate_energy[i]*vol_bar;

            }
    } // end parallel region
    outfile16 << total_energy_of_system << " ";
    outfile18 << total_kinetic_energy_of_system << " ";
    outfile19 << total_damping_dissipation_energy_of_system << " ";
    outfile20 << total_fracture_energy_of_system << " ";
    outfile21 << total_strain_energy_of_system << " ";
    outfile22 << total_substrate_energy_of_system << " ";

   /* for (int i = 0; i < totnode; ++i) {
    total_energy_of_system += total_energy_density[i] * vol_bar;
    }*/

//!To check numerical noise-------------------------------
if (track_numerical_noise) {
    // Write diagnostic data every 100 ADR steps OR if a bond breaks
    bool any_bond_broke = false;
    for (int i = 0; i < totnode && !any_bond_broke; ++i) {
        for (int j = 0; j < numfam[i]; ++j) {
            if (fail[i][j] == 0) { // just broke
                any_bond_broke = true;
                break;
            }
        }
    }

    if (tt_mech % 100 == 0 || any_bond_broke) {
        diag.diag_stress << loop_counter << ' ' << tt_mech;
        diag.diag_pforce << loop_counter << ' ' << tt_mech;
        diag.diag_disp << loop_counter << ' ' << tt_mech;
        diag.diag_stretch << loop_counter << ' ' << tt_mech;

        for (size_t k = 0; k < diag.diagnostic_nodes.size(); ++k) {
            int ni = diag.diagnostic_nodes[k];
            diag.diag_stress << ' ' << stress[ni][0];
            diag.diag_pforce << ' ' << pforce_mechanical[ni][0];
            diag.diag_disp << ' ' << disp_bar[ni][0];

            // Find maximum stretch for this node
            double max_stretch = 0.0;
            for (int j = 0; j < numfam[ni]; ++j) {
                if (fail[ni][j] == 1) {
                    const int cnode = nodefam[pointfam[ni] + j - 1] - 1;
                    const double nx = (coord_bar[cnode][0] + disp_bar[cnode][0]) - (coord_bar[ni][0] + disp_bar[ni][0]);
                    const double ny = (coord_bar[cnode][1] + disp_bar[cnode][1]) - (coord_bar[ni][1] + disp_bar[ni][1]);
                    const double nlength = std::sqrt(nx*nx + ny*ny);
                    const double stretch = (nlength - bond_idist[ni][j]) * bond_inv_idist[ni][j];
                    if (stretch > max_stretch) max_stretch = stretch;
                }
            }
            diag.diag_stretch << ' '  << max_stretch;
        }

        diag.diag_stress << '\n';
        diag.diag_pforce << '\n';
        diag.diag_disp << '\n';
        diag.diag_stretch << '\n';

        // Flush after writing
        diag.diag_stress.flush();
        diag.diag_pforce.flush();
        diag.diag_disp.flush();
        diag.diag_stretch.flush();
    }
  // Write ALL ADR steps for snapshot increments
    if (do_snap) {
        diag.diag_stress_allsteps << tt_mech;
        diag.diag_pforce_allsteps << tt_mech;
        diag.diag_disp_allsteps << tt_mech;
        diag.diag_stretch_allsteps << tt_mech;

        for (size_t k = 0; k < diag.diagnostic_nodes.size(); ++k) {
            int ni = diag.diagnostic_nodes[k];
            diag.diag_stress_allsteps << ' ' << stress[ni][0];
            diag.diag_pforce_allsteps << ' ' << pforce_mechanical[ni][0];
            diag.diag_disp_allsteps << ' ' << disp_bar[ni][0];

            // Find maximum stretch for this node
            double max_stretch = 0.0;
            for (int j = 0; j < numfam[ni]; ++j) {
                if (fail[ni][j] == 1) {
                    const int cnode = nodefam[pointfam[ni] + j - 1] - 1;
                    const double nx = (coord_bar[cnode][0] + disp_bar[cnode][0]) - (coord_bar[ni][0] + disp_bar[ni][0]);
                    const double ny = (coord_bar[cnode][1] + disp_bar[cnode][1]) - (coord_bar[ni][1] + disp_bar[ni][1]);
                    const double nlength = std::sqrt(nx*nx + ny*ny);
                    const double stretch = (nlength - bond_idist[ni][j]) * bond_inv_idist[ni][j];
                    if (stretch > max_stretch) max_stretch = stretch;
                }
            }
            diag.diag_stretch_allsteps << ' ' << max_stretch;
        }

        diag.diag_stress_allsteps << '\n';
        diag.diag_pforce_allsteps << '\n';
        diag.diag_disp_allsteps << '\n';
        diag.diag_stretch_allsteps << '\n';

          // Flush after writing
        diag.diag_stress_allsteps.flush();
        diag.diag_pforce_allsteps.flush();
        diag.diag_disp_allsteps.flush();
        diag.diag_stretch_allsteps.flush();
    }
}
//!End:To check numerical noise---------------------------

        // Nodal strain definition matching ../pery/pery.py:
        // endpoint values use the adjacent element, interior values average neighbors.
        if (totnode > 1) {
            strain[0] = (disp_bar[1][0] - disp_bar[0][0]) / dx_bar;
            for (int i = 1; i < totnode - 1; ++i) {
                const double left_strain = (disp_bar[i][0] - disp_bar[i - 1][0]) / dx_bar;
                const double right_strain = (disp_bar[i + 1][0] - disp_bar[i][0]) / dx_bar;
                strain[i] = 0.5 * (left_strain + right_strain);
            }
            strain[totnode - 1] =
                (disp_bar[totnode - 1][0] - disp_bar[totnode - 2][0]) / dx_bar;
        }

    //!Writing ADR (inner)loop results
    // per-ADR snapshots
// Always buffer every ADR step; will write to disk at end of increment if needed
    {
        std::vector<double> row_stress(ndivx), row_strain(ndivx),
            row_se(ndivx), row_ke(ndivx), row_fe(ndivx),
            row_de(ndivx), row_sub(ndivx), row_tot(ndivx);
        for (int i = 0; i < ndivx; ++i) {
            row_stress[i] = stress[i][0];
            row_strain[i] = strain[i];
            row_se[i]     = strain_energy[i];
            row_ke[i]     = kinetic_energy[i];
            row_fe[i]     = fracture_energy[i];
            row_de[i]     = damping_dissipation_energy[i];
            row_sub[i]    = substrate_energy[i];
            row_tot[i]    = total_energy_density[i];
        }
        buf_stress.push_back(row_stress);
        buf_strain.push_back(row_strain);
        buf_strain_energy.push_back(row_se);
        buf_kinetic.push_back(row_ke);
        buf_fracture.push_back(row_fe);
        buf_damp.push_back(row_de);
        buf_sub.push_back(row_sub);
        buf_total.push_back(row_tot);
    }
// Track stress convergence for the target node
    if (track_stress_convergence) {
        convergence_file << stress[target_node][0] << ' ';
    }


    counterADR++;

    // Convergence check matching ../pery/pery.py: stop when max velocity is small.
    if (!fixed_nt_mechanical && maxval <= tolerance) {
        break; // converged
    }
    //outfile14 <<cn<<" ";
}//!end of ADR (inner)loop--------------------------------------------------------
outfile16 << "\n";
outfile18 << "\n";
outfile19 << "\n";
outfile20 << "\n";
outfile21 << "\n";
outfile22 << "\n";
//!In case of Transient simulation------------------------------------------------
    /* for (int i = ndivx*ndivy_substrate+1;i<totnode+1;i++){
      //!Calculation of acceleration of material point i
        acc[i-1][0] = (pforce_mechanical[i-1][0] + bforce_mechanical[i-1][0]) / density;
       // acc[i-1][1] = (pforce_mechanical[i-1][1] + bforce_mechanical[i-1][1]) / density;

        //!Calculation of vel_barocity of material point i
        //!by integrating the acceleration of material point i
        vel_bar[i-1][0] =  vel_bar[i-1][0] + acc[i-1][0]*dt_mechanical_bar;
        //vel_bar[i-1][1] = vel_bar[i-1][1] + acc[i-1][1]*dt_mechanical_bar;

        disp_bar[i-1][0] = disp_bar[i-1][0] + vel_bar[i-1][0] * dt_mechanical_bar;
       // disp_bar[i-1][1] = 0.0;//disp_bar[i-1][1] + vel_bar[i-1][1] * dt_mechanical_bar;

    }*/
//!--------------------------------------------------------------------
//!Writing last ADR loop results for each increment /only for the sake of efficient post processing!
//   outfile14<<endl;
for (int i=0;i<ndivx;i++){
    outfile2 <<stress[i][0]<<" ";
    outfile3 <<strain[i]<<" ";
    outfile4 <<strain_energy[i]<<" ";
    outfile5 <<kinetic_energy[i]<<" ";
    outfile6 <<fracture_energy[i]<<" ";
    outfile7 <<damping_dissipation_energy[i]<<" ";
    outfile17 <<substrate_energy[i]<<" ";

}
   outfile2 <<'\n';
   outfile3 <<'\n';
   outfile4 <<'\n';
   outfile5 <<'\n';
   outfile6 <<'\n';
   outfile7 <<'\n';
   outfile17 <<'\n';

   {
       std::ofstream conf_file(result_output_dir + "/conf" + std::to_string(loop_counter + 1) + ".dat");
       conf_file << "# position displacement strain damage\n";
       for (int i = 0; i < totnode; ++i) {
           conf_file << coord_bar[i][0] << ' '
                     << disp_bar[i][0] << ' '
                     << strain[i] << ' '
                     << dmg[i] << '\n';
       }
   }

   // --- Early-stop: detect first crack ---
if (stop_after_first_crack && first_crack_increment < 0) {
    for (int i = 0; i < totnode; ++i) {
        if (dmg[i] > 0.0) {
            first_crack_increment = loop_counter;
            early_stop_increment  = loop_counter + 5;
            std::cout << ">>> First crack detected at increment " << loop_counter
                      << ". Will stop at increment " << early_stop_increment << ".\n";
            break;
        }
    }
}

if (do_snap || (stop_after_first_crack && loop_counter == first_crack_increment)) {
    const std::string tag = std::to_string(loop_counter);
    auto flush = [&](const std::string& name, const std::vector<std::vector<double>>& buf) {
        std::ofstream f(name + tag + ".txt", std::ios::trunc);
        for (const auto& row : buf) {
            for (double v : row) f << v << ' ';
            f << '\n';
        }
    };
    flush("stress",                     buf_stress);
    flush("strain",                     buf_strain);
    flush("strain_energy",              buf_strain_energy);
    flush("kinetic_energy",             buf_kinetic);
    flush("fracture_energy",            buf_fracture);
    flush("damping_dissipation_energy", buf_damp);
    flush("substrate_energy",           buf_sub);
    flush("total_energy_density",       buf_total);
}



//!Writing last ADR loop results for ovito files!
  /*  for (int i=1;i<=number_of_write_incr;i++){ //For writing specific increments due to limited storage
         if (counterADR==write_increments[i-1]){
            #pragma omp parallel for
            for (int i=0;i<totnode;i++){
                disp_bar_u_wr[loop_counter][i]=disp_bar[i][0];
                //disp_bar_v_wr[loop_counter][i]=disp_bar[i][1];
                vel_bar_u_wr[loop_counter][i]=vel_bar[i][0];
                //vel_bar_v_wr[loop_counter][i]=vel_bar[i][1];
                stress_xx[loop_counter][i]=stress[i][0];
                //stress_xy[loop_counter][i]=stress[i][1];
                //stress_yx[loop_counter][i]=stress[i][2];
                //stress_yy[loop_counter][i]=stress[i][3];
                strain_energy_wr[loop_counter][i]=strain_energy[i];
                strain_ave[loop_counter][i]=strain[i];
                dmg_wr[loop_counter][i]=dmg[i];
            }
        }
    }*/
//!writing average stress and strain over the interface to plot stress-strain graph
   /* double sum_interface_stress = 0.0, sum_interface_strain = 0.0;
#pragma omp parallel for reduction(+:sum_interface_stress,sum_interface_strain) schedule(static)
for (int i = 0; i < ndivx; ++i) {
    sum_interface_stress += stress[n_substrate + i][0];
    sum_interface_strain += strain[n_substrate + i];
}
outfile1 << (sum_interface_strain/ndivx) << ' ' << (sum_interface_stress/ndivx) << '\n';
*/

//!Ovito result file-the whole results-------------------------------------------------
// Build all lines in parallel, then write once on the main thread
/*{
    std::vector<std::string> ovito_lines(totint);

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < totint; ++i) {
        std::ostringstream ss;
        ss << coord_bar[i][0] << ' ' << coord_bar[i][1] << ' '
           << disp_bar[i][0]  << ' ' << disp_bar[i][1]  << ' '
           << (coord_bar[i][0] + disp_bar[i][0]) << ' '
           << (coord_bar[i][1] + disp_bar[i][1]) << ' '
           << dmg[i] << ' ' << stress[i][0] << ' ' << strain_energy[i] << ' '
           << strain[i] << '\n';
        ovito_lines[i] = ss.str();
    }

    // single threaded write
    for (int i = 0; i < totint; ++i) {
        outfile15 << ovito_lines[i];
    }
}*/
//!Ovito result file-write a section-------------------------------------------------
// Build all lines in parallel, then write once on the main thread
    // Open outfile15 HERE (not at top of loop) so it's always open when do_snap is true
    std::ofstream outfile15(ovito_output_dir + "/for_ovito_" + std::to_string(loop_counter + 1) + ".xyz");
    outfile15 << totint << '\n'
              << "# x y ux uy x_def y_def damage stress strain_energy strain\n";

    std::vector<std::string> ovito_lines(totint);

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < totint; ++i) {
        std::ostringstream ss;
        ss << coord_bar[i][0] << ' ' << 0.0 << ' '
           << disp_bar[i][0]  << ' ' << 0.0  << ' '
           << (coord_bar[i][0] + disp_bar[i][0]) << ' '
           << 0.0 << ' '
           << dmg[i] << ' ' << stress[i][0] << ' ' << strain_energy[i] << ' '
           << strain[i] << '\n';
        ovito_lines[i] = ss.str();
    }

    for (int i = 0; i < totint; ++i) {
        outfile15 << ovito_lines[i];
    }
    // outfile15 closes automatically here (RAII)


    // --- Early-stop: break outer loop if past the threshold ---
    if (stop_after_first_crack && early_stop_increment >= 0 && loop_counter >= early_stop_increment) {
        std::cout << ">>> Early stop triggered at increment " << loop_counter << ".\n";
        break;
    }

        if (track_stress_convergence) {
            convergence_file << '\n';
        }
        // Close diagnostic allsteps files if they were opened for this increment


close_allsteps_files(diag);

    // At the very end of the outer loop body, before the closing }
    if (stop_after_first_crack && first_crack_increment < 0 && loop_counter == effective_num_loop - 1) {
        std::cout << ">>> WARNING: Reached max_extended_loop=" << max_extended_loop
                  << " with no crack detected.\n";
    }
}//! End of outer loop

first_crack_increment_out = first_crack_increment;

outfile16.close();
outfile18.close();
outfile19.close();
outfile20.close();
outfile21.close();
outfile22.close();


//!To check numerical noise---------------------------
if (track_stress_convergence) {
    convergence_file.close();
}
if (track_numerical_noise) {
    if (diag.diag_stress.is_open()) diag.diag_stress.close();
    if (diag.diag_pforce.is_open()) diag.diag_pforce.close();
    if (diag.diag_disp.is_open()) diag.diag_disp.close();
    if (diag.diag_stretch.is_open()) diag.diag_stretch.close();
}
//!End:To check numerical noise---------------------------
}
