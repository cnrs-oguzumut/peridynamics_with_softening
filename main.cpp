#include <iostream>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <vector>
#include <array>
#include <cmath>
#include <stdio.h>
#include <algorithm>
#include <omp.h>
#include <cstddef>
#include <new>
#include <sstream>
using namespace std;
using std::cout;
using std::endl;
#include <iomanip>
using std::setw;
using std::max_element;

double heaviside(double t){
    double heav=0;
    if (t>=0){
        heav=1.0;
    }
    else{
        heav=0.0;
    }
    return heav;
}
#include "Geometry.h"
#include "surface_correction_factors.h"
#include "matrices.h"
#include "only_mechanical.h"
/////////////////////////////////////////////////////
int main()
{

#include "variable_initialization.h"
#include "time_settings.h"
//#include "write_settings.h"


//!Geometrical dimensions and PD discretization parameters
double length=1;
int ndivx=300;
//double dx=length/ndivx;
//cout <<"dx is "<<dx<<endl;
//double delta=3.015*dx;//PD horizon
//double area=1;
//double vol=area * dx;
//double radij = dx / 2.0;//radij: Material point radius
int maxfam=20;// Maximum number of material points inside a horizon of a material point
double maxval=1.0;//for dynamic relaxation loop
//int ndivy_substrate=1;
//int ndivy_film=1;
//double hf=ndivy_film*dx;//height of the film
//int ndivy=ndivy_substrate+ndivy_film;
int totint=ndivx;
//int nbnd_mechanical=3;//Number of divisions in the fixed disp_bar region for imaginary BC
int totnode=ndivx;//totnode: Total number of material points
//int n_substrate=ndivx*ndivy_substrate;
//double substrate_width=ndivy_substrate*dx;
//double width=substrate_width+hf;
double thickness=1;
//double y_film=(width/2.0-hf);




//!Material properties
double pratio=1.0/3.0;// Bond-based PD is limited to 1/3
double rho_substrate = 1;//1200.0;
double rho_film= 1;//8000.0;
double emod_substrate=1;//500e9;//192.0e9;
double emod_film=1;//10e9;
double emod_bar_film=emod_film/(1-pow(pratio,2));
//double density=1200;//8000;//7870.0;//1200.0;
//double cv=472.0;//1466.0;
double mu_substrate=emod_substrate/(2*(1.0+pratio));
double mu_film=emod_film/(2*(1.0+pratio));
//double b_substrate=6*mu_substrate/(pi*pow(delta,4)*thickness);
//double b_film=6*mu_film/(pi*pow(delta,4)*thickness);
//double bc_mechanical_substrate =4*b_substrate*delta;
//double bc_mechanical_film = 2*emod_film/(area*pow(delta,2));// 4*b_film*delta;//;//9*emod/(pi*pow(delta,3)*thickness);Bond constant



double stress_coeff=1.0/2.0;





//!BCs and time discretisation parameters
double epsilonn=0.0;//1.0e-5;//initial strain applied as BC
double target_total_strain = 0.004;//0.004;  // 0.1%
int num_loop=80;//80;//number of (outer loop) time increments
double del_epsilonn=target_total_strain / num_loop;//1.0e-3;//strain increment
double tolerance=1.0e-13;//criteria to stop Dynamic relaxation loop
int max_iteration=100000;
//double force_coeff=1.0e-2;
double cn = 10;
//double appres=10e6;
int counterADR=0;//inner loop time step


//!Dimensionless Params
// Geometric
//double L_bar = 1.0;  // By definition
double dx_bar = length/ndivx;
double delta_bar = 3.015*dx_bar;
double radij_bar = dx_bar / 2.0;
double area_bar=1;
double vol_bar=area_bar * dx_bar;
double m = delta_bar / dx_bar;



int first_crack_increment = -1;

//for transient sims:
double wave_speed=sqrt(emod_bar_film/rho_film);

//double critical_time_step=0.3*dx_bar/wave_speed;


// Time scale
double tau = length / std::sqrt(emod_film / rho_film);
double dt_mechanical_bar=dt_mechanical/tau;
//dt_bar = total_time / (tau * 10000.0); // Assuming nt=10000
//total_time_bar = total_time / tau;
//nt = static_cast<int>(total_time_bar / dt_bar);

// Material - micromodulus
double discrete_correction = m / (m + 1.0);
double bc_mechanical_film_bar = 2*emod_film/(area_bar*pow(delta_bar,2));//*discrete_correction;
double force_coeff_bar=1;////////////?
double lambda2 = 1e3;
double dt_pd = sqrt(dx_bar / (delta_bar*delta_bar));
double dt_lambda = 1.0 / sqrt(lambda2);
double dt_bar = 0.3 * min(dt_pd, dt_lambda);
cout<<"critical time step is "<<dt_pd<<endl;
double h_interface_bar = 1.0 * dx_bar; // Vertical spacing to substrate

double Gc=1e-8;//60e-10;

int degradation_model = 1;
//! Degradation model flag
//! 0 = brittle (original sudden bond break)
//! 1 = linear softening
//! 2 = bilinear softening (trilinear)
//! 3 = Cornelissen softening
//! 4 = exponential softening

// Softening rule
double s_initial = sqrt(3.0*Gc/(emod_film*delta_bar)); //1D

// choose Rd from fracture energy
double Rd = Gc;

// compute s_c depending on softening law
double s_c;

if (degradation_model == 1) {
    // Linear softening
    s_c = 1.5*(s_initial + 2.0 * Rd / (bc_mechanical_film_bar * s_initial));
}
else if (degradation_model == 2) {
    // Bilinear (use same total energy as linear for now)
    s_c = s_initial + 2.0 * Rd / (bc_mechanical_film_bar * s_initial);
}
else if (degradation_model == 3) {
    // Cornelissen
    s_c = s_initial + 5.1361 * Rd / (bc_mechanical_film_bar * s_initial);
}
else if (degradation_model == 4) {
    // Exponential → no finite s_c
    // use large cutoff just for numerical purposes
    s_c = 10.0 * s_initial;
}
else {
    // fallback
    s_c = 4.0 * s_initial;
}

// bilinear parameter
double beta_softening = 0.25;

// define s_k AFTER s_c
double s_k = s_initial + 0.15 * (s_c - s_initial);

//double Rd = 0.5 * bc_mechanical_film_bar * s_initial * (s_c - s_initial);

if (s_c <= s_initial) {
    std::cerr << "ERROR: s_c must be greater than s_initial\n";
    return 1;
}
if (s_k <= s_initial || s_k >= s_c) {
    std::cerr << "ERROR: s_k must satisfy s_initial < s_k < s_c\n";
    return 1;
}
if (Rd <= 0.0) {
    std::cerr << "ERROR: Rd must be positive\n";
    return 1;
}


cout << "s_initial = " << s_initial << endl;
cout << "s_c = " << s_c << endl;
cout << "s_k = " << s_k << endl;
cout << "Rd = " << Rd << endl;


// choose which outer increments to snapshot
std::vector<int> write_snaps = {100000};// will be filled after first crack is detected
int ovito_start_increment = 0;  // start time of writing
double target_x_coord_bar = 0.3;  // Specify the x-coord_barinate of the film node to track

//! Flag to enable/disable numerical noise tracking
bool track_numerical_noise = false;  // Set to false to disable diagnostic output
//! Flag to enable/disable stress convergence tracking
bool track_stress_convergence = true;  // Set to false to disable convergence_of_stress.txt output
//! Flag to use fixed number of ADR steps vs convergence-based stopping
bool fixed_nt_mechanical = false;  // Set to false to use convergence criterion instead of fixed steps
//! Flag to use fixed damping coefficient
bool fixed_damping_coefficient = true;  // Set to false to use convergence criterion instead of fixed steps
//! Flag to stop simulation 5 increments after first crack
bool stop_after_first_crack = true;  // Set to false to run until end of num_loop
//! Flag for force degradation rule


if (degradation_model == 0) cout << "Degradation model: brittle" << endl;
else if (degradation_model == 1) cout << "Degradation model: linear softening" << endl;
else if (degradation_model == 2) cout << "Degradation model: bilinear softening (trilinear)" << endl;
else if (degradation_model == 3) cout << "Degradation model: Cornelissen softening" << endl;
else if (degradation_model == 4) cout << "Degradation model: exponential softening" << endl;
else cout << "Degradation model: unknown" << endl;

//! Defect parameters - to nucleate crack at center
bool introduce_defect = false;              // Enable/disable defect
double defect_x_center = 0.0;              // x-coord_barinate of defect center (center of film)
double defect_y_center = 0.0;              // y-coord_barinate of defect center (at interface)
double defect_radius = dx_bar;              // Radius of defect zone in meters (2mm)
double defect_strength_factor = 0.6;       // Reduce critical stretch to 60% in defect zone

//!parameters needed for calculation of surface correction factors
//double sedload1_mechanical=0.0;
double sedload2_mechanical=0.0;
double sedload1_mechanical_substrate= emod_substrate/(2*(1-pow(pratio,2)))* 1.0e-6;//strain energy density based on classical continuum mechanics
double sedload1_mechanical_film= 0.5*emod_film* 1.0e-6;
double sedload2_mechanical_substrate= emod_substrate/(2*(1-pow(pratio,2)))* 1.0e-6;// 9.0/16.0 * emod * 1.0e-6 ;
double sedload2_mechanical_film= emod_film/(2*(1-pow(pratio,2)))* 1.0e-6;

// declare variables with the same names matrices.h uses
std::vector<std::array<double,2>> stendens_mechanical;
std::vector<std::array<double,2>> fncst_mechanical;
std::vector<double> strain_energy, kinetic_energy, fracture_energy, dissipation_energy, substrate_energy, total_energy;
double **scr_mechanical;
double **bond_dix   = nullptr;
double **bond_diy   = nullptr;
double **bond_idist = nullptr;
double **bond_inv_idist = nullptr;
double **bond_weight = nullptr;
double **bond_s_history = nullptr;
int    **bond_is_sub = nullptr;   // NOTE: int** (mask)

// SoA versions of hot 2D fields
std::vector<double> disp_bar_x, disp_bar_y;
std::vector<double> vel_x, vel_y;
std::vector<double> velhalf_x, velhalf_y;
std::vector<double> velhalfold_x, velhalfold_y;


std::vector<std::array<double,2>> coord_bar;
std::vector<std::array<double,2>> disp_bar, acc;
std::vector<std::array<double,2>> disp_barold;
std::vector<std::array<double,2>> vel;
std::vector<std::array<double,2>> velhalf;
std::vector<std::array<double,2>> velhalfold;
std::vector<std::array<double,2>> pforce_mechanical;
std::vector<std::array<double,2>> bforce_mechanical;
std::vector<std::array<double,2>> pforceold;
std::vector<std::array<double,2>> massvec_mechanical;
std::vector<std::array<double,4>> stress;
std::vector<double> strain;
std::vector<int> numfam, pointfam;
int *nodefam;
std::vector<double> sum_ac_of_stress;
double **fac, **scx_mechanical, **scy_mechanical;
int **fail;
std::vector<double> dmg;

matrices(/*input:*/totnode, maxfam,/*output:*/
    stendens_mechanical, fncst_mechanical, pforce_mechanical,
    strain_energy, kinetic_energy, fracture_energy, dissipation_energy, substrate_energy, total_energy,
    bforce_mechanical, coord_bar, scr_mechanical, bond_s_history,
    pointfam, numfam, nodefam, strain,
    disp_bar, disp_barold, vel, sum_ac_of_stress, fac,
    scx_mechanical, scy_mechanical, acc, stress, fail, dmg,
    massvec_mechanical, velhalf, velhalfold, pforceold
);



build_Geometry(ndivx, length, dx_bar,delta_bar,coord_bar,pointfam,numfam,nodefam,totint,totnode);


surface_correction_factors(totnode, delta_bar, radij_bar, vol_bar, pi,
    bc_mechanical_film_bar,
    h_interface_bar,
    sedload1_mechanical_substrate, sedload1_mechanical_film,
    sedload2_mechanical_substrate, sedload2_mechanical_film,
    coord_bar, disp_bar, numfam, pointfam, nodefam, fac,
    stendens_mechanical, fncst_mechanical, scx_mechanical, scy_mechanical, scr_mechanical);

//for the sake of fast computations
precompute_bond_invariants(
    totnode, maxfam,
    pointfam, numfam, nodefam, coord_bar,
    vol_bar, fac, scr_mechanical,
    bond_dix, bond_diy, bond_idist, bond_inv_idist, bond_weight, bond_is_sub
);


//!Initialization of disp_barlacements and velocities
for (int i = 0;i<totnode;i++){
    vel[i][0]=0.0;
    vel[i][1]=0.0;
    disp_bar[i][0] = 0.0;
    disp_bar[i][1] = 0.0;
}

//!Stable mass vector computation for ADR
/*for (int i = 0;i<totnode;i++){
   //5 is a safety factor

       massvec_mechanical[i][0] = 0.25 * dt_mechanical_bar * dt_mechanical_bar * pi*pow(delta,2)*thickness * bc_mechanical_film / dx * 5.0;//0.25 * dt_mechanical_bar * dt_mechanical_bar * (2.0 * area * delta) * bc_mechanical_film / dx * 5.0;//
       massvec_mechanical[i][1] = 0.25 * dt_mechanical_bar * dt_mechanical_bar * pi*pow(delta,2)*thickness * bc_mechanical_film  / dx * 5.0;//0.25 * dt_mechanical_bar * dt_mechanical_bar * (2.0 * area * delta) * bc_mechanical_film / dx * 5.0;//

}*/

//! Physical lumped mass per node
for (int i = 0; i < totnode; ++i) {
    double m = 1;  // 1=substrate, else film -- kg/m it is density

    massvec_mechanical[i][0] = m;
    massvec_mechanical[i][1] = m;
}



mechanical(
    /* geometry / sizes */
    dx_bar, ndivx, totnode, totint,
    /* time/load */
    num_loop, epsilonn, del_epsilonn, dt_mechanical_bar, nt_mechanical, tolerance, max_iteration,
    /* material/const */
    stress_coeff, s_initial, s_c, s_k, beta_softening, Rd, degradation_model, vol_bar, bc_mechanical_film_bar, cn,
    lambda2, h_interface_bar, delta_bar,
    /* fields (arrays) */
    coord_bar, disp_bar, vel, velhalf, velhalfold, pforce_mechanical, massvec_mechanical,
    bforce_mechanical, pforceold, strain, dmg,
    /* topology */
    numfam, pointfam, nodefam, fail, fac, scr_mechanical, bond_s_history, bond_dix, bond_diy, bond_idist, bond_inv_idist, bond_weight, bond_is_sub,
    /* energies/outputs */
    strain_energy, kinetic_energy, fracture_energy, dissipation_energy, substrate_energy, total_energy,
    sum_ac_of_stress, stress,
    /* snapshot selection */
    write_snaps.data(), static_cast<int>(write_snaps.size()),
    /* ovito start gate */
    ovito_start_increment,
    target_x_coord_bar,
    /* noise tracking flag */
    track_numerical_noise,
     /*  stress convergence tracking flag */
    track_stress_convergence,
    /* fixed nt_mechanical flag */
    fixed_nt_mechanical,
    /* fixed damping_coefficient flag */
    fixed_damping_coefficient,
    /* defect parameters */
    introduce_defect,
    defect_x_center,
    defect_y_center,
    defect_radius,
    defect_strength_factor,
    stop_after_first_crack,
    first_crack_increment
);

// Update write_snaps to the cracking increment
if (first_crack_increment >= 0) {
    write_snaps = {first_crack_increment};
    std::cout << ">>> write_snaps set to cracking increment: " << first_crack_increment << "\n";
} else {
    std::cout << ">>> No crack detected. write_snaps remains empty.\n";
}
///////////////////////////////////////////////////////////////////
//!Post Processing
///////////////////////////////////////////////////////////////////
/*
//!writing average stress and strain over the interface to plot stress-strain graph
std::ofstream outfile1;
outfile1.open("Average_stress_of_interface.txt");
outfile1<<"epsilon"<<" xx-stress"<<endl;
for (int j=0;j<num_loop;j++){
    sum_interface_stress=0.0;
    sum_interface_strain=0.0;
    for (int i=0;i<ndivx;i++){
            sum_interface_stress+=stress_xx[j][n_substrate+i];
            sum_interface_strain+=strain_ave[j][n_substrate+i];
    }
    outfile1<<sum_interface_strain/ndivx<<" "<<sum_interface_stress/ndivx<<endl;
}*/

// free memory
free_matrices(totnode, maxfam,
    stendens_mechanical, fncst_mechanical, pforce_mechanical,
    strain_energy, kinetic_energy, fracture_energy, dissipation_energy, substrate_energy, total_energy,
    bforce_mechanical, coord_bar, scr_mechanical, bond_s_history, bond_dix, bond_diy, bond_idist, bond_inv_idist, bond_weight, bond_is_sub,
    pointfam, numfam, nodefam, strain,
    disp_bar, disp_barold, sum_ac_of_stress, vel, fac,
    scx_mechanical, scy_mechanical, acc, stress, fail, dmg,
    massvec_mechanical, velhalf, velhalfold, pforceold
);


cout<<"Done"<<endl;
printf("\a\a");
    return 0;
}
