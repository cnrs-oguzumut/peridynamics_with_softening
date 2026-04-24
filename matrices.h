// Function template (helper): allocate rows x cols and fill with init
template <typename T>
T** alloc2D(int rows, int cols, T init) {
    T** a = new T*[rows];
    for (int i = 0; i < rows; ++i) {
        a[i] = new T[cols];
        for (int j = 0; j < cols; ++j) a[i][j] = init;
    }
    return a;
}

inline void matrices(
    int totnode,
    int maxfam,
    // --- outputs (allocated & initialized) ---
    std::vector<std::array<double,2>>& stendens_mechanical,
    std::vector<std::array<double,2>>& fncst_mechanical,
    std::vector<std::array<double,2>>& pforce_mechanical,
    std::vector<double>&  strain_energy,
    std::vector<double>&  kinetic_energy,
    std::vector<double>&  fracture_energy,
    std::vector<double>&  damping_dissipation_energy,
    std::vector<double>& substrate_energy,
    std::vector<double>& total_energy,
    std::vector<std::array<double,2>>& bforce_mechanical,
    std::vector<std::array<double,2>>& coord_bar,
    double**& scr_mechanical,
    double**& bond_s_history,
    std::vector<int>&     pointfam,
    std::vector<int>&     numfam,
    int*&     nodefam,
    std::vector<double>& strain,
    std::vector<std::array<double,2>>& disp,
    std::vector<std::array<double,2>>& dispold,
    std::vector<std::array<double,2>>& vel,

    std::vector<double>&  sum_ac_of_stress,
   // double*&  sum_bc_of_stress,
   // double*&  sum_bd_of_stress,
   // double*&  sum_ad_of_stress,
    double**& fac,
    double**& scx_mechanical,
    double**& scy_mechanical,
    std::vector<std::array<double,2>>& acc,
    std::vector<std::array<double,4>>& stress,
    int**&    fail,
    std::vector<double>& dmg,
    std::vector<std::array<double,2>>& massvec_mechanical,
    std::vector<std::array<double,2>>& velhalf,
    std::vector<std::array<double,2>>& velhalfold,
    std::vector<std::array<double,2>>& pforceold
) {
        // 2D arrays
    stendens_mechanical.assign(totnode, std::array<double,2>{0.0, 0.0});
    fncst_mechanical.assign(totnode, std::array<double,2>{0.0, 0.0});
    pforce_mechanical.assign(totnode, std::array<double,2>{0.0, 0.0});
    bforce_mechanical.assign(totnode, std::array<double,2>{0.0, 0.0});
    coord_bar.assign(totnode, std::array<double,2>{0.0, 0.0});
    scr_mechanical      = alloc2D<double>(totnode, maxfam, 0.0);
    bond_s_history = alloc2D<double>(totnode, maxfam, 0.0);
    disp.assign(totnode, std::array<double,2>{0.0, 0.0});
    dispold.assign(totnode, std::array<double,2>{0.0, 0.0});
    vel.assign(totnode, std::array<double,2>{0.0, 0.0});
    fac                 = alloc2D<double>(totnode, maxfam, 0.0);
    scx_mechanical      = alloc2D<double>(totnode, maxfam, 0.0);
    scy_mechanical      = alloc2D<double>(totnode, maxfam, 0.0);
    acc.assign(totnode, std::array<double,2>{0.0, 0.0});
    stress.assign(totnode, std::array<double,4>{0.0, 0.0});
    massvec_mechanical.assign(totnode, std::array<double,2>{0.0, 0.0});
    velhalf.assign(totnode, std::array<double,2>{0.0, 0.0});
    velhalfold.assign(totnode, std::array<double,2>{0.0, 0.0});
    pforceold.assign(totnode, std::array<double,2>{0.0, 0.0});
    fail                = alloc2D<int>(totnode, maxfam, 1);



    // 1D arrays
    strain_energy.assign(totnode, 0.0);         // zero-initialized
    kinetic_energy.assign(totnode, 0.0);
    fracture_energy.assign(totnode, 0.0);
    damping_dissipation_energy.assign(totnode, 0.0);
    substrate_energy.assign(totnode, 0.0);
    total_energy.assign(totnode, 0.0);
    sum_ac_of_stress.assign(totnode, 0.0);
    //sum_ad_of_stress   = new double[totnode]();
    //sum_bc_of_stress   = new double[totnode]();
    //sum_bd_of_stress   = new double[totnode]();
    strain.assign(totnode, 0.0);
    dmg.assign(totnode, 0.0);

    pointfam.assign(totnode, 0);   // zeros
    numfam.assign(totnode, 0);   // zeros
    nodefam            = new int[totnode * maxfam]();


};
//!================================================================================
// for the sake of fast computations
inline void precompute_bond_invariants(
    int totnode,
    int maxfam,                                  // <— new
    const std::vector<int>& pointfam,
    const std::vector<int>& numfam,
    const int* nodefam,
    const std::vector<std::array<double,2>>& coord_bar,
    double vol,
    double** fac,
    double** scr_mechanical,
    // outputs
    double**& bond_dix,
    double**& bond_diy,
    double**& bond_idist,
    double**& bond_inv_idist,
    double**& bond_weight,
    int**&    bond_is_sub
) {
    // allocate with maxfam columns
    bond_dix       = alloc2D<double>(totnode, maxfam, 0.0);
    bond_diy       = alloc2D<double>(totnode, maxfam, 0.0);
    bond_idist     = alloc2D<double>(totnode, maxfam, 0.0);
    bond_inv_idist = alloc2D<double>(totnode, maxfam, 0.0);
    bond_weight    = alloc2D<double>(totnode, maxfam, 0.0);
    bond_is_sub    = alloc2D<int>   (totnode, maxfam, 0);

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < totnode; ++i) {
        const auto& ci = coord_bar[i];
        const int nf   = numfam[i];

        // optional sanity guard (helpful during debugging)
        // if (nf > maxfam) { fprintf(stderr,"numfam[%d]=%d > maxfam=%d\n", i, nf, maxfam); abort(); }

        const int base = pointfam[i];
        for (int j = 0; j < nf; ++j) {
            const int cnode = nodefam[base + j - 1] - 1;
            const auto& cj  = coord_bar[cnode];

            const double dix = cj[0] - ci[0];
            const double diy = cj[1] - ci[1];
            const double r2  = dix*dix;// + diy*diy;
            const double r   = std::sqrt(r2);
            const double inv = (r > 0.0 ? 1.0/r : 0.0);

            bond_dix[i][j]       = dix;
            bond_diy[i][j]       = diy;
            bond_idist[i][j]     = r;
            bond_inv_idist[i][j] = inv;
            bond_weight[i][j]    = vol * scr_mechanical[i][j] * fac[i][j];
            //bond_is_sub[i][j]    = ((cj[1] - y_film) <= 0.0) ? 1 : 0;
        }
    }
}
//!================================================================================
//! Free memory
// helpers
template <typename T>
void free2D(T**& a, int rows) {
    if (!a) return;
    for (int i = 0; i < rows; ++i) delete[] a[i];
    delete[] a;
    a = nullptr;
}

// free everything allocated by matrices(...)
inline void free_matrices(
    int totnode,
    int maxfam,
    // same pointer set you passed out of matrices(...)
    std::vector<std::array<double,2>>& stendens_mechanical,
    std::vector<std::array<double,2>>& fncst_mechanical,
    std::vector<std::array<double,2>>& pforce_mechanical,
    std::vector<double>&  strain_energy,
    std::vector<double>& kinetic_energy,
    std::vector<double>&  fracture_energy,
    std::vector<double>& damping_dissipation_energy,
    std::vector<double>& substrate_energy,
    std::vector<double>& total_energy,
    std::vector<std::array<double,2>>& bforce_mechanical,
    std::vector<std::array<double,2>>& coord_bar,
    double**& scr_mechanical,
    double**& bond_s_history,
    double**& bond_dix,
    double**& bond_diy,
    double**& bond_idist,
    double**& bond_inv_idist,
    double**& bond_weight,
    int**& bond_is_sub,

    std::vector<int>&      pointfam,
    std::vector<int>&     numfam,
    int*&     nodefam,
    std::vector<double>& strain,
    std::vector<std::array<double,2>>& disp,
    std::vector<std::array<double,2>>& dispold,
    std::vector<double>&  sum_ac_of_stress,
    std::vector<std::array<double,2>>& vel,
    double**& fac,
    double**& scx_mechanical,
    double**& scy_mechanical,
    std::vector<std::array<double,2>>& acc,
    std::vector<std::array<double,4>>& stress,
    int**&    fail,
    std::vector<double>&  dmg,
    std::vector<std::array<double,2>>& massvec_mechanical,
    std::vector<std::array<double,2>>& velhalf,
    std::vector<std::array<double,2>>& velhalfold,
    std::vector<std::array<double,2>>& pforceold
) {
    // 2D double
    stendens_mechanical.clear();
    stendens_mechanical.shrink_to_fit();
    fncst_mechanical.clear();
    fncst_mechanical.shrink_to_fit();
    pforce_mechanical.clear();
    pforce_mechanical.shrink_to_fit();
    bforce_mechanical.clear();
    bforce_mechanical.shrink_to_fit();
    coord_bar.clear();
    coord_bar.shrink_to_fit();
    free2D(scr_mechanical,      totnode); // rows=totnode, each row size=maxfam
    free2D(bond_s_history, totnode);
    disp.clear();
    disp.shrink_to_fit();
    dispold.clear();
    dispold.shrink_to_fit();
    vel.clear();
    vel.shrink_to_fit();
    free2D(fac,                 totnode);
    free2D(scx_mechanical,      totnode);
    free2D(scy_mechanical,      totnode);
    acc.clear();
    acc.shrink_to_fit();
    stress.clear();
    stress.shrink_to_fit();
    massvec_mechanical.clear();
    massvec_mechanical.shrink_to_fit();
    velhalf.clear();
    velhalf.shrink_to_fit();
    velhalfold.clear();
    velhalfold.shrink_to_fit();
    pforceold.clear();
    pforceold.shrink_to_fit();
    free2D(bond_dix,       totnode);
    free2D(bond_diy,       totnode);
    free2D(bond_idist,     totnode);
    free2D(bond_inv_idist, totnode);
    free2D(bond_weight,    totnode);
    free2D(bond_is_sub,    totnode);

    // 2D int
    free2D(fail, totnode);

    // 1D arrays
    strain_energy.clear();
    strain_energy.shrink_to_fit();
    kinetic_energy.clear();
    kinetic_energy.shrink_to_fit();
    fracture_energy.clear();
    fracture_energy.shrink_to_fit();
    damping_dissipation_energy.clear();
    damping_dissipation_energy.shrink_to_fit();
    substrate_energy.clear();
    substrate_energy.shrink_to_fit();
    total_energy.clear();
    total_energy.shrink_to_fit();
    sum_ac_of_stress.clear();
    sum_ac_of_stress.shrink_to_fit();
    strain.clear();
    strain.shrink_to_fit();

    dmg.clear();
    dmg.shrink_to_fit();

    pointfam.clear();
    pointfam.shrink_to_fit();
    numfam.clear();
    numfam.shrink_to_fit();
    delete[] nodefam;   nodefam = nullptr; // size = totnode * maxfam

}


