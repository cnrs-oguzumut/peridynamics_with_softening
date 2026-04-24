void surface_correction_factors(
    // scalars
    int    totnode,
    double delta_bar,
    double radij_bar,
    double vol_bar,
    double pi,

    double bc_mechanical_film_bar,

    double h_interface_bar,
    double sedload1_mechanical_substrate,
    double sedload1_mechanical_film,
    double sedload2_mechanical_substrate,
    double sedload2_mechanical_film,
    // arrays (pre-allocated)
    std::vector<std::array<double,2>>& coord_bar,    // [totnode][2]
    std::vector<std::array<double,2>>& disp_bar,     // [totnode][2]

    std::vector<int>&     numfam,                // [totnode]
    std::vector<int>&     pointfam,              // [totnode] (1-based offsets)
    int*     nodefam,               // big flat array (uses 1-based indices)
    double** fac,                   // [totnode][maxfam]
    // outputs
    std::vector<std::array<double,2>>& stendens_mechanical,   // [totnode][2]  (x-load, y-load)
    std::vector<std::array<double,2>>&  fncst_mechanical,      // [totnode][2]
    double** scx_mechanical,        // [totnode][maxfam]
    double** scy_mechanical,        // [totnode][maxfam]
    double** scr_mechanical         // [totnode][maxfam]
) {
    //!Determination of surface correction factors using equating
    //!strain energy obtained by classical continuum mechanics and PD

    std::vector<double> geometric_correction(totnode, 1.0);
    for (int i = 0; i < totnode; ++i) {
        const double horizon_left = std::min(delta_bar, coord_bar[i][0]);
        const double horizon_right = std::min(delta_bar, coord_bar[totnode - 1][0] - coord_bar[i][0]);
        const double visible_horizon = horizon_left + horizon_right;
        const double theoretical_horizon = 2.0 * delta_bar;

        if (visible_horizon > 1.0e-12) {
            geometric_correction[i] = theoretical_horizon / visible_horizon;
        }
    }

    const int prescribed_boundary_nodes = 3;
    for (int k = 0; k < prescribed_boundary_nodes && k < totnode; ++k) {
        const int left = k;
        const int right = totnode - 1 - k;
        geometric_correction[left] = 1.0;
        if (right != left) geometric_correction[right] = 1.0;
    }

    //!Surface correction factors in every direction
     for (int i = 0;i<totnode;i++){
        for (int j = 0;j<numfam[i];j++){
                const int cnode = nodefam[pointfam[i] + j - 1];
                const int cj    = cnode - 1;

              //  const double dxij = coord_bar[cj][0] - coord_bar[i][0];
             //   const double dyij = coord_bar[cj][1] - coord_bar[i][1];

                //idist = sqrt(pow((coord_bar[cnode-1][0] - coord_bar[i][0]),2)+pow((coord_bar[cnode-1][1] - coord_bar[i][1]),2));
                //!vol_barume correction
               /* double theta;
                if (std::abs(dyij) < 1e-10) {
                    theta = 0.0;
                } else if (std::abs(dxij) < 1e-10) {
                    theta = 0.5 * pi; // 90 degrees in radians
                } else {
                    theta = std::atan(std::abs(dyij) / std::abs(dxij));
                }*/
                /*if (abs(coord_bar[cnode-1][1]-coord_bar[i][1])<1e-10){
                    double theta=0.0;
                }
                    else if (abs(coord_bar[cnode-1][0]-coord_bar[i][0])<1e-10){
                        double theta=90.0*pi/180.0;
                    }
                        else{
                        double theta=atan(abs(coord_bar[cnode-1][1]-coord_bar[i][1])/abs(coord_bar[cnode-1][0]-coord_bar[i][0]));
                    }*/

                //!Determination of the surface correction between two material points
                 scr_mechanical[i][j] = 0.5 * (geometric_correction[i] + geometric_correction[cnode-1]);
                 fac[i][j] = 1.0;
                 /*scy_mechanical[i][j]  = (fncst_mechanical[i][1] + fncst_mechanical[cnode-1][1]) / 2.0;

                 // anisotropic correction
                const double c = std::cos(theta);
                const double s = std::sin(theta);

                const double denom = ( (c * c) / (scx_mechanical[i][j] * scx_mechanical[i][j]) )
                                   + ( (s * s) / (scy_mechanical[i][j] * scy_mechanical[i][j]) );

                scr_mechanical[i][j] = std::sqrt(1.0 / denom);*/
               //  scr_mechanical[i][j]  = 1.0 / ((pow(cos(theta),2.0) / pow(scx_mechanical[i][j],2.0)) + (pow(sin(theta),2.0) / pow(scy_mechanical[i][j],2.0)));
               //  scr_mechanical[i][j] = sqrt(scr_mechanical[i][j]);
        }
    }
}
