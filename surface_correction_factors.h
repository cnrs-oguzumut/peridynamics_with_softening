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

    //!Loading 1:constant strain deformation in x direction
    for (int i = 0;i<totnode;i++){
        disp_bar[i][0] = 0.001 * coord_bar[i][0];
        disp_bar[i][1]=0;
    }

    // Virtual substrate displacement (follows same affine motion)
    /*auto get_substrate_disp = [&](double x_coord) -> double {
        return 0.001 * x_coord;
    };*/

    //std::ofstream outfile15;
    //outfile15.open("scr.txt");
    for (int i = 0;i<totnode;i++){
        stendens_mechanical[i][0] = 0.0;

        //sedload1_mechanical:strain energy in x dir obtained by analytical classical continuum mechanics
        double sedload1_mechanical =sedload1_mechanical_film;

        //loop over neighbors
        // ========================================
        // UNIFIED LOOP: Film-Film + Film-Substrate
        // ========================================
        for  (int j = 0;j<numfam[i];j++){
            const int cnode = nodefam[pointfam[i] + j - 1];
            const int cj    = cnode - 1;

            // ==============================
            // 1. FILM-FILM BOND
            // ==============================
            const double dx0 = coord_bar[cj][0] - coord_bar[i][0];
            const double dy0 = coord_bar[cj][1] - coord_bar[i][1];
            const double idist = std::sqrt(dx0 * dx0 + dy0 * dy0);
            //idist = sqrt(pow((coord_bar[cnode-1][0] - coord_bar[i][0]),2)+pow((coord_bar[cnode-1][1] - coord_bar[i][1]),2));
            const double dx1 = (coord_bar[cj][0] + disp_bar[cj][0]) - (coord_bar[i][0] + disp_bar[i][0]);
            const double dy1 = (coord_bar[cj][1] + disp_bar[cj][1]) - (coord_bar[i][1] + disp_bar[i][1]);
            const double nlength = std::sqrt(dx1 * dx1 + dy1 * dy1);
            //nlength = sqrt(pow((coord_bar[cnode-1][0] + disp_bar[cnode-1][0] - coord_bar[i][0] - disp_bar[i][0]),2)+pow((coord_bar[cnode-1][1] + disp_bar[cnode-1][1] - coord_bar[i][1] - disp_bar[i][1]),2));

            //fac determines what fraction of each vol_barume is inside the horizon
            if (idist<delta_bar-radij_bar) {
                fac[i][j] = 1.0;
                }
            else if (idist<delta_bar+radij_bar) {
                fac [i][j]= (delta_bar+radij_bar-idist)/(2.0*radij_bar);
                }
            else{
                fac[i][j] = 0.0;
            }

            //bc_mechanical is PD bond constant: force=bc_mechanical*stretch
            // Film-film bond stiffness
            double bc_mechanical=bc_mechanical_film_bar;


            //strain energy obtained by PD
            const double stretch = (nlength - idist) / idist;

            // Film-film strain energy
            //double se_film = 0.25 * bc_mechanical * (stretch * stretch) * idist * vol_bar * fac[i][j];

            // ==============================
            // 2. FILM-SUBSTRATE INTERFACE ENERGY
            // ==============================
           /* double se_interface = 0.0;

            // Virtual substrate displacements
            double u_sub_i = get_substrate_disp(coord_bar[i][0]);
            double u_sub_j = get_substrate_disp(coord_bar[cj][0]);

            // Initial 2D distance to virtual substrate
            double dx_init_sub = coord_bar[cj][0] - coord_bar[i][0];
            double dist_init_sub = std::sqrt(dx_init_sub * dx_init_sub + h_interface_bar * h_interface_bar);

            // Only compute if within horizon
            if (dist_init_sub <= delta_bar) {
                // Deformed distance
                double dx_def_sub = (coord_bar[cj][0] + u_sub_j) - (coord_bar[i][0] + disp_bar[i][0]);
                double dist_def_sub = std::sqrt(dx_def_sub * dx_def_sub + h_interface_bar * h_interface_bar);

                // Interface stretch
                double stretch_sub = (dist_def_sub - dist_init_sub) / dist_init_sub;

                // Interface stiffness
                double c_interface = bc_mechanical_film_bar * c_int_ratio;

                // Interface strain energy
                se_interface = 0.25 * c_interface * (stretch_sub * stretch_sub) * dist_init_sub * vol_bar;
            }*/

            // ==============================
            // 3. TOTAL STRAIN ENERGY
            // ==============================
            //stendens_mechanical[i][0] += se_film;
            stendens_mechanical[i][0] += 0.25 * bc_mechanical
                                       * (stretch * stretch) * idist * vol_bar * fac[i][j];
            //stendens_mechanical[i][0] += 0.5* 0.5 * bc_mechanical * pow(((nlength - idist) / idist),2) * idist * vol_bar * fac[i][j];
    }
        /*Calculation of surface correction factor in x direction
        by finding the ratio of the analytical strain energy density value
        to the strain energy density value obtained from PD Theory*/
        fncst_mechanical[i][0] = sedload1_mechanical / stendens_mechanical[i][0];
    }

    //!Loading 2:constant strain deformation in y direction
   /* for (int i = 0;i<totnode;i++){
        disp_bar[i][1] = 0.001 * coord_bar[i][1];
        disp_bar[i][0]=0;
    }
    //std::ofstream outfile15;
    //outfile15.open("scr.txt");
    for (int i = 0;i<totnode;i++){
        stendens_mechanical[i][1] = 0.0;
        //sedload2_mechanical:strain energy in y dir obtained by classical continuum mechanics

        double sedload2_mechanical =
            (mat_id[i] == 1) ? sedload2_mechanical_substrate
                             : sedload2_mechanical_film;

        for  (int j = 0;j<numfam[i];j++){
            const int cnode = nodefam[pointfam[i] + j - 1];
            const int cj    = cnode - 1;

            const double dx0 = coord_bar[cj][0] - coord_bar[i][0];
            const double dy0 = coord_bar[cj][1] - coord_bar[i][1];
            const double idist = std::sqrt(dx0 * dx0 + dy0 * dy0);

            const double dx1 = (coord_bar[cj][0] + disp_bar[cj][0]) - (coord_bar[i][0] + disp_bar[i][0]);
            const double dy1 = (coord_bar[cj][1] + disp_bar[cj][1]) - (coord_bar[i][1] + disp_bar[i][1]);
            const double nlength = std::sqrt(dx1 * dx1 + dy1 * dy1);

            double bc_mechanical;
            if (mat_id[i] == 1 && mat_id[cj] == 1) {
                bc_mechanical = bc_mechanical_substrate;
            } else {
                bc_mechanical = bc_mechanical_film_bar;
            }

            const double stretch = (nlength - idist) / idist;
            stendens_mechanical[i][1] += 0.25 * bc_mechanical
                                       * (stretch * stretch) * idist * vol_bar * fac[i][j];

            //stendens_mechanical[i][1] += 0.5* 0.5 * bc_mechanical * pow(((nlength - idist) / idist),2) * idist * vol_bar * fac[i][j];
    }

        fncst_mechanical[i][1] = sedload2_mechanical / stendens_mechanical[i][1];
     }*/
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
                 scr_mechanical[i][j] = (fncst_mechanical[i][0] + fncst_mechanical[cnode-1][0]) / 2.0;
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
