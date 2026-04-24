void build_Geometry(
    int ndivx,
    double length,
    double dx_bar,
    double delta_bar,
    // outputs (already allocated by matrices.h)
    std::vector<std::array<double,2>>& coord_bar,

    std::vector<int>&     pointfam,
    std::vector<int>&     numfam,
    int*&     nodefam,
    // outputs (scalars)
    int     totint,
    int     totnode
) {
    // Derived dims to mirror main.cpp
    // ── Precompute invariants
//    const double width  = ndivy * dx;
    const double delta2_bar = delta_bar * delta_bar;
    const double x0     = 0.0;
//    const double y0     = -0.5 * width  + 0.5 * dx;

    //!Material points------------------------
    //internal points
    //Generate grid coords + material IDs
    int nnum = 0;
//    for (int i = 0; i < ndivy; ++i) {
//        const double yi = y0 + i * dx;
        for (int j = 0; j < ndivx; ++j) {
            const double xi = x0 + j * dx_bar;
            ++nnum;                         // 1-based counter
            coord_bar[nnum - 1][0] = xi;
            coord_bar[nnum - 1][1] = 0.0;

            // match your intent but fix off-by-one index
          //  mat_id[nnum - 1] = (nnum < n_substrate) ? 1 : 2;
        }
   // }

    totint=nnum;

    //!bottom boundary
    /*for  (int i = 1;i<nbnd_mechanical+1;i++){
            for (int j=1;j<ndivx+1;j++){
                nnum = nnum + 1;
                coord[nnum-1][0] = -1.0 /2.0 * length + (dx / 2.0) + (j - 1) * dx;
                coord[nnum-1][1] = -1.0 /2.0 * width - (dx / 2.0) - (i - 1) * dx;
            }
    }
    int totbottom = nnum;*/

    //!top boundary
    /*for  (int i = 1;i<nbnd_mechanical+1;i++){
            for (int j=1;j<ndivx+1;j++){
                nnum = nnum + 1;
                coord[nnum-1][0] = -1.0 /2.0 * length + (dx / 2.0) + (j - 1) * dx;
                coord[nnum-1][1] = 1.0 /2.0 * width + (dx / 2.0) + (i - 1) * dx;
            }
    }
    int tottop = nnum;*/

    //!left boundary
    /*for (int i=0;i<ndivy;i++){
        for (int j=0; j<nbnd_mechanical;j++){
            nnum = nnum + 1;
            coord[nnum-1][0]=-1.0/2*length-(dx / 2.0)-2*dx + j  * dx;
            coord[nnum-1][1]=-1.0/2*width+(dx / 2.0) + i  * dx;
        }
    }
    int totleft = nnum;*/

    //!right boundary
    /*for (int i=0;i<ndivy;i++){
        for (int j=0; j<nbnd_mechanical;j++){
                nnum = nnum + 1;
        coord[nnum-1][0]=(1.0/2.0)*length+(dx / 2.0) + j  * dx;
        coord[nnum-1][1]=-1.0/2*width+(dx / 2.0) + i  * dx;
        }
    }
    int totright = nnum;*/

    totnode=totint;
    //---------------------------------------------------------------------------
    //!Determination of material points inside the horizon of each material point

    for  (int i = 0;i<totnode;i++){
        if (i==0){
            pointfam[i] = 1;
        }else{
            pointfam[i] = pointfam[i-1]+ numfam[i-1];
        }
        for  (int j = 0; j<totnode;j++){
                if (i == j) continue;

            const double dxij_bar = coord_bar[j][0] - coord_bar[i][0];
            const double dyij_bar = coord_bar[j][1] - coord_bar[i][1];
            const double idist = std::sqrt(dxij_bar * dxij_bar + dyij_bar * dyij_bar);

          // double idist = sqrt(pow((coord[j][0] - coord[i][0]),2)+pow((coord[j][1] - coord[i][1]),2));
            if (idist <= delta_bar + 1.0e-12) {
                ++numfam[i];
                nodefam[pointfam[i] + numfam[i] - 2] = j + 1;
            }
        }
    }
}
