int counter=0;
double pi = 2*acos(0.0);


double ctime_mech=0.0;//Current time
double idist = 0.0;//idist: Initial distance
int nnum = 0;//nnum: Material point number
double theta=0.0;
double nlength  = 0.0;//Length of deformed bond
int cnode = 0;//cnode: Current material point
double stretch=0.0;
double dforce1_mechanical= 0.0;
double dforce2_mechanical= 0.0;
double d_strain_energy=0.0;
double d_fracture_energy=0.0;
int write_inc_index=0;

double sum_interface_stress=0.0;
double sum_interface_strain=0.0;
int index=0;
double bc_mechanical =0;
double ac_of_stress=0.0;
//double bc_of_stress=0.0;
//double ad_of_stress=0.0;
//double bd_of_stress=0.0;
