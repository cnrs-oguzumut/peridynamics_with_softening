//!Ovito result file--------------------------------------------------
int ovito_write_start_increment=1;
//int number_of_wrt_incr_per_file=100;
int number_of_ovito_files=num_loop-ovito_write_start_increment+1;
//int number_of_write_incr=number_of_wrt_incr_per_file*number_of_files+1;
int *ovito_write_increments;
ovito_write_increments=new int [number_of_ovito_files];
for (int i=0;i<number_of_ovito_files;i++){
    ovito_write_increments[i]=0.0;
}
for (int i=0; i<number_of_ovito_files;i++){
    ovito_write_increments[i]=ovito_write_start_increment+(num_loop-ovito_write_start_increment)/(number_of_ovito_files-1)*i;
    cout <<ovito_write_increments[i]<<endl;
}
for (int ff=0;ff<number_of_ovito_files;ff++){
    index=ovito_write_increments[ff];
    std::ofstream outfileff2;
    outfileff2.open("for_ovito_"+ std::to_string(index) +".xyz");
    outfileff2 <<totint<<endl<<endl;

    for (int i=0;i<totint;i++){
        outfileff2 << coord[i][0]<<" "<<coord[i][1]<<" "<<disp_u_wr[index-1][i]<<" "<<disp_v_wr[index-1][i]<<" "<<coord[i][0]+disp_u_wr[index-1][i]<<" "<<coord[i][1]+disp_v_wr[index-1][i]<<" "<< dmg_wr[index-1][i] <<" "<<stress_xx[index-1][i]<<" "<<strain_energy_wr[index-1][i];
        outfileff2 <<endl;
    }
}



