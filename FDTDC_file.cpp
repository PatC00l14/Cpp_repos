#include <iostream>
#include <math.h>
#include <fstream>

using namespace std;



float pulse(float n_step, float dtt){
    float exponent = ((n_step+1)*dtt - 5*16*dtt)/(0.0032*dtt);
    float exponent2 = exponent* exponent;

    float exponent_trial = pow((n_step+1) / (100),2);
    return exp(-exponent_trial);
}


int main(){
    float C = 1.0; //speed of light
    int Nx = 201; //spatial array size
    float dx = 0.001; //spatial step size
    int Nt = 1000; // numner of time iterations
    float dt = dx / (1.0 * C); //time step
    float f = 1 / (dt * 10);
    float dtdx = dt / dx;

    float E_y[Nx] = {0}; //E-field
    float H_z[Nx-1] = {0}; //H-field

    int source_ind = 35;

    float E_yl = 0;
    float E_yr = 0;

    ofstream outfile("ftdt_save.csv");
    if (!outfile.is_open()) {
        cerr << "Failed to open file for writing.\n";
        return 1;
    }

    cout << "FDTD starting"<< endl;

    for (int i =0; i <Nt; i++){
        cout << "Timestep: " << i << endl;

        outfile << E_y[0] << ",";
        for (int j=1; j<Nx-1; j++){
            E_y[j] = E_y[j] - dtdx * (H_z[j]-H_z[j-1]);
            outfile << E_y[j] << ",";
        }
        outfile << E_y[Nx] << ",";
        outfile << " \n";
        outfile << H_z[0] << ",";
        for (int j=0; j<Nx-2; j++){
            H_z[j] = H_z[j] - dtdx * (E_y[j+1] - E_y[j]);
            outfile << H_z[j] << ",";            
        }
        outfile << H_z[Nx-1] << ",";
        outfile << "\n";

        //set 0D source - for 1D simulation
        float pulse0 = pulse(i,dt);
        E_y[source_ind] = E_y[source_ind] + pulse0;

        //Boundary conditions
        E_y[0] =  E_yl + (C*dt - dx)*(E_y[1] - E_y[0]);
        E_yl = E_y[1];
        E_y[Nx] = E_yr + (C*dt - dx)*(E_y[Nx-1] - E_y[Nx]);
        E_yr = E_y[Nx-1];

    }
    // Closing the file
    outfile.close();
    cout << "FDTD finished"<< endl;
    return 0;
}