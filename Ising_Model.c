// Code adapted from Mathworks Ising Model and Metropolis Algorithm code, Copyright 2017 The MathWorks, Inc. 
// Based of the adaption by Lewis Mosby 02/03/2020-19/08/2020 for PX428
// uses the GNU Scientific Libary for set up of random number generator

// My aim for this code is not to create a carefully optimized implementation of the metropolis algorithm
// but to create a version in c fast enough to generate large amounts of data at low temperatures within a reasonable timeframe

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h> // to allow for boolean flags to be toggled when running tasks with many different requirements
#include <time.h> // for setting a new seed when run
#include <gsl/gsl_rng.h>
// best libary i could find for well distributed random numbers, creating for large simulations such as this one

gsl_rng *r; // global generator declaration

int index_2_1(int x,int y,int N_spins){// returns the 1D index equivelent of a given 2D index, named since index is an existing inbuilt function
    return (y * N_spins) + x; // returns the 1D array index equivalent from a standard 2d array indexing method
}

double rand_x(int x){ // generates a random number from 0 to x
    double rand = gsl_rng_uniform(r);
    return rand*x;
}

double OLD(int x){ //replaced by a difference tand more uniformly distributed generator
    return ((double)rand())/(RAND_MAX)*x; // creates scaled random numbers, since there is no need for the numbers to be secure rand is sufficient
}

double mean_energy_calc(int* spin, int J, double mu, double H, int N_spins){ // calculates the mean energy // UPDATE
    double T_energy = 0;
    double Mag_energy = 0;
    double M_energy = 0;
    for ( int i = 0; i< (N_spins*N_spins); i++){
        int x = i / N_spins;
        int y = i % N_spins;

        int up = (x + 1)%(N_spins);// the modulo accounts for the periodic boundary conditions
        int down = (x - 1 + N_spins)%(N_spins);
        int left = (y - 1 + N_spins)%(N_spins);
        int right = (y + 1)%(N_spins);

        int up_1 = index_2_1(up, y, N_spins); 
        int down_1 = index_2_1(down, y, N_spins); 
        int left_1 = index_2_1(x, left, N_spins); 
        int right_1 = index_2_1(x, right, N_spins); 

        int sum_of_neighbours = spin[up_1] + spin[down_1] + spin[left_1] + spin[right_1];
        T_energy += (-0.5 * J * spin[i] * sum_of_neighbours) - (spin[i] * mu * H);
    } 
    M_energy = T_energy / (N_spins*N_spins);
    return M_energy;
}

double total_energy_calc(int* spin, int J, double mu, double H, int N_spins){
        double T_energy = 0;
    double Mag_energy = 0;
    double M_energy = 0;
    for ( int i = 0; i< (N_spins*N_spins); i++){
        int x = i / N_spins;
        int y = i % N_spins;

        int up = (x + 1)%(N_spins);// the modulo accounts for the periodic boundary conditions
        int down = (x - 1 + N_spins)%(N_spins);
        int left = (y - 1 + N_spins)%(N_spins);
        int right = (y + 1)%(N_spins);

        int up_1 = index_2_1(up, y, N_spins); 
        int down_1 = index_2_1(down, y, N_spins); 
        int left_1 = index_2_1(x, left, N_spins); 
        int right_1 = index_2_1(x, right, N_spins); 

        int sum_of_neighbours = spin[up_1] + spin[down_1] + spin[left_1] + spin[right_1];
        T_energy += (-0.5 * J * spin[i] * sum_of_neighbours) - (spin[i] * mu * H);
    }
    return T_energy;
}

double mean_square_energy_calc(int* spin, int J, double mu, double H, int N_spins){ // calculates the mean of the energy squared to approximate heat capcity
    double T_energy = 0;
    double Mag_energy = 0;
    double M_energy = 0;
    for ( int i = 0; i< (N_spins*N_spins); i++){
        int x = i / N_spins;
        int y = i % N_spins;

        int up = (x + 1)%(N_spins);// the modulo accounts for the periodic boundary conditions
        int down = (x - 1 + N_spins)%(N_spins);
        int left = (y - 1 + N_spins)%(N_spins);
        int right = (y + 1)%(N_spins);

        int up_1 = index_2_1(up, y, N_spins); 
        int down_1 = index_2_1(down, y, N_spins); 
        int left_1 = index_2_1(x, left, N_spins); 
        int right_1 = index_2_1(x, right, N_spins); 

        int sum_of_neighbours = spin[up_1] + spin[down_1] + spin[left_1] + spin[right_1];
        T_energy += pow((-0.5 * J * spin[i] * sum_of_neighbours) - (spin[i] * mu * H),2);
    } 
    M_energy = T_energy / (N_spins*N_spins);
    return M_energy;
}

double magnetization_calc(int* spin, int N_spins){ // calculates the magnetization //UPDATE
    int sum = 0;
    for ( int i = 0; i<(N_spins*N_spins); i++){
        sum += spin[i];
    }
    double M_mean = (double)sum/(double)(N_spins*N_spins);
    return M_mean;
}

void print_current_spins(int* spin, int N_spins, FILE* output_file_spins){//function to output the current spin matrix to a file
    fprintf(output_file_spins,"%d",spin[0]);
    for (unsigned int i=1; i<(N_spins*N_spins); i++){
        if (i%(N_spins) == 0){
            fprintf(output_file_spins,"\n%d",spin[i]); // resets to a new line every N_spins to create a neat output file
        } else {
            fprintf(output_file_spins,",%d",spin[i]);
        }
    }
    fprintf(output_file_spins,"\n\n"); // first \n ends the line, the second creates a blank line in the file to separate data from different iterations
}

void print_current_emi(double E, double Mag, int interation, FILE* output_file_emi){
    fprintf(output_file_emi,"%lf,%lf,%d\n", E, Mag, interation);
}

void print_current_emt(double E, double Mag, double kT, double kTc, FILE* output_file_emt){
    fprintf(output_file_emt,"%lf,%lf,%lf\n", E, Mag, kT/kTc);
}

void print_current_emht(double E, double Mag, double H, double kT, double kTc, FILE* output_file_emht){
    fprintf(output_file_emht,"%lf,%lf,%lf,%lf\n", E, Mag, H, kT/kTc);
}

void print_current_emm(double E, double Mag, double H, FILE* output_file_emi){
    fprintf(output_file_emi,"%lf,%lf,%lf\n", E, Mag, H);
}

void Metropolis(int* spin, double kT, int J, double mu, double H, int N_iter, int N_spins, int power, FILE* output_file_spins, FILE* output_file_emi){ // implementation of the Metropolis algorithm
    unsigned int h_count = 0;
    unsigned int m_count = 0;
    double current_energy = 0;
    double current_magnetization = 0;

    for (int count = 0; count < N_iter; count++){
        int current_index_1 = trunc(rand_x((N_spins*N_spins) - 1)); // the 1 indicates it is an index for the 1D implementation of the array
        //printf("%d\n", current_index_1);

        int x = current_index_1 / N_spins;
        int y = current_index_1 % N_spins;

        int up = (x + 1)%(N_spins);// the modulo accounts for the periodic boundary conditions
        int down = (x - 1 + N_spins)%(N_spins);
        int left = (y - 1 + N_spins)%(N_spins);
        int right = (y + 1)%(N_spins);

        int up_1 = index_2_1(up, y, N_spins); 
        int down_1 = index_2_1(down, y, N_spins); 
        int left_1 = index_2_1(x, left, N_spins); 
        int right_1 = index_2_1(x, right, N_spins); 

        int sum_of_neighbours = spin[up_1] + spin[down_1] + spin[left_1] + spin[right_1];

        double dE = 2*(( J * spin[current_index_1] * sum_of_neighbours) + (mu * H * spin[current_index_1])); // calculation of the change in energy if it flips

        double Boltzmann_P = exp(-dE / kT);

        if (dE <= 0 || rand_x(1) <= Boltzmann_P){  
            spin[current_index_1] *= -1;
            h_count ++;
        } else {
            m_count ++;
        }

        if (count % (power*power) == 0){
            print_current_spins(spin,N_spins,output_file_spins);
            current_energy = mean_energy_calc(spin, J, mu, H, N_spins);
            current_magnetization = magnetization_calc(spin, N_spins);
            print_current_emi(current_energy, current_magnetization, count, output_file_emi);
        } 
    }
    printf("%d hits, %d rejects\n",h_count, m_count);

}

void Metropolis_T(int* spin, double kT, double kTc, int J, double mu, double H, int N_iter, int N_spins, int power, FILE* output_file_emt){ // implementation of the Metropolis algorithm
    unsigned int h_count = 0;
    unsigned int m_count = 0;
    double current_energy = 0;
    double current_magnetization = 0;

    for (int count = 0; count < N_iter; count++){
        int current_index_1 = trunc(rand_x((N_spins*N_spins) - 1)); // the 1 indicates it is an index for the 1D implementation of the array
        //printf("%d\n", current_index_1);

        int x = current_index_1 / N_spins;
        int y = current_index_1 % N_spins;

        int up = (x + 1)%(N_spins);// the modulo accounts for the periodic boundary conditions
        int down = (x - 1 + N_spins)%(N_spins);
        int left = (y - 1 + N_spins)%(N_spins);
        int right = (y + 1)%(N_spins);

        int up_1 = index_2_1(up, y, N_spins); 
        int down_1 = index_2_1(down, y, N_spins); 
        int left_1 = index_2_1(x, left, N_spins); 
        int right_1 = index_2_1(x, right, N_spins); 

        int sum_of_neighbours = spin[up_1] + spin[down_1] + spin[left_1] + spin[right_1];
        //printf("selected node has spin %d\n",spin[current_index_1]);
        //printf("sum of neighbours is %d\n",sum_of_neighbours);

        double dE = 2 *(( J * spin[current_index_1] * sum_of_neighbours) + (mu * H * spin[current_index_1])); // calculation of the change in energy if it flips

        double Boltzmann_P = exp(-dE / kT);

        if (dE <= 0 || rand_x(1) <= Boltzmann_P){  
            spin[current_index_1] = -spin[current_index_1];
            h_count ++;
        } else {
            m_count ++;
        }
    }
    current_energy = mean_energy_calc(spin, J, mu, H, N_spins);
    current_magnetization = magnetization_calc(spin, N_spins);
    print_current_emt(current_energy, current_magnetization, kT, kTc, output_file_emt);
}

void Metropolis_H(int* spin, double kT, double kTc, int J, double mu, double H, int N_iter, int N_spins, int power, int equi_threshold, FILE* output_file_emht){ // implementation of the Metropolis algorithm
    unsigned int h_count = 0;
    unsigned int m_count = 0;
    double current_energy = 0;
    double current_magnetization = 0;
    double mean_squared_energy = 0;
    double kB = 1.380649 * pow(10,-23);
    double average_energy = 0;
    double average_square_energy = 0;

    for (int count = 0; count < N_iter; count++){
        int current_index_1 = trunc(rand_x((N_spins*N_spins) - 1)); // the 1 indicates it is an index for the 1D implementation of the array
        //printf("%d\n", current_index_1);

        int x = current_index_1 / N_spins;
        int y = current_index_1 % N_spins;

        int up = (x + 1)%(N_spins);// the modulo accounts for the periodic boundary conditions
        int down = (x - 1 + N_spins)%(N_spins);
        int left = (y - 1 + N_spins)%(N_spins);
        int right = (y + 1)%(N_spins);

        int up_1 = index_2_1(up, y, N_spins); 
        int down_1 = index_2_1(down, y, N_spins); 
        int left_1 = index_2_1(x, left, N_spins); 
        int right_1 = index_2_1(x, right, N_spins); 

        int sum_of_neighbours = spin[up_1] + spin[down_1] + spin[left_1] + spin[right_1];
        //printf("selected node has spin %d\n",spin[current_index_1]);
        //printf("sum of neighbours is %d\n",sum_of_neighbours);

        double dE = 2 *(( J * spin[current_index_1] * sum_of_neighbours) + (mu * H * spin[current_index_1])); // calculation of the change in energy if it flips

        double Boltzmann_P = exp(-dE / kT);

        if (dE <= 0 || rand_x(1) <= Boltzmann_P){  
            spin[current_index_1] = -spin[current_index_1];
            h_count ++;
        } else {
            m_count ++;
        }
        if (count > (N_iter - (equi_threshold) -1)){ // can toggle the fraction fo N?iter, its is fractionally how far though the simulation you go before yo believe the simulation to be in equilibrium
            average_energy += total_energy_calc(spin, J, mu, H, N_spins);
            average_square_energy += pow(total_energy_calc(spin, J, mu, H, N_spins),2);
        }    
    }

    current_energy = mean_energy_calc(spin, J, mu, H, N_spins);
    current_magnetization = magnetization_calc(spin, N_spins); // the new thing also needs changing here
    double heat_capcity = (((average_square_energy)/((equi_threshold)*N_spins*N_spins)) - (pow(average_energy,2)/(pow((equi_threshold),2)*N_spins*N_spins)))/pow(kT,2);
    print_current_emht(current_energy, current_magnetization, heat_capcity, kT, kTc, output_file_emht);
}

void Metropolis_M(int* spin, double kT, double kTc, int J, double mu, double H, int N_iter, int N_spins, int power, FILE* output_file_emm){ // implementation of the Metropolis algorithm
    unsigned int h_count = 0;
    unsigned int m_count = 0;
    double current_energy = 0;
    double current_magnetization = 0;

    for (int count = 0; count < N_iter; count++){
        int current_index_1 = trunc(rand_x((N_spins*N_spins) - 1)); // the 1 indicates it is an index for the 1D implementation of the array
        //printf("%d\n", current_index_1);

        int x = current_index_1 / N_spins;
        int y = current_index_1 % N_spins;

        int up = (x + 1)%(N_spins);// the modulo accounts for the periodic boundary conditions
        int down = (x - 1 + N_spins)%(N_spins);
        int left = (y - 1 + N_spins)%(N_spins);
        int right = (y + 1)%(N_spins);

        int up_1 = index_2_1(up, y, N_spins); 
        int down_1 = index_2_1(down, y, N_spins); 
        int left_1 = index_2_1(x, left, N_spins); 
        int right_1 = index_2_1(x, right, N_spins); 

        int sum_of_neighbours = spin[up_1] + spin[down_1] + spin[left_1] + spin[right_1];
        //printf("selected node has spin %d\n",spin[current_index_1]);
        //printf("sum of neighbours is %d\n",sum_of_neighbours);

        double dE = 2 *(( J * spin[current_index_1] * sum_of_neighbours) + (mu * H * spin[current_index_1])); // calculation of the change in energy if it flips

        double Boltzmann_P = exp(-dE / kT);

        if (dE <= 0 || rand_x(1) <= Boltzmann_P){  
            spin[current_index_1] = -spin[current_index_1];
            h_count ++;
        } else {
            m_count ++;
        }
    }
    current_energy = mean_energy_calc(spin, J, mu, H, N_spins);
    current_magnetization = magnetization_calc(spin, N_spins);
    print_current_emm(current_energy, current_magnetization, H, output_file_emm);
}

void task_1(FILE* output_file_spins){

    FILE* output_file_emi = fopen("output_spins_emi.txt","w");
    FILE* output_file_emt = fopen("output_spins_emt.txt","w");

    int J = 1; // Constant, >0 for a ferromagnet
    int N_dim = 2;
    int N_spins =20; // Number of spins simulated along one axis
    float P_up = 0.5; // average probability of a partial being initaliased with up spin
    double H = 0; //magnetic field strength
    double mu = 1;// magnetic moment associated with each spin
    double kTc = 2 * J * N_dim; // Curie temperature
    double kT = 0.1 * kTc;
    int N_temp = 200; //  number of considered temperatures
    double kT_max = 2 * kTc;
    double dkT = (double)kT_max /(double)(N_temp - 1); // difference between temperatures being considered

    int power = 8;
    int N_iter = pow(2,power) * pow(N_spins,N_dim); // this ensures on average each spin is considered 2^pow times

    int* spin = malloc(N_spins*N_spins*sizeof(int)); // allocated memory for the array of spins                                                                                                                                                                                                          

    if (!spin){
        printf("Memory allocation failed\n");
    }

    int sum = 0;
    for (int i = 0;i < (N_spins*N_spins); i++){
        double r = P_up - rand_x(1); // subtracts a random number between 0 and 1 from probability it is initalized with spin up
        spin[i] = copysign(1,r); // copys the sign of the number and put is on a one, generating a disribution of +1 or -1 spins
        sum += spin[i];
    }

    double current_temp = 0;

    printf("Now simulating and saving all relevent data at specified temperature\n");
    printf("for sanity the sum of all elements in spin is %d\n",sum);

    double mean_energy_init = mean_energy_calc(spin,J,mu,H,N_spins);
    printf("The mean energy of the system initially is %lf.\n",mean_energy_init);
    double magnetization_init = magnetization_calc(spin, N_spins);
    printf("The magnetization of the system initially is %lf.\n",magnetization_init);

    Metropolis(spin, kT, J, mu, H, N_iter, N_spins, power, output_file_spins, output_file_emi);

    double mean_energy_end = mean_energy_calc(spin,J,mu,H,N_spins);
    printf("The mean energy of the system by the end is %lf.\n",mean_energy_end);
    double magnetization_end = magnetization_calc(spin, N_spins);
    printf("The magnetization of the system by the end is %lf.\n",magnetization_end);

    printf("Now simulating over many temperatures\n");

        for (int i = 0; i<N_temp; i++){

            for (int i = 0;i < (N_spins*N_spins); i++){
                double r = P_up - rand_x(1); // subtracts a random number between 0 and 1 from probability it is initalized with spin up
                spin[i] = copysign(1,r); // copys the sign of the number and put is on a one, generating a disribution of +1 or -1 spins
            }
            Metropolis_T(spin, current_temp, kTc, J, mu, H, N_iter, N_spins, power, output_file_emt);
            current_temp += dkT;
        }

    free(spin);
    fclose(output_file_spins);
    fclose(output_file_emi);
    gsl_rng_free(r);

    printf("All Done :)\n");
}

void task_2(FILE* output_file_spins){
    FILE* output_file_emi = fopen("output_spins_emi.txt","w");
    FILE* output_file_emht = fopen("output_spins_emht.txt","w");

    int J = 1; // Constant, >0 for a ferromagnet
    int N_dim = 2;
    int N_spins = 50; // Number of spins simulated along one axis
    float P_up = 0.5; // average probability of a partial being initaliased with up spin
    double H = 0; //magnetic field strength
    double mu = 1;// magnetic moment associated with each spin
    //double kTc = 2 * J * N_dim; // Curie temperature
    double kTc = (4 * abs(J))/(acosh(3)); //Onsager L. 1944. Phys. Rev. 65(3 & 4) 117-149.
    double kT = 0.1 * kTc;

    int N_temp = 200; //  number of considered temperatures
    double kT_max = 2 * kTc;
    double dkT = (double)kT_max /(double)(N_temp - 1); // difference between temperatures being considered
    int equi_threshold = 100; // needs careful consideration when running at low powers as the system may not actually reach equilibrium, maybe just keep power high
    // is the number of threshold before the end where you can be pretty certain the system is in equilibrium

    int power = 16;
    int N_iter = pow(2,power) * pow(N_spins,N_dim); // this ensures on average each spin is considered 2^pow times

    int* spin = malloc(N_spins*N_spins*sizeof(int)); // allocated memory for the array of spins                                                                                                                                                                                                          

    if (!spin){
        printf("Memory allocation failed\n");
    }

    int sum = 0;
    for (int i = 0;i < (N_spins*N_spins); i++){
        double r = P_up - rand_x(1); // subtracts a random number between 0 and 1 from probability it is initalized with spin up
        spin[i] = copysign(1,r); // copys the sign of the number and put is on a one, generating a disribution of +1 or -1 spins
        sum += spin[i];
    }

    double current_temp = 0;

    printf("Now simulating and saving all relevent data at specified temperature\n");
    printf("for sanity the sum of all elements in spin is %d\n",sum);

    double mean_energy_init = mean_energy_calc(spin,J,mu,H,N_spins);
    printf("The mean energy of the system initially is %lf.\n",mean_energy_init);
    double magnetization_init = magnetization_calc(spin, N_spins);
    printf("The magnetization of the system initially is %lf.\n",magnetization_init);

    Metropolis(spin, kT, J, mu, H, N_iter, N_spins, power, output_file_spins, output_file_emi);

    double mean_energy_end = mean_energy_calc(spin,J,mu,H,N_spins);
    printf("The mean energy of the system by the end is %lf.\n",mean_energy_end);
    double magnetization_end = magnetization_calc(spin, N_spins);
    printf("The magnetization of the system by the end is %lf.\n",magnetization_end);

    printf("Now simulating over many temperatures\n");

    for (int i = 0; i<N_temp; i++){

        for (int i = 0;i < (N_spins*N_spins); i++){
            double r = P_up - rand_x(1); // subtracts a random number between 0 and 1 from probability it is initalized with spin up
            spin[i] = copysign(1,r); // copys the sign of the number and put is on a one, generating a disribution of +1 or -1 spins
        }
        Metropolis_H(spin, current_temp, kTc, J, mu, H, N_iter, N_spins, power, equi_threshold, output_file_emht);
        current_temp += dkT;
    }

    free(spin);
    fclose(output_file_spins);
    fclose(output_file_emi);
    fclose(output_file_emht);
    gsl_rng_free(r);

    printf("All Done :)\n");
}

void task_3(FILE* output_file_spins){
    FILE* output_file_emi = fopen("output_spins_emi.txt","w");
    FILE* output_file_emht = fopen("output_spins_emht.txt","w");
    FILE* output_file_emm = fopen("ouput_spins_emm.txt","w");

    int J = 1; // Constant, >0 for a ferromagnet
    int N_dim = 2;
    int N_spins = 50; // Number of spins simulated along one axis
    float P_up = 0.5; // average probability of a partial being initaliased with up spin
    double H = 2; //magnetic field strength
    double mu = 1;// magnetic moment associated with each spin
    //double kTc = 2 * J * N_dim; // Curie temperature
    double kTc = (4 * abs(J))/(acosh(3)); //Onsager L. 1944. Phys. Rev. 65(3 & 4) 117-149.
    printf("kTc = %lf\n",kTc);
    double kT = 1 * kTc;
    //double kT = 1.5;

    bool Temp_variation = false; // flag to run over a simulation over many temperatures
    int N_temp = 20; //  number of considered temperatures
    double kT_max = 2 * kTc;
    double dkT = (double)kT_max /(double)(N_temp - 1); // difference between temperatures being considered
    int equi_threshold = 10000; 

    bool Mag_variation = true; // Flag to run a simulation over many values of magnetic field strength
    bool Strictly_increasing = false; // Flag to cinside the effects of both increasing and decreasing the magnetic field
    int N_mag = 1500; //  number of considered magnetic field strengths
    double mag_max = 0.1;
    double dmag = ((double)mag_max*2) /(double)(N_mag - 1); // for the 'two sided' version
    //double dmag = (double)mag_max / (double)(N_mag - 1); // for the strictly increasing magnetic field version

    int power = 4;
    int N_iter = pow(2,power) * pow(N_spins,N_dim); // this ensures on average each spin is considered 2^pow times

    int* spin = malloc(N_spins*N_spins*sizeof(int)); // allocated memory for the array of spins                                                                                                                                                                                                          

    if (!spin){
        printf("Memory allocation failed\n");
    }

    int sum = 0;
    for (int i = 0;i < (N_spins*N_spins); i++){
        double r = P_up - rand_x(1); // subtracts a random number between 0 and 1 from probability it is initalized with spin up
        spin[i] = copysign(1,r); // copys the sign of the number and put is on a one, generating a disribution of +1 or -1 spins
        sum += spin[i];
    }

    double current_temp = 0;
    double current_mag = -mag_max;

    printf("Now simulating and saving all relevent data at specified temperature\n");
    printf("for sanity the sum of all elements in spin is %d\n",sum);

    double mean_energy_init = mean_energy_calc(spin,J,mu,H,N_spins);
    printf("The mean energy of the system initially is %lf.\n",mean_energy_init);
    double magnetization_init = magnetization_calc(spin, N_spins);
    printf("The magnetization of the system initially is %lf.\n",magnetization_init);

    Metropolis(spin, kT, J, mu, H, N_iter, N_spins, power, output_file_spins, output_file_emi);

    double mean_energy_end = mean_energy_calc(spin,J,mu,H,N_spins);
    printf("The mean energy of the system by the end is %lf.\n",mean_energy_end);
    double magnetization_end = magnetization_calc(spin, N_spins);
    printf("The magnetization of the system by the end is %lf.\n",magnetization_end);

    if (Temp_variation){
        printf("Now simulating over many temperatures\n");

        for (int i = 0; i<N_temp; i++){

            for (int i = 0;i < (N_spins*N_spins); i++){
                double r = P_up - rand_x(1); // subtracts a random number between 0 and 1 from probability it is initalized with spin up
                spin[i] = copysign(1,r); // copys the sign of the number and put is on a one, generating a disribution of +1 or -1 spins
            }
            Metropolis_H(spin, current_temp, kTc, J, mu, H, N_iter, N_spins, power, equi_threshold, output_file_emht);
            current_temp += dkT;
        }
    }

    if (Mag_variation){
        printf("Now simulating over many magnetic field strengths\n");

        for (int i = 0; i<N_mag; i++){

            for (int i = 0;i < (N_spins*N_spins); i++){
                double r = P_up - rand_x(1); // subtracts a random number between 0 and 1 from probability it is initalized with spin up
                spin[i] = copysign(1,r); // copys the sign of the number and put is on a one, generating a disribution of +1 or -1 spins
            }
            Metropolis_M(spin, kT, kTc, J, mu, current_mag , N_iter, N_spins, power, output_file_emm);
            current_mag += dmag;
        }

        if (!Strictly_increasing){
            current_mag = mag_max;

            for (int i = 0; i<N_mag; i++){

                for (int i = 0;i < (N_spins*N_spins); i++){
                    double r = P_up - rand_x(1); // subtracts a random number between 0 and 1 from probability it is initalized with spin up
                    spin[i] = copysign(1,r); // copys the sign of the number and put is on a one, generating a disribution of +1 or -1 spins
                }
                Metropolis_M(spin, kT, kTc, J, mu, current_mag , N_iter, N_spins, power, output_file_emm);
                current_mag -= dmag;
            }
        }
    }

    free(spin);
    fclose(output_file_spins);
    fclose(output_file_emi);
    fclose(output_file_emht);
    fclose(output_file_emm);
    gsl_rng_free(r);

    printf("All Done :)\n");

}

void task_4(FILE* output_file_spins){// information about autocorrelation can be found in origin under analysis - signal processing
    // information about magnetic susceptibility can be found by using the linear fit function in origin to find a gradient
    FILE* output_file_emi = fopen("output_spins_emi.txt","w");
    FILE* output_file_emht = fopen("output_spins_emht.txt","w");
    FILE* output_file_emm = fopen("ouput_spins_emm.txt","w");

    int J = 1; // Constant, >0 for a ferromagnet
    int N_dim = 2;
    int N_spins = 50; // Number of spins simulated along one axis
    float P_up = 0.5; // average probability of a partial being initaliased with up spin
    double H = 2; //magnetic field strength
    double mu = 1;// magnetic moment associated with each spin
    //double kTc = 2 * J * N_dim; // Curie temperature
    double kTc = (4 * abs(J))/(acosh(3)); //Onsager L. 1944. Phys. Rev. 65(3 & 4) 117-149.
    double kT = 10 * kTc;
    //double kT = 1.5;

    bool Temp_variation = true; // Flag to run a simulation over many temperatures and outputing, average energy, average magnetization, heat capacity and temperature over constant H
    int N_temp = 20; //  number of considered temperatures
    double kT_max = 2 * kTc;
    double dkT = (double)kT_max /(double)(N_temp - 1); // difference between temperatures being considered
    int equi_threshold = 10000; 

    bool Mag_variation = true; // flag to run a simulation over many magnetic field strengths and outputing : average energy, average magnetization and magnetic field strength over a constant temp
    bool Strictly_increasing = false;
    int N_mag = 50; //  number of considered magnetic field strengths
    double mag_max = 5;
    double dmag = ((double)mag_max*2) /(double)(N_mag - 1); // for the 'two sided' version
    //double dmag = (double)mag_max / (double)(N_mag - 1); // for the strictly increasing magnetic field version

    int power = 4;
    int N_iter = pow(2,power) * pow(N_spins,N_dim); // this ensures on average each spin is considered 2^pow times

    int* spin = malloc(N_spins*N_spins*sizeof(int)); // allocated memory for the array of spins                                                                                                                                                                                                          

    if (!spin){
        printf("Memory allocation failed\n");
    }

    int sum = 0;
    for (int i = 0;i < (N_spins*N_spins); i++){
        double r = P_up - rand_x(1); // subtracts a random number between 0 and 1 from probability it is initalized with spin up
        spin[i] = copysign(1,r); // copys the sign of the number and put is on a one, generating a disribution of +1 or -1 spins
        sum += spin[i];
    }

    double current_temp = 0;
    double current_mag = -mag_max;

    printf("Now simulating and saving all relevent data at specified temperature\n");
    printf("for sanity the sum of all elements in spin is %d\n",sum);

    double mean_energy_init = mean_energy_calc(spin,J,mu,H,N_spins);
    printf("The mean energy of the system initially is %lf.\n",mean_energy_init);
    double magnetization_init = magnetization_calc(spin, N_spins);
    printf("The magnetization of the system initially is %lf.\n",magnetization_init);

    Metropolis(spin, kT, J, mu, H, N_iter, N_spins, power, output_file_spins, output_file_emi);

    double mean_energy_end = mean_energy_calc(spin,J,mu,H,N_spins);
    printf("The mean energy of the system by the end is %lf.\n",mean_energy_end);
    double magnetization_end = magnetization_calc(spin, N_spins);
    printf("The magnetization of the system by the end is %lf.\n",magnetization_end);

    if (Temp_variation){
        printf("Now simulating over many temperatures\n");

        for (int i = 0; i<N_temp; i++){

            for (int i = 0;i < (N_spins*N_spins); i++){
                double r = P_up - rand_x(1); // subtracts a random number between 0 and 1 from probability it is initalized with spin up
                spin[i] = copysign(1,r); // copys the sign of the number and put is on a one, generating a disribution of +1 or -1 spins
            }
            Metropolis_H(spin, current_temp, kTc, J, mu, H, N_iter, N_spins, power, equi_threshold, output_file_emht);
            current_temp += dkT;
        }
    }

    if (Mag_variation){
        printf("Now simulating over many magnetic field strengths\n");

        for (int i = 0; i<N_mag; i++){

            for (int i = 0;i < (N_spins*N_spins); i++){
                double r = P_up - rand_x(1); // subtracts a random number between 0 and 1 from probability it is initalized with spin up
                spin[i] = copysign(1,r); // copys the sign of the number and put is on a one, generating a disribution of +1 or -1 spins
            }
            Metropolis_M(spin, kT, kTc, J, mu, current_mag , N_iter, N_spins, power, output_file_emm);
            current_mag += dmag;
        }

        if (!Strictly_increasing){
            current_mag = mag_max;

            for (int i = 0; i<N_mag; i++){

                for (int i = 0;i < (N_spins*N_spins); i++){
                    double r = P_up - rand_x(1); // subtracts a random number between 0 and 1 from probability it is initalized with spin up
                    spin[i] = copysign(1,r); // copys the sign of the number and put is on a one, generating a disribution of +1 or -1 spins
                }
                Metropolis_M(spin, kT, kTc, J, mu, current_mag , N_iter, N_spins, power, output_file_emm);
                current_mag -= dmag;
            }
        }
    }

    free(spin);
    fclose(output_file_spins);
    fclose(output_file_emi);
    fclose(output_file_emht);
    fclose(output_file_emm);
    gsl_rng_free(r);

    printf("All Done :)\n");

}

int main(){

    FILE* output_file = fopen("output_spins_file.txt","w");
    
    //srand(time(NULL)); // used to ensure the program generates different random number each time it is run, not super secure but is probably sufficient
    const gsl_rng_type *T;

    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, (unsigned long)time(NULL));

    int task_num;

    printf("Please enter the number corresponding to the task you are attempting to simulate.\n");
    scanf("%d", &task_num);

    if (task_num < 1 || task_num > 6){
        printf("This is not a task 1-6");
        return 0;
    } else if (task_num == 1){
        task_1(output_file);
    } else if (task_num == 2){
        task_2(output_file);
    } else if (task_num == 3){
        task_3(output_file);
    } else if (task_num == 4){
        task_4(output_file);
    }
    return 0;
}