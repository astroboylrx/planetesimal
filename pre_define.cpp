//
//  pre_define.cpp
//  planetesimal
//
//  Created by Rixin Li on 9/3/13.
//  Copyright (c) 2013 李 日新. All rights reserved.
//

#include "pre_define.h"


/** member functions for class Parameter **/
int C_Parameter::DeriveParameter()
{
    R_preset = R_inau * AU;
    Mass_index = pow(R_max/R_min, 3./Totalbin);
    Truncate_bin = 3*log(R_trun/R_min)/log(Mass_index);
    
    Totaltime = Totaltime_inKyr * Kyr;
    Timestep = Timestep_inyr * yr;
    Record_times = Totaltime_inKyr*1000/Record_interval_inyr+2;     // the additional one records initial conditions
    
    
    Omega_Kepler = sqrt(Grav*M_solar/pow(R_preset,3.));
    V_Kepler = sqrt(Grav*M_solar/R_preset);
    pMass_tot_i *= M_earth;
    pMass_i *= M_earth;
    pMass_back_i = pMass_tot_i - pMass_i;
    Acc_rate = Acc_rate * M_solar / yr;
    //gMass_tot_i = pMass_tot_i / GPratio_i;
    gMass_tot_i = gMass_d;
    Sigma_gas = (2200 * Fa * pow(R_inau, -3./2)) * 10.0; // BD2008: -1.5 is for MMSN, 0.5 is for observation
    Sigma_p = (33 * Fa * Z_rel * pow(R_inau, -3./2)) * 10.0; // g/cm2 to kg/m2
    T_gas = (120 * pow(R_inau, -3./7));
    Rho_gas = (2.7e-6*pow(R_inau, -39./14));
    Scale_h = (0.022 * R_preset * pow(R_inau, 2./7));
    
    miu = 2.3; // BD2008: using 2.3 is due to consider helium
    C_gas = sqrt(k_B * T_gas / (miu* M_H));
    // Scale_h = C_gas/Omega_Kepler // this definition won't differ much with the former one
    // Rho_gas = Sigma_gas / Scale_h * Exp(-z^2 / 2H^2) / Sqrt(2*Pi)
    sigma_mcs = (2. * pow(10, -19.));
    Lambda_gas = (5.e-3*pow(R_inau, 39./14));
    Eta_sK = (C_gas*C_gas/V_Kepler/V_Kepler);
    V_gas_sK = (V_Kepler*(1-Eta_sK));
    
    Volume_factor_earth = (Rho_water*4.*Pi/3.);
    Volume_factor_water = (Rho_water*4.*Pi/3.);
    Volume_factor_solid = (Rho_solid*4.*Pi/3.);
    
    pAcc_rate = (Acc_rate * Timestep * GPratio_d);
    gVolume = gMass_tot_i / Rho_gas;
    //pVolume = gVolume * sqrt(para->Alpha);
    pVolume = Pi * R_preset * R_preset * Scale_h;
    
    R_in_grav_regime = 10000.0;
    R_in_ep_regime = 1000000.0;
    R_well_coupled = 0.001;
    //mean_H_ratio = 0.01;
    mean_H_ratio = sqrt(Alpha);
    B_totseg = 10;
    
    c_star = 5;
    miu_bar = 0.37;
    countnumber = (long *)malloc(sizeof(long)*100);
    for (int ii = 0; ii < 100; ii++) {
        countnumber[ii] = 0;
    }
    
    return 0;
}

/** member functions for class Delta **/
int C_FileOp::openinput()
{
    if ((fi = fopen("input.txt", "r")) != NULL) {
        return 0;
    }
    if ((fi = fopen("/Users/isaac/Documents/Dropbox/ProgramOS/Xcodestudy/planetesimal_mpi/planetesimal_mpi/input.txt", "r")) != NULL) {
        return 0;
    }
    if ((fi = fopen("input", "r")) != NULL) {
        return 0;
    }
    if ((fi = fopen("INPUT", "r")) != NULL) {
        return 0;
    }
    if ((fi = fopen("input.dat", "r")) != NULL) {
        return 0;
    }
    if ((fi = fopen("input.h", "r")) != NULL) {
        return 0;
    }
    
    printf("Open input file failed. \n");
    exit(1);
}

int C_FileOp::readinput(C_Parameter * para)
{
    char temp[500];
    if (!feof(fop->fi)) {
        
        fgets(temp, 500, fop->fi);
        while (temp[0] == '#') {
            if (!feof(fop->fi)) {
                fgets(temp, 500, fop->fi);
            } else {
                return 1;
            }
        }
        printf("%s\n", temp);
        sscanf(temp, "%d%ld%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf", &para->Order_number, &para->Totalbin, &para->R_min, &para->R_max, &para->Totaltime_inKyr, &para->Timestep_inyr, &para->Record_interval_inyr, &para->v_frag_th, &para->v_disr_th, &para->pMass_tot_i, &para->pMass_i, &para->Rho_solid, &para->Acc_rate, &para->GPratio_i, &para->R_inau, &para->R_trun, &para->Index_init, &para->Q_lr, &para->Alpha, &para->Fa, &para->Z_rel);
        
    } else {
        return 1;
    }
    
    return 0;
}

int C_FileOp::strcatname(char *tail)
{
    strcat(f_output, tail);
    strcat(f_distribution, tail);
    strcat(f_mass, tail);
    strcat(f_coag, tail);
    strcat(f_frag, tail);
    strcat(f_vd, tail);
    
    return 0;
}


int C_FileOp::filename_gen(C_Parameter * para)
/** generate output filename **/
{
    
    /*f_output: file to save initial condition.
     *f_distribution: file to output distributions.
     *f_mass: file to output masses in each component at different time step.
     *f_coag: file to output coagulation contribution in each bin.
     *f_frag: file to output fragmentation contribution in each bin. */
    
    int i;
    char temp[10];
    
    strcpy(f_output, outputpath);
    strcpy(f_distribution, outputpath);
    strcpy(f_mass, outputpath);
    strcpy(f_coag, outputpath);
    strcpy(f_frag, outputpath);
    strcpy(f_vd, outputpath);
    
    i = 0;
    do {
        i++;
    } while(f_output[i] != '.');
    do {
        f_output[i] = '\0';
        f_distribution[i] = '\0';
        f_mass[i] = '\0';
        f_coag[i] = '\0';
        f_frag[i] = '\0';
        f_vd[i] = '\0';
        i++;
    } while(f_output[i] != '\0');
    
    sprintf(temp, "-%d%c", para->Order_number, '\0');
    strcatname(temp);
    
    strcat(f_output, "-output.txt");
    strcat(f_distribution, "-distribution.txt");
    strcat(f_mass, "-mass.txt");
    strcat(f_coag, "-coagulation.txt");
    strcat(f_frag, "-fragmentation.txt");
    strcat(f_vd, "-v_dispersion.txt");
    
    printf("outputfiles will be:\n %s\n %s\n %s\n %s\n %s\n %s\n", f_output, f_distribution, f_mass, f_coag, f_frag, f_vd);
    return 0;
}

int C_FileOp::openoutput()
{
    if (fop->conti == 0) {
        if ((fo = fopen(f_output, "w")) == NULL) {
            printf("Fail to open output file: %s \n", f_output);
            exit(1);
        }
        if ((fd = fopen(f_distribution, "w")) == NULL) {
            printf("Fail to open output file: %s \n", f_distribution);
            exit(1);
        }
        if ((fm = fopen(f_mass, "w")) == NULL) {
            printf("Fail to open output file: %s \n", f_mass);
            exit(1);
        }
        if ((fc = fopen(f_coag, "w")) == NULL) {
            printf("Fail to open output file: %s \n", f_coag);
            exit(1);
        }
        if ((ff = fopen(f_frag, "w")) == NULL) {
            printf("Fail to open output file: %s \n", f_frag);
            exit(1);
        }
        if ((fv = fopen(f_vd, "w")) == NULL) {
            printf("Fail to open output file: %s \n", f_vd);
            exit(1);
        }
    } else {
        if ((fo = fopen(f_output, "a+")) == NULL) {
            printf("Fail to open output file: %s \n", f_output);
            exit(1);
        }
        if ((fd = fopen(f_distribution, "a+")) == NULL) {
            printf("Fail to open output file: %s \n", f_distribution);
            exit(1);
        }
        if ((fm = fopen(f_mass, "a+")) == NULL) {
            printf("Fail to open output file: %s \n", f_mass);
            exit(1);
        }
        if ((fc = fopen(f_coag, "a+")) == NULL) {
            printf("Fail to open output file: %s \n", f_coag);
            exit(1);
        }
        if ((ff = fopen(f_frag, "a+")) == NULL) {
            printf("Fail to open output file: %s \n", f_frag);
            exit(1);
        }
        if ((fv = fopen(f_vd, "w")) == NULL) {
            printf("Fail to open output file: %s \n", f_vd);
            exit(1);
        }
    }
    
    return 0;
}

/*
int C_FileOp::openoutput()
{
    if ((fo = fopen(f_output, "w")) == NULL) {
        printf("Fail to open output file: %s \n", f_output);
        exit(1);
    }
    if ((fd = fopen(f_distribution, "w")) == NULL) {
        printf("Fail to open output file: %s \n", f_distribution);
        exit(1);
    }
    if ((fm = fopen(f_mass, "w")) == NULL) {
        printf("Fail to open output file: %s \n", f_mass);
        exit(1);
    }
    if ((fc = fopen(f_coag, "w")) == NULL) {
        printf("Fail to open output file: %s \n", f_coag);
        exit(1);
    }
    if ((ff = fopen(f_frag, "w")) == NULL) {
        printf("Fail to open output file: %s \n", f_frag);
        exit(1);
    }
    if ((fv = fopen(f_vd, "w")) == NULL) {
        printf("Fail to open output file: %s \n", f_vd);
        exit(1);
    }
    return 0;
}
//*/

int C_FileOp::print_paras()
{
    fprintf(fo, "This file will present the initial parameters we used.\n");
    
    fprintf(fo, "Total mass bins: %ld\n", para->Totalbin);
    fprintf(fo, "Radius of min mass bin: %f\n", para->R_min);
    fprintf(fo, "Radius of max mass bin: %f\n", para->R_max);
    
    fprintf(fo, "Fragmentation threshold velocity: %f\n", para->v_frag_th);
    fprintf(fo, "Disruption threshold velocity: %f\n", para->v_disr_th);
    
    fprintf(fo, "Initial total solid mass: %f M_earth\n", para->pMass_tot_i/M_earth);
    fprintf(fo, "Initial planetesimal mass excluding background dust: %f M_earth\n", para->pMass_i/M_earth);
    fprintf(fo, "Solid materials' density: %f kg/m^3\n", para->Rho_solid);
    fprintf(fo, "Accretion rate: %f M_solar/yr\n", para->Acc_rate * yr / M_solar);
    fprintf(fo, "Initial dust-to-gas ratio: %f\n", para->GPratio_i);
    
    fprintf(fo, "Pre-set radius from star: %fAU\n", para->R_inau);
    fprintf(fo, "Pre-set radius for the largest boulder: %f m\n", para->R_trun);
    fprintf(fo, "Initial index: %f\n", para->Index_init);
    fprintf(fo, "While disruption: Largest remnant mass/Totalmass = Q_lr: %f\n", para->Q_lr);
    fprintf(fo, "The Shakura–Sunyaev alpha parameter: %f\n", para->Alpha);
    
    fprintf(fo, "Mass of center star: %f M_solr\n", M_star / M_solar);
    fprintf(fo, "Radius of center star: %f R_solar\n", R_star / R_solar);
    fprintf(fo, "Disk Mass/M_MMSN = Fa: %f\n", para->Fa);
    fprintf(fo, "Disk Metalicity/Z_MMSN = Z_rel: %f\n", para->Z_rel);
    
    fprintf(fo, "\n");
    fprintf(fo, "Designed evolution time: %f Kyr\n", para->Totaltime_inKyr);
    fprintf(fo, "Timestep: %fyr\n", para->Timestep_inyr);
    fprintf(fo, "Recording interval: %fyr\n", para->Record_interval_inyr);
    
    fprintf(fop->fm, "# This file will record the mass of the following in certain time, in unit of M_earth.\n");
	fprintf(fop->fm, "#time(Kyr)\tplanetesimal\tbackground\tghost_mass\ttotal_mass\tgas_mass\tdust-to-gas ratio\n");
    
    
    return 0;
}

int C_FileOp::dataoutput(long count_column)
/** print data result **/
{
    long i = 0, j = 0;
    
    fprintf(fd, "%e\t%e\t", 0.0, 0.0);
    fprintf(fc, "%e\t%e\t", 0.0, 0.0);
    fprintf(ff, "%e\t%e\t", 0.0, 0.0);
    for (i = 0; i < count_column; i++) {
        fprintf(fd, "%e\t", storage->Evotime[i]);
        fprintf(fc, "%e\t", storage->Evotime[i]);
        fprintf(ff, "%e\t", storage->Evotime[i]);
    }
    fprintf(fd, "\n");
    fprintf(fc, "\n");
    fprintf(ff, "\n");
    fprintf(fd, "%e\t%e\t", 0.0, 0.0);
    fprintf(fc, "%e\t%e\t", 0.0, 0.0);
    fprintf(ff, "%e\t%e\t", 0.0, 0.0);
    for (i = 0; i < count_column; i++) {
        fprintf(fd, "%e\t", storage->GPratio[i]);
        fprintf(fc, "%e\t", storage->GPratio[i]);
        fprintf(ff, "%e\t", storage->GPratio[i]);
    }
    fprintf(fd, "\n");
    fprintf(fc, "\n");
    fprintf(ff, "\n");
    for (i = 0; i < para->Totalbin; i++) {
        fprintf(fd, "%e\t%e\t", r[i], m[i]);
        fprintf(fc, "%e\t%e\t", r[i], m[i]);
        fprintf(ff, "%e\t%e\t", r[i], m[i]);
        
        for (j = 0 ; j < count_column; j++) {
            fprintf(fd, "%e\t", storage->result[j][i]);
            fprintf(fc, "%e\t", storage->coag[j][i]);
            fprintf(ff, "%e\t", storage->frag[j][i]);
        }
        fprintf(fd, "\n");
        fprintf(fc, "\n");
        fprintf(ff, "\n");
        
    }
    
    fprintf(fv, "%e\t%e\t", 0.0, 0.0);
    for (i = 0; i < para->Totalbin; i++) {
        fprintf(fv, "%e\t", r[i]);
    }
    fprintf(fv, "\n");
    fprintf(fv, "%e\t%e\t", 0.0, 0.0);
    for (i = 0; i < para->Totalbin; i++) {
        fprintf(fv, "%e\t", m[i]);
    }
    fprintf(fv, "\n");
    for (i = 0; i < para->Totalbin; i++) {
        fprintf(fv, "%e\t%e\t", r[i], m[i]);
        for (j = 0; j < para->Totalbin; j++) {
            fprintf(fv, "%e\t", v_d_relative[i][j]);
        }
        fprintf(fv, "\n");
    }
    
    
    return 0;
}


char *C_FileOp::readlastline(FILE *readfile)
{
    char *lastline;
    long length;
    fseek(readfile, -1, SEEK_END);
    length = ftell(readfile);
    if (length == 0) {
        printf("The original data file is an empty file.");
        exit(2);
    }
    int keeplooping = 1;
    lastline = (char *)malloc(sizeof(char)*65536);
    while (keeplooping) {
        char ch;
        ch = fgetc(readfile);
        //printf("0 %d\n", (int)ftell(readfile));
        if (ch == '\n') {
            if ((int)ftell(readfile) == length+1) {
                fseek(readfile, -2, SEEK_CUR);
                //printf("1 %d\n", (int)ftell(readfile));
                continue;
            }
            //printf("2 %d\n", (int)ftell(readfile));
            keeplooping = 0;
        } else {
            //printf("3 %d\n", (int)ftell(readfile));
            if((int)ftell(readfile) <= 1 ) {
                printf("This file only has one line of data or text with another empty line.\n");
                fseek(readfile, 0, SEEK_SET);
                //printf("4 %d\n", (int)ftell(readfile));
                keeplooping = 0;
                continue;
            }
            //printf("5 %d\n", (int)ftell(readfile));
            fseek(readfile, -2, SEEK_CUR);
            //printf("6 %d\n", (int)ftell(readfile));
        }
    }
    
    fgets(lastline, 65536, readfile);
    return lastline;
}

string C_FileOp::readlastline_cpp(char *filename)
{
    ifstream fin;
    long length;
    string lastLine;
    fin.open(filename);
    if(fin.is_open()) {
        fin.seekg(-1,ios_base::end);                // go to one spot before the EOF
        length = fin.tellg();
        if(length == -1) {
            cout << filename << "is an empty file." << endl;
            return 0;
        }
        //cout << "length is " << length << endl;
        bool keepLooping = true;
        while(keepLooping) {
            char ch;
            fin.get(ch);                            // Get current byte's data
            if(ch == '\n') {                   // If the data was a newline
                //cout << "1" << (long)fin.tellg() << endl;
                if(fin.tellg() == length+1) {
                    fin.seekg(-2, ios_base::cur);
                    //cout << "2" << (long)fin.tellg() << endl;
                    continue;
                }
                keepLooping = false;                // Stop at the current position.
            } else {                                  // If the data was neither a newline nor at the 0 byte
                //cout << "4" << (long)fin.tellg() << endl;
                fin.seekg(-2,ios_base::cur);        // Move to the front of that data, then to the front of the data before it
                //cout << "5" << (long)fin.tellg() << endl;
                if((long)fin.tellg() < 0) {
                    cout << "This file only has one line of data or text with another empty line." << endl;
                    fin.close();
                    fin.open(filename);
                    fin.seekg(0, ios_base::beg);
                    //cout << "6" <<(long)fin.tellg() << endl;
                    keepLooping = false;
                }
            }
        }
        //cout << "7" <<(long)fin.tellg() << endl;
        getline(fin,lastLine);                      // Read the current line
        //cout << "Result: " << lastLine << '\n';     // Display it
        
        fin.close();
    }
    return lastLine;
    
}



int C_FileOp::closeoutput()
{
    fclose(fd);
    fclose(fo);
    fclose(fm);
    fclose(fc);
    fclose(ff);
    fclose(fv);
    return 0;
}

C_FileOp::~C_FileOp()
{
    fclose(fi);
}

/** member functions for class Delta **/

double *C_Delta::allocate_1d_array(long n)
/* allocate 1d array */
{
    double *temp_array;
    if ((temp_array=(double *)malloc(n*sizeof(double))) == NULL) {
        printf("Fail to allocate space for 1d array.\n");
        exit(2);
    }
    return temp_array;
}

int C_Delta::AllocateSpace(long n)
{
    delta_n = allocate_1d_array(n);
    //delta_Ek_rel = allocate_1d_array(Totalbin);
    
    delta_n_coag = allocate_1d_array(n);
    delta_n_frag = allocate_1d_array(n);
    delta_n_accr = allocate_1d_array(n);
    delta_n_cond = allocate_1d_array(n);
    /*
     delta_n_disr = allocate_1d_array(n);
     delta_n_hitr = allocate_1d_array(n);
     
     */
    
    return 0;
}

int C_Delta::InitializeDelta(long n)
{
    for (long i = 0; i< n; i++) {
        delta_n[i] = 0;
        //delta_Ek_rel[i] = 0;
        delta_n_coag[i] = 0;
        delta_n_frag[i] = 0;
        //delta_n_accr[i] = 0; this is preset
        delta_n_cond[i] = 0;
        /*
         delta_n_disr[i] = 0;
         delta_n_hitr[i] = 0;
         
         */
    }
    
    return 0;
}

int C_Delta::sum_delta_n(long i)
{
    //if (i == 41) {
    //printf("%ld\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", i, n[i], delta_n_coag[i], delta_n_frag[i], delta_n_accr[i], delta_n_cond[i], delta_n_disr[i], delta_n_hitr[i]);
    //}
    delta_n[i] += delta_n_coag[i] + delta_n_frag[i] + delta_n_accr[i] + delta_n_cond[i];
    // + delta_n_disr[i] + delta_n_hitr[i];
    
    //for experiment
    if (n[i] + delta_n[i] < 0) {
        pMass_back -= (n[i] + delta_n[i]) * m[i];
        printf("supply mass from back %e\n", (n[i] + delta_n[i]) * m[i]);
        //n[i] = 0;
    }
    //*/
    
    return 0;
}

C_Delta::~C_Delta()
{
    free(delta_n);
    //free(delta_Ek_rel);
    free(delta_n_accr);
    free(delta_n_coag);
    free(delta_n_frag);
    free(delta_n_cond);
    /*
     free(delta_n_disr);
     free(delta_n_hitr);
     */
}

/** member functions for class Storage **/

double **C_Storage::allocate_2d_array(long n, long m)
/* allocate 2d array */
{
    double **temp_array;
    if ((temp_array=(double **)malloc(n*sizeof(double*))) == NULL) {
        printf("Fail to allocate space for 2d array.\n");
        exit(2);
    }
    for (long i = 0; i < n; i++) {
        if ((temp_array[i]=(double *)malloc(m*sizeof(double))) == NULL) {
            printf("Fail to allocate space for 2d array.\n");
            exit(2);
        }
    }
    return temp_array;
}

int C_Storage::AllocateSpace(long Record_times, long Totalbin)
{
    result = allocate_2d_array(ceil(Record_times), Totalbin);   /* result[0] stores the initial distribution */
    
    coag = allocate_2d_array(ceil(Record_times), Totalbin);     /* others store the separate contributions */
    frag = allocate_2d_array(ceil(Record_times), Totalbin);
    //disr = allocate_2d_array(ceil(Record_times), Totalbin);
    //hitr = allocate_2d_array(ceil(Record_times), Totalbin);
    cond = allocate_2d_array(ceil(Record_times), Totalbin);;
    accr = allocate_2d_array(ceil(Record_times), Totalbin);;
    GPratio = (double *)malloc(Record_times * sizeof(double));
    Evotime = (double *)malloc(Record_times * sizeof(double));
    
    //Ek_evo = allocate_2d_array(ceil(Record_times), Totalbin);
    return 0;
}

int C_Storage::InitializeStorage(long Record_times, long Totalbin)
{
    for (long i = 0; i< Record_times; i++) {
        for (long j = 0; j < Totalbin; j++) {
            result[i][j] = 0;
            coag[i][j] = 0;
            frag[i][j] = 0;
            //disr[i][j] = 0;
            //hitr[i][j] = 0;
            cond[i][j] = 0;
            accr[i][j] = 0;
            //Ek_evo[i][j] = 0;
        }
        GPratio[i] = 0;
        Evotime[i] = 0;
    }
    
    return 0;
}

C_Storage::~C_Storage()
{
    free(result);
    free(coag);
    free(frag);
    //free(disr);
    //free(hitr);
    free(cond);
    free(accr);
    free(GPratio);
    free(Evotime);
    //free(Ek_evo);
}

double Heaviside(double x)
/** define a Heaviside function **/
{
    if (x > 0) return 1;
    if (x < 0) return 0;
    if (x == 0) return 0.5;
    return 0;
}

long findbin(double m_to_locate)
/** define a function to find the bin for a specific mass **/
{
    long temp;
    
    temp = round(log(m_to_locate/m[0])/log(para->Mass_index));
    return temp;
}



