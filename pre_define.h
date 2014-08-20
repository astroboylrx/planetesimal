//
//  pre_define.h
//  planetesimal
//
//  Created by Rixin Li on 9/3/13.
//  Copyright (c) 2013 李 日新. All rights reserved.
//

#ifndef __planetesimal__pre_define__
#define __planetesimal__pre_define__

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <mpi.h>
using namespace std;


//#include book: Astrophysics of Planet Formation by PHILIP J. ARMITAGE, hereafter A's Book
//#include paper: Chiang, Youdin 2010, hereafter CY2010
//                F. Brauer et al. 2008: Particle growth around the snow line, hereafter B2008
//                Brauer, F., Dullemond, C. P. & Henning, T. 2008, hereafter BD2008
//                Leinhardt, Z. M. & Stewart, S. T. 2012,
//                Stewart, S. T. & Leinhardt, Z. M. 2012, hereafter LS2012

#define Asbook_B2008_settling
//#define BD2008_settling

/***** astrophysical constants: *****/

/** fundamental **/
#define Pi              (double)(3.141592653)
#define Grav            (double)(6.6738e-11)        /****the constant of gravity****/
#define k_B             (double)(1.3807e-23)        /****the Boltzman constant****/

/** time **/
#define Day             (double)(8.6400e04)         /****s****/
#define yr              (double)(3.1557e07)         /****s****/ /* for an average Gregorian year */
//#define yr            (double)(3.1536e07)         /****s****/ /* for a common year */
#define Kyr             (double)(1.0000e03*yr)
#define Myr             (double)(1.0000e06*yr)

/** length **/
#define AU              (double)(1.4960e11)         /****m****/

/** solar **/
#define M_solar         (double)(1.9885e30)         /****kg***/ /* http://nssdc.gsfc.nasa.gov/planetary/factsheet/sunfact.html */
#define R_solar         (double)(6.9550e08)         /****m****/

/** terrestrial **/
#define M_earth         (double)(5.9722e24)         /****kg***/
#define R_earth         (double)(6.3710e06)         /****m****/
#define Rho_earth       (double)(5.5515e03)         /*kg/m^3**/
#define Rho_water       (double)(1.0000e03)         /*kg/m^3**/

/** Gas and Hydrogen **/
#define R_H             (double)(1.0000e-10)        /****m****/
#define M_H             (double)(1.6736e-27)        /****kg***/


/***** initial values: *****/

// expressions below are from YS2010
#define M_star          (M_solar)
#define R_star          (1.7*R_solar)
#define T_eff           (4350)
#define gMass_d         (100*M_earth)        // default total gas mass
#define GPratio_d       (0.005)               // default dust-to-gas ratio

/***** global classes: *****/
class C_Parameter {
public:
    /*************************/
    /***** input values: *****/
    /*************************/

    int Order_number;                       // the order if run many times
    
    long Totalbin;
    double R_min;                           // radius for the first mass bin
    double R_max;                           // radius for the last mass bin
    
    double Totaltime_inKyr;
    double Timestep_inyr;
    double Record_interval_inyr;            // time interval of output record;
    
    double v_frag_th;                       // unit: m/s, this is critical fragmentation velocity
    double v_disr_th;                       // unit: m/s, this is critical disruption velocity
    
    double pMass_tot_i;                     // initial total mass of particles and plantesimals
    double pMass_i;                         // initial mass of plantesimals excluding background dust
    double Rho_solid;                       // kg/m^3

    double Acc_rate;
    double GPratio_i;                       // if GPratio_i isn't GPratio_d, then we'll assume the total gas mass
                                            // is gMass_d, and then calculate the corresponding pMass_tot_i
                                            // pMass_i will be set according to the given ratio of (pMass_i/pMass_tot_i)

    double R_inau;
    double R_trun;                          // truncate radius for the largest planetesimals existing initially
    double Index_init;                      // parameter for initial mass distribution and velocity damping factor (relate to inertia): index number
    double Q_lr;                            // the largest post-collision remnant mass ratio relative to the total mass in non-gravity regime
    double Alpha;                           // the Shakura–Sunyaev alpha parameter

    // below from YS2010
    double Fa;                              // how much total mass the disk has relative to MMSN
    double Z_rel;                           // how metal rich the disk is compared with a gas of solar abundances
    // above from YS2010
    
    /************************************/
    /***** derived or other values: *****/
    /************************************/
    double R_preset;                        // the pre-set radius where the gas density culminates near inner boundary
    long Truncate_bin;                      // parameter for initial mass distribution: truncation bin number;
    double Mass_index;                      // index for mass logarithmic distribution
    double Totaltime;
    double Timestep;
    double Omega_Kepler;
    double V_Kepler;
    double pMass_back_i;                    // initial background mass, consists of particles smaller than 0.1m
    double gMass_tot_i;
    // below from YS2010
    double Sigma_gas;                       // Sigma denotes surface density
    double Sigma_p;                         // p denotes particles
    double T_gas;
    double Rho_gas;
    double Scale_h;                         // scale_hight
    double miu;                             // mean molecular weight in unit of M_H
    double C_gas;                           // this one is from A's book
    double sigma_mcs;                       // constant molecular cross section of Hydrogen
    double Lambda_gas;                      // mean free path of gas molecular
    // if we need to consider sub-Kepler velocity of gas
    double Eta_sK;                          // parameter for the sub-Kepler velocity coefficient
    double V_gas_sK;
    // above from YS2010
    double Record_times;
    double Volume_factor_water;
    double Volume_factor_earth;
    double Volume_factor_solid;
    double pAcc_rate;
    double pAcc_intoback;
    double gVolume;
    double pVolume;
    double R_well_coupled;                  // unit: m, used in calculating velocity disperison
    double R_in_grav_regime;                // unit: m, used in calculating velocity disperison
    double R_in_ep_regime;                  // unit: m, used in calculating velocity dispersion, ep means energy equipartition
    double mean_H_ratio;                    // average scale height ratio of gas and dust
    // below from LS2012
    double B_totseg;                        // how many segments to divide
    double c_star;                          // value of the equal-mass head-on disruption crition
    double miu_bar;
    // above from LS2012
    long *countnumber;
    /***************************/
    /**** member functions: ****/
    /***************************/
    C_Parameter();
    ~C_Parameter();
    
    int InitializeParameter();
    int DeriveParameter();
    
private:

    
}; // for all the paramters

class C_FileOp {
public:

    FILE *fi;
    
    int conti;
    char outputpath[250], cont[10];
    char f_output[250], f_distribution[250], f_mass[250], f_coag[250], f_frag[250], f_vd[250];
    FILE *fo, *fd, *fm, *fc, *ff, *fv;
    
    int openinput();
    int readinput(C_Parameter * para);
    
    int strcatname(char *tail);
    /** generate output filename **/
    int filename_gen(C_Parameter * para);
    
    int openoutput();
    
    /** print initial parameters **/
    int print_paras();
    
    /** print data result **/
    int dataoutput(long count_column);
    
    /** read the last line of file **/
    char *readlastline(FILE *readfile);
    std::string readlastline_cpp(char *filename);
    
    /** read the data file into 2-d array **/
    int read_2d_array(FILE *readfile);
    
    int closeoutput();
    ~C_FileOp();
    
}; // for all file operation




class C_Delta {
public:
    double *delta_n;
    //double *delta_Ek_rel;
    double *delta_n_coag, *delta_n_frag, *delta_n_accr;
    //double *delta_n_disr, *delta_n_hitr;
    double *delta_n_cond;
    
    C_Delta(long n);
    ~C_Delta();
    
    int AllocateSpace(long n);
    int InitializeDelta(long n);
    double *allocate_1d_array(long n);
    
    int sum_delta_n(long i);  // calculate the total delta_n
    
private:
    
}; // record the changing things




class C_Storage {
public:
    double **result;
    double **coag, **frag, **cond, **accr; //**disr, **hitr,
    double *GPratio, *Evotime;
    //double **Ek_evo;
    
    C_Storage();
    ~C_Storage();
    
    int AllocateSpace(long Record_times, long Totalbin);
    int InitializeStorage(long Record_times, long Totalbin);
    double **allocate_2d_array(long n, long m);
    
private:
  	 
}; // store the changing things


/** define a Heaviside function **/
double Heaviside(double x);

/** define a function to find the bin for a specific mass **/
long findbin(double m_to_locate);

/***** external variables: *****/
#ifdef MAIN_C

/** main solution array: **/
double *n;
//double Ek_rel[Totalbin];

/** constant array: **/
double *m, *m_inearth;                                          // auxiliary array for m[]
double *v_d_set, *v_d_tur, **v_d_relative;                      // velocity dispersion for each bin
double *r;                                                      // radius for m[]
double *dens;                                                   // density for m[]
double **GCS;                                                   // geometrical crossing section
double **MEV;                                                   // mutual escape velocity in direct impact

/** need to update every timestep!!!!!!!!: **/
C_Parameter *para;
C_FileOp *fop;
C_Delta *delta;
C_Storage *storage;
double pMass_tot, pMass_back, pMass, pMass_ghost, GPratio, pMass_con;
//double totEk;

/** timer: **/
clock_t start_t, end_t;

/** for MPI **/
int myid, master, numprocs;
long myregime;
double *tempdelta_coag, *tempdelta_frag, tempghost, tempback;
MPI_Status status;

#else

/** main solution array: **/
extern double *n;
//extern double Ek_rel[Totalbin];

/** constant array: **/
extern double *m, *m_inearth;                                   // auxiliary array for m[]
extern double *v_d_set, *v_d_tur, **v_d_relative;               // velocity dispersion for each bin
extern double *r;                                               // radius for m[]
extern double *dens;                                            // density for m[]
extern double **GCS;                                            // geometrical crossing section
extern double **MEV;                                            // mutual escape velocity in direct impact


/** need to update every timestep!!!!!!!!: **/
extern C_Parameter *para;
extern C_FileOp *fop;
extern C_Delta *delta;
extern C_Storage *storage;
extern double pMass_tot, pMass_back, pMass, pMass_ghost, GPratio, pMass_con;
//extern double totEk;

/** timer: **/
extern clock_t start_t, end_t;

/** for MPI **/
extern int myid, master, numprocs;
extern long myregime;
extern double *tempdelta_coag, *tempdelta_frag, tempghost, tempback;
extern MPI_Status status;

#endif /* #ifdef MAIN_C */


#endif /* defined(__planetesimal__pre_define__) */
