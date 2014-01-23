//
//  main.cpp
//  planetesimal
//
//  Created by Rixin Li on 8/20/13.
//  Copyright (c) 2013 李 日新. All rights reserved.
//

#ifndef MAIN_C
#define MAIN_C

#include <iostream>
#include "pre_define.h"
#include "condensation.h"
#include "coagulation.h"
#include "fragmentation.h"
#include "disruption.h"
#include "cog_frag_ingrav.h"

// declare some sub functions
int initialize();
int update_mass();
int planetesimal();
int GI();


int main(int argc, const char * argv[])
/* main function */
{
    fop = (C_FileOp *)malloc(sizeof(C_FileOp));
    para = (C_Parameter *)malloc(sizeof(C_Parameter));
    
    /*************************************/
    /***** read the output file path *****/
    /*************************************/
    if (argc < 2) {
        printf("USAGE: %s <data_output_path>\n", argv[0]);
        exit(1);
    } else {
        switch (argc)
        {
            case 2:
                strcpy(fop->outputpath, argv[1]);
                break;
            default:
                printf("Only one arguments: USAGE: %s, <data_output_path>\n", argv[0]);
                exit(1);
        }
    }
    /********** end of block **********/
    
    /***********************************/
    /****** execute the main body ******/
    /***********************************/
    
    printf("Begin the whole program.\n");
    fop->openinput();
    while ((fop->readinput(para)) == 0) {
        printf("\nNow processing task %d.\n\n", para->Order_number);
        fop->filename_gen(para);
        fop->openoutput();
        para->DeriveParameter();
        planetesimal();
        fop->closeoutput();
        printf("\nTask %d done!\n\n", para->Order_number);
    }
    
    printf("The whole program ends.\n");
    
    /********** end of block ***********/
    return 0;
    
}

int planetesimal()
/* main body */
{
    
    /***********************************************/
    /***** define variables and allocate space *****/
    /***********************************************/
    start_t = clock();                          /* program timing */
	long i = 0, j = 0, k = 0, count_column = 0;
	double mytime = 0, mytime_inkyr = 0;        /* set up time */
	
    long count_iter = 0, count_iter2 = 0;
    
    delta = (C_Delta *)malloc(sizeof(C_Delta));
    storage = (C_Storage *)malloc(sizeof(C_Storage));
    delta->AllocateSpace(para->Totalbin);
    delta->InitializeDelta(para->Totalbin);
    storage->AllocateSpace(para->Record_times, para->Totalbin);
    storage->InitializeStorage(para->Record_times, para->Totalbin);
    
    /********** end of block **********/
    
    /********************************************/
    /***** initialize and open output files *****/
    /********************************************/
	printf("Task %d: initialize...\n", para->Order_number);
    initialize();
    
	mytime_inkyr = (double)mytime/Kyr;
    
    for(i = 0; i < para->Totalbin; i++) {
    /* storage the initial values */
        storage->result[count_column][i] = n[i];
    }
    storage->GPratio[0] = GPratio;
    storage->Evotime[0] = 0;
    count_column++;
    
	//printf("Output file path: %s \n",outputpath);

	fop->print_paras();
    fprintf(fop->fm, "%e\t%e\t%e\t%e\t%e\t%e\t%e\t\n", mytime_inkyr, pMass/M_earth, pMass_back/M_earth, pMass_ghost/M_earth,
            pMass_tot/M_earth, para->gMass_tot_i/M_earth, GPratio);
	/********** end of block **********/
    
    /********************************************/
    /***** begin the iteration of main body *****/
    /********************************************/
	printf("Task %d: begin iteration...\n", para->Order_number);
    
    do {
		count_iter2++;  /*count_iter, count_iter2 are two counters in charge of when to dump data in the files*/
        delta->InitializeDelta(para->Totalbin);

        /** condensation from background **/
        condensation();
        
        for (i = 0; i < para->Totalbin; i++) {
            
			for (j = i; j < para->Totalbin; j++) {
                //printf("i = %ld, j = %ld, Totalbin = %d\n", i, j, Totalbin);
                
                
                if (r[i] > para->R_in_grav_regime || r[j] > para->R_in_grav_regime) {
                    
                    // if in gravity regime, calculate the coagulation and fragmentation
                    cog_frag_ingrav(i, j);
                    
                } else {
                    
                    // coagulate and fragment each bin
                    coagulation(i, j);
                    
                    // fragment each bin
                    fragmentation(i, j);
                }
                //*/
                
                /*
                // coagulate and fragment each bin
                coagulation(i, j);
                
                // fragment each bin
                fragmentation(i, j);
                //*/
			}
		}
        
        for (i = 0; i < para->Totalbin; i++) {
            delta->sum_delta_n(i);
            n[i] += delta->delta_n[i];
            if (n[i] < 0) {
                printf("Task %d: n[%ld] became negative, chabie %e, qianjiu_coag %ld, qianjiu_frag %ld\n", para->Order_number, i, pMass_tot-pMass_con, para->countnumber[0], para->countnumber[1]);
                if (pMass_tot - pMass_con > 10 * M_earth) {
                    printf("Total mass deviates from original too much.\n");
                    exit(0);
                }
                //n[i] = 0; continue; //!@#
                printf("%ld\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", i, n[i], delta->delta_n_coag[i], delta->delta_n_frag[i], delta->delta_n_accr[i], delta->delta_n_cond[i], delta->delta_n_disr[i], delta->delta_n_hitr[i]);
                fprintf(fop->fo, "Error: n[%ld] became negative.\n", i);
                end_t = clock();
                printf("Task %d: finalize...\n", para->Order_number);
                printf("Task %d: Time_cost: %fs\n", para->Order_number, (double)(end_t - start_t)/CLOCKS_PER_SEC);
                fprintf(fop->fo, "Real Evolution time: %fKyr\n", mytime/Kyr);
                fprintf(fop->fo, "Time_cost: %fs\n", (double)(end_t - start_t)/CLOCKS_PER_SEC);
                
                delete delta;
                delete storage;
        
                return 0;
            }
        }
        update_mass();
        //printf("%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", mytime_inkyr, pMass/M_earth, pMass_back/M_earth, pMass_ghost/M_earth, pMass_tot/M_earth, para->gMass_tot_i/M_earth, GPratio, pMass_tot-pMass_con);
        /*
        if (pMass_back >= 0.5*para->gMass_tot_i) {
            GI();
            update_mass();
            printf("%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", mytime_inkyr, pMass/M_earth, pMass_back/M_earth, pMass_ghost/M_earth,
                    pMass_tot/M_earth, para->gMass_tot_i/M_earth, GPratio, pMass_tot-pMass_con);

        }
        //*/
        /** the ending of one loop, recording **/
		mytime += para->Timestep;
		count_iter++;
        
        //if (count_iter2 % 10 == 0) {
        //   printf("another 10 yrs\n");
        //}
        
        if(count_iter * para->Timestep_inyr >= para->Record_interval_inyr)
		{
			mytime_inkyr = mytime/Kyr;
            
            //fprintf(fop->fm, "%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", mytime_inkyr, pMass/M_earth, pMass_back/M_earth, pMass_ghost/M_earth,                   pMass_tot/M_earth, para->gMass_tot_i/M_earth, GPratio, pMass_tot-pMass_con);

			for(k = 0; k < para->Totalbin; k++) {
				storage->result[count_column][k] = n[k];
                storage->coag[count_column][k] = delta->delta_n_coag[k];
                storage->frag[count_column][k] = delta->delta_n_frag[k];
                storage->accr[count_column][k] = delta->delta_n_accr[k];
                storage->cond[count_column][k] = delta->delta_n_cond[k];
                storage->disr[count_column][k] = delta->delta_n_disr[k];
                storage->hitr[count_column][k] = delta->delta_n_hitr[k];
			}
            storage->GPratio[count_column] = GPratio;
            storage->Evotime[count_column] = mytime_inkyr;
			count_column++;
			count_iter = 0;
		}
        if (count_iter2 % (long)(10*para->Record_interval_inyr) == 0)
            printf("Task %d: Now evolves to %f yr\n", para->Order_number, mytime/yr);
            fprintf(fop->fm, "%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", mytime_inkyr, pMass/M_earth, pMass_back/M_earth, pMass_ghost/M_earth,
                pMass_tot/M_earth, para->gMass_tot_i/M_earth, GPratio, pMass_tot-pMass_con);
	} while(mytime <= para->Totaltime);
    
	/********** end of block **********/
    
    /*****************************/
    /***** finalize the code *****/
    /*****************************/
    
	end_t = clock();
	printf("Task %d: finalize...\n", para->Order_number);
    printf("Task %d: Time_cost: %fs\n", para->Order_number, (double)(end_t - start_t)/CLOCKS_PER_SEC);
    fprintf(fop->fo, "Real Evolution time: %fKyr\n", mytime/Kyr);
    fprintf(fop->fo, "Time_cost: %fs\n", (double)(end_t - start_t)/CLOCKS_PER_SEC);
    for (int ii = 0; ii < 10; ii++) {
        printf("countnumber[%d] = %ld\n", ii, para->countnumber[ii]);
    }
    fop->dataoutput(count_column);

    delete delta;
    delete storage;
    
    return 0;
    /********** end of block **********/
}


/**************sub functions***************/




int r_gen()
/* compute radius for each mass bin */
{
	double temp;
	r = (double *)malloc(sizeof(double) * para->Totalbin);
	long i;
	for(i = 0; i < para->Totalbin; i++)
	{
		temp = 3 * m[i]/4/Pi/dens[i];
		r[i] = pow(temp, (double)1/3);
	}
	return 0;
}

int dens_gen()
/* density for planetesimals, right now we are assume uniform solid-density profile
 * independent of planetesimal masses. */
{
	dens = (double *)malloc(sizeof(double) * para->Totalbin);
	long i;
	for(i = 0; i < para->Totalbin; i++) {
        //we may find a more accurate planetesimal density profile later on
		dens[i] = para->Rho_solid;
	}
    return 0;
}

int v_d_gen()
/* velocity dispersion for planetesimals, since the region we consider is located around
 * the density and pressure maxima, sandwiched between gas with sub-Keplerian velocity 
 * and gas with super-Keplerian velocity, which means gas drag is not the main source of
 * velocity dispertion. So at this time, we all take account of the relative turbulent velocities
 * and differential settling.
 */
{
    
    /* *****************************************************************************************************************************
     * From A's book:
     * v_settle can be calculated by F_grav = |F_drag|
     *
     * where F_grav = m * Omega_Kepler^2 * z = 4Pi/3 * r^3 * Rho_solid * Omega_Kepler^2 * z
     * here we use the dust scale height as z instead of gas scale height
     * and (H_dust/H_gas) ~ 0.01, which means z = 0.01 * Scale_h
     * 
     * F_drag is in Stokes region
     * so, F_drag = -C/2 * Pi * r^2 * Rho_gas * v^2, where C ~ 24 Re^(-0.6)
     * Re = 2 * r * v / niu_m, and niu_m (i.e. molecular viscosity) ~ Lambda_gas * C_gas
     * Then from CY2010: Re is defined as r * v / (Lambda_gas * C_gas), thus we decide to use this definition.
     *
     * From all above, we can get
     *  ConstA = Lambda_gas * C_gas;
     *  ConstB = Rho_solid/Rho_gas * Omega_Kepler*Omega_Kepler * 0.01*Scale_h / 9;
     * for Re > 1
     *    v_settle = (ConstB * ConstA^(-0.6) * r^(1.6))^(5./7)
     * for Re < 1
     *    v_settle = ConstB/ConstA * r^2
     *
     * We assume that meter-size boulders has the fastest settling velocity: v_settle_max
     * Once v_settle exceed it, since large planetesimals will strongly
     * concentrate on the midplane, which means they won't have too large v_settle
     * we'll make v_settle decrease with the increasing radius by times a damping factor
     *
     * Besides, from B2008 or CY2010, v_turbulence ~ sqrt(alpha) * C_gas
     * *****************************************************************************************************************************
     *
     * On the other hand, according to BD2008, Stokes number -- St = Omega_Kepler * r * Rho_solid / (C_gas * Rho_gas) * Alpha^(2q-1)
     * here, parameter q determines whether turbulent diffusion in the disk is realized by big turbulent eddies moving slow (q = 1)
     * or by small turbulent eddies moving fast (q = 0). Throughout this code we will assume that q = 1/2 following Cuzzi et al.
     * (2001) and Schräpler & Henning (2004) unless otherwise stated.
     * thence, v_turbulence = Alpha^q * C_gas, Alpha ~ 0.001
     *
     * for those particles whose St less than 1, v_settle = z * St * Omega_Kepler
     * but in order to incorporate those whose St larger than 1 and not exceed the vertically projected Kepler velocity
     * z * V_Kepler / r, BD2008 use v_settle = z * Omega_Kepler * (St / 1 + St)
     * here, we will use the scale height of dust Scale_h_d = Scale_h * sqrt(Alpha / min(St, 0.5) / (1+St))
     * so v_settle = Scale_h * Alpha^0.5 * Omega_Kepler * St / (1+St)^1.5 / min(St, 0.5)
     * the advantage of this method is such velocity automatically decreses while r become large
     *
     * don't know why, such method will lead to relatively small velocity dispersion
     *
     * *****************************************************************************************************************************
     *
     * As for the v_turbulence, since large planetesimals also won't significantly affected by turbulence,
     * so we handle it by adding a quickly damping factor pow(r/0.001, -0.2), for particles with radius of ~1mm will couple
     * with the turbulent eddy very well
     *
     * *****************************************************************************************************************************
     *
     *
     */
    
    long i = 0, j = 0, i_max = 0;
    double v_turbulence, v_settle;
    double ConstA, ConstB, r_cutoff = 1.0, v_settle_max, Re_number; // seems r_cutoff = 1.75 will lead to the most plausible v_d distribution
    //double St, Scale_h_d, loser;
    
    v_d_set = (double *)malloc(sizeof(double) * para->Totalbin);
    v_d_tur = (double *)malloc(sizeof(double) * para->Totalbin);
    v_d_relative = (double **)malloc(sizeof(double *) * para->Totalbin);
    for (i = 0; i < para->Totalbin; i++) {
        v_d_relative[i] = (double *)malloc(sizeof(double) * para->Totalbin);
    }
    v_turbulence = sqrt(para->Alpha) * para->C_gas;
#ifdef  Asbook_B2008_settling
    ConstA = para->Lambda_gas * para->C_gas;
    ConstB = para->Rho_solid / para->Rho_gas * para->Omega_Kepler * para->Omega_Kepler * para->mean_H_ratio * para->Scale_h / 9;
#endif
   
    for (i = 0; i < para->Totalbin; i++) {

#ifdef  Asbook_B2008_settling
        v_settle = pow(ConstB  * pow(ConstA, -0.6) * pow(r[i], 1.6), 5./7);
        Re_number = r[i] * v_settle / ConstA;
        
        if (Re_number < 1) {
            v_settle = ConstB /ConstA * r[i] * r[i];
        }

        //printf("r[%ld] = %f, v_set = %f\n", i, r[i], v_settle);
#endif

#ifdef  BD2008_settling
        St = para->Omega_Kepler * r[i] * para->Rho_solid / (para->C_gas * para->Rho_gas);
        loser = St;
        if (loser > 0.5) loser = 0.5;
        Scale_h_d = para->Scale_h * sqrt(para->Alpha / (1+St) / loser );
        v_settle = Scale_h_d * para->Omega_Kepler * St / (1+St);
        //printf("r[%ld] = %f, scale_h_d = %e, v_set = %f\n", i, r[i], Scale_h_d, v_settle);
#endif
        
        v_d_set[i] = v_settle;
        v_d_tur[i] = v_turbulence * pow(r[i]/para->R_well_coupled, para->Index_init-1);

#ifdef  Asbook_B2008_settling
        if (r[i] > r_cutoff) {
            v_settle_max = v_settle;
            break;
        }
    }
    
    i_max = ++i;
    for ( ; i < para->Totalbin; i++) {
        //v_settle = v_settle_max * pow(r[i]/r[i_max], -0.2);
        v_settle = v_d_set[i_max-1] * pow(r[i]/r_cutoff, para->Index_init-1);
        //printf("r[%ld] = %f, m = %e, v_set = %f\n", i, r[i], m[i], v_settle);
        v_d_set[i] = v_settle;
        v_d_tur[i] = v_turbulence * pow(r[i]/para->R_well_coupled, para->Index_init-1);
    
#endif
    }
    for (i = 0; i < para->Totalbin; i++) {
        //printf("%ld %f %f %f\n", i, r[i], v_d_set[i], v_d_tur[i]);
        for (j = 0; j < para->Totalbin; j++) {
    
            if (r[i] > para->R_in_grav_regime) {
                if (r[j] < para->R_in_grav_regime) {
                    // only r[i] in gravity regime
                    v_d_set[i] = 0.0; v_d_tur[i] = 0.0;
                    v_d_relative[i][j] = pow(v_d_set[j]*v_d_set[j]+v_d_tur[j]*v_d_tur[j]+2*Grav*m[i]/r[i], 0.5);
                } else {
                    // both in gravity regime
                    v_d_set[i] = 0.0; v_d_tur[i] = 0.0; v_d_set[j] = 0.0; v_d_tur[j] = 0.0;
                    v_d_relative[i][j] = pow(2*Grav*(m[i]/r[i]+m[j]/r[j]), 0.5);
                }
            } else {
                if (r[j] < para->R_in_grav_regime) {
                    // both not
                    v_d_relative[i][j] = pow((v_d_set[i]-v_d_set[j])*(v_d_set[i]-v_d_set[j])+(v_d_tur[i]*v_d_tur[i])+(v_d_tur[j]*v_d_tur[j]), 0.5);
                } else {
                    // only r[j] in gravity regime
                    v_d_set[j] = 0.0; v_d_tur[j] = 0.0;
                    v_d_relative[i][j] = pow(v_d_set[i]*v_d_set[i]+v_d_tur[i]*v_d_tur[i]+2*Grav*m[j]/r[j], 0.5);
                }
            }
            //printf("v_d_relative[%ld][%ld] = %f\n", i, j, v_d_relative[i][j]);
        }
    }
    return 0;
    
}

int initialize()
/*Initialize the zero time profile*/
{
	long i, j;
    /* a mass-distribution with log-profile will be set up
     * at first, the density is three-times to water
     * M_min is the mass of a 1-decimeter particle
     * M_max is the mass of a 1-km plantesimal
     * but initially, there would be a truncation among mass bins
     * that is, only small planetesimals exist at the beginning
     * m_j = m_k * a^(j-k+1) --> for example: a = (10^9)^(1/250) = 1.096478
     */

	double A, B, temp=0;
    
    /** constants initialize: **/
    n = (double *)malloc(sizeof(double)*para->Totalbin);
    m = (double *)malloc(sizeof(double)*(para->Totalbin+1));
    m_inearth = (double *)malloc(sizeof(double)*para->Totalbin);
    m[0] = para->Volume_factor_solid * pow(para->R_min, 3);
    m[para->Totalbin] = para->Volume_factor_solid * pow(para->R_max, 3);
    /* CY2010:
     * "As the dust-aggregate size increase, the mean
     * collision velocity also increase, leading to stalling of the growth
     * and possibly to fragmentation, once the dust aggregates has reached
     * decimeter size."
     * "At largest sizes, >= 1km, gravity promotes growth"
     */

    GPratio = para->GPratio_i;
    if (GPratio != GPratio_d) {
        temp = para->pMass_i / para->pMass_tot_i;
        //printf("temp (gpmass ratio) = %f\n", temp);
        //temp = 0.5;
        para->pMass_tot_i = gMass_d * GPratio;
        para->pMass_i = para->pMass_tot_i * temp;
        para->pMass_back_i = para->pMass_tot_i - para->pMass_i;
        // if GPratio_i isn't GPratio_d, then we'll assume the total gas mass
        // is gMass_d, and then calculate the corresponding pMass_tot_i
        // pMass_i will be set according to the given ratio of (pMass_i/pMass_tot_i)
    }
    
    pMass_tot = para->pMass_tot_i;
    pMass_back = para->pMass_back_i;
    pMass = para->pMass_i;
    pMass_ghost = 0;
    pArea_tot = 0;
    pArea = (double *)malloc(sizeof(double)*para->Totalbin);
    //printf("pVolume = %e\n", para->pVolume);
    pMass_con = pMass_tot;
    
	m_inearth[0] = (double)m[0]/M_earth;    /* might be useless */
    
    temp += pow(m[0], 1.0-para->Index_init);      //the initial mass distribution is a power-law, this is integrating
	for( i = 1; i < para->Totalbin; i++) {
        /*compute masses for each bin*/
		m[i] = m[i-1] * para->Mass_index;
        m_inearth[i] = (double)m[i]/M_earth;
        
        if(i < para->Truncate_bin){
            temp += pow(m[i],1.0-para->Index_init); /* compute the power-law coeffs. */
        }
    }
    A = para->pMass_i/temp; /* pp_Mass = A * Integrate[m^(-Index_init)]. */
    B = para->pAcc_rate / 2 / temp;
    temp = 0;
    
    for(i = 0; i < para->Truncate_bin; i++) {
        n[i] = A * pow(m[i], -para->Index_init);
        temp += n[i]*m[i];
        delta->delta_n_accr[i] = B * pow(m[i], -para->Index_init);
    }
    //printf("pMass is %f, tempsum is %f.\n", pMass, temp);
    
    for(i = para->Truncate_bin; i < para->Totalbin; i++) {
        n[i] = 0;
        delta->delta_n_accr[i] = 0;
    }
	
    /** initialize the physical parameters: **/
	dens_gen();
    r_gen();
    v_d_gen();
    
    CC = (double **)malloc(sizeof(double *) * para->Totalbin);
    Cij = (double **)malloc(sizeof(double *) * para->Totalbin);
    Fij = (double **)malloc(sizeof(double *) * para->Totalbin);
    Cz = (double **)malloc(sizeof(double *) * para->Totalbin);
    Fz = (double **)malloc(sizeof(double *) * para->Totalbin);
    for (i = 0; i < para->Totalbin; i++) {
        CC[i] = (double *)malloc(sizeof(double) * para->Totalbin);
        Cij[i] = (double *)malloc(sizeof(double) * para->Totalbin);
        Fij[i] = (double *)malloc(sizeof(double) * para->Totalbin);
        Cz[i] = (double *)malloc(sizeof(double) * para->Totalbin);
        Fz[i] = (double *)malloc(sizeof(double) * para->Totalbin);
    }
    for (i = 0; i < para->Totalbin; i++) {
        for (j = 0; j < para->Totalbin; j++) {
            // using gravitional focusing from A's book
            temp = pow(2 * Grav * (m[i] + m[j]) / (r[i] + r[j]) / v_d_relative[i][j], 2);
            CC[i][j] = Pi * (r[i]+r[j]) * (r[i]+r[j]) * v_d_relative[i][j] * (1 + temp);
            
            temp = para->v_frag_th - v_d_relative[i][j];
            temp = (v_d_relative[i][j]/para->v_frag_th)*Heaviside(temp) + Heaviside(0-temp);
            Cij[i][j] = CC[i][j] * (1 - temp);
            Fij[i][j] = CC[i][j] * temp;
            //printf("CC[%ld][%ld] = %e, Cij[%ld][%ld] = %e, Fij[%ld][%ld] = %e\n", i, j, CC[i][j], i, j, Cij[i][j], i, j, Fij[i][j]);
            
        }
    }
    
    for (i = 0; i < para->Totalbin; i++) {
        pArea[i] = Pi * r[i] * r[i];
        pArea_tot += n[i] * pArea[i];
    }
    pArea_ratio = (double *)malloc(sizeof(double)*para->Totalbin);
    partn = (double **)malloc(sizeof(double *) * para->Totalbin);
    for (i = 0; i < para->Totalbin; i++) {
        partn[i] = (double *)malloc(sizeof(double) * para->Totalbin);
    }
    
    for (temp = 0, i = 0; i < para->Totalbin; i++) {
        pArea_ratio[i] = n[i] * pArea[i] / pArea_tot;
        //temp += pArea_ratio[i];
    }
    /*
    if (temp > 0.9999) {
        for (i = 0; i < para->Totalbin; i++) {
            pArea_ratio[i] = pArea_ratio[i] * 0.9999 / temp;
        }
    }
    */
    
    for (i = 0; i < para->Totalbin; i++) {
        for (temp = 0, j = 0; j < para->Totalbin; j++) {
            // here remove 0.5: partn[i][j] = 0.5 * n[i] * pArea_ratio[j];
            partn[i][j] = n[i] * pArea_ratio[j];
            if (partn[i][j] < 1) partn[i][j] = 0;
            temp += partn[i][j];
        }
        if (temp > 1.01*n[i] || temp < 0.99*n[i]) {
            printf("%ld temp deviates from n[%ld]\n", i, i);
        }
        //if (temp * 2 > n[i]) {
        //    n[i] = temp;
        //}
    }
    
    /*
    // security mechanism
    for (i = 0; i < para->Totalbin; i++) {
        for (j = 0; j < para->Totalbin; j++) {
            Cz[i][j] = partn[i][j] * partn[j][i] * Cij[i][j] * (double)(para->Timestep) / para->pVolume;
            Fz[i][j] = partn[i][j] * partn[j][i] * Fij[i][j] * (double)(para->Timestep) / para->pVolume;
        }
    }
    // Cz[i][j] = Cz[j][i]
    double coag_z = 0, frag_z = 0;
    for (i = 0; i < para->Totalbin; i++) {
        for (j = 0; j < para->Totalbin; j++) {
            coag_z += Cz[i][j]; frag_z += Fz[i][j];
        }
        temp = coag_z + frag_z;
        if (temp > n[i]) {
            for (j = 0; j < para->Totalbin; j++) {
                Cz[i][j] = Cz[i][j] * n[i] / temp;
                Fz[i][j] = Fz[i][j] * n[i] / temp;
            }
            para->countnumber[6]++;
        }
    }
    //*/
    
    //printf("accretion: %f\n", para->pAcc_rate/M_earth);
	return 0;
}


int GI()
/* GI happen */
{
    long i, j, k;
    double temp, A;
    temp = para->Volume_factor_solid * 500 * 500 * 500;
    j = findbin(temp);
    temp = para->Volume_factor_solid * 1000 * 1000 * 1000;
    k = findbin(temp);
    temp = 0;
    for(i = j; i < k; i++) {
        temp += pow(m[i],1.0-para->Index_init); /* compute the power-law coeffs. */
    }
    A = pMass_back/2.0/temp;
    pMass_back /= 2.0;
    temp = 0;
    for(i = j; i < k; i++) {
        n[i] += A * pow(m[i], -para->Index_init);
        temp+=n[i]*m[i];
        //printf("n[%ld] = %e\n", i, n[i]);
    }
    printf("temp/pMass_back = %f.\n", temp/pMass_back);
    return 0;
}


int update_mass()
{
	long i, j;
    //double temp = 0;
    //double checkn = 0;
    pMass = 0;
    pArea_tot = 0;
	for(i = 0; i < para->Totalbin; i++) {
		pMass += n[i] * m[i];
        pArea_tot += n[i] * pArea[i];
    }
    
    for (i = 0/*, checkn = 0*/; i < para->Totalbin; i++) {
        pArea_ratio[i] = n[i] * pArea[i] / pArea_tot;
        //checkn += pArea_ratio[i];
    }
    
    /*
    if (checkn > 0.9999) {
        for (i = 0; i < para->Totalbin; i++) {
            pArea_ratio[i] = pArea_ratio[i] * 0.9999 / checkn;
        }
    }
     */
    for (i = 0; i < para->Totalbin; i++) {
        //checkn = 0;
        for (j = 0; j < para->Totalbin; j++) {
            //here remove 0.5 : partn[i][j] = 0.5 * n[i] * pArea_ratio[j];
            partn[i][j] = n[i] * pArea_ratio[j];
            if (partn[i][j] < 1) partn[i][j] = 0;
            //checkn += partn[i][j];
            //else if (i > 500 || j > 500)
                //printf("partn[%ld][%ld] = %e\tn[%ld] = %e\tn[%ld] = %e\n", i, j, partn[i][j], i, n[i], j, n[j]);
        }
        /*
        if (checkn * 2 > n[i]) {
            //printf("error here, checkn = %e, n[%ld] = %e, diff = %e.\n", checkn, i, n[i], checkn-n[i]);
            n[i] = checkn;
        }
        */
    }
    
    /*
    // security mechanism
    for (i = 0; i < para->Totalbin; i++) {
        for (j = 0; j < para->Totalbin; j++) {
            Cz[i][j] = partn[i][j] * partn[j][i] * Cij[i][j] * (double)(para->Timestep) / para->pVolume;
            Fz[i][j] = partn[i][j] * partn[j][i] * Fij[i][j] * (double)(para->Timestep) / para->pVolume;
        }
    }
    // Cz[i][j] = Cz[j][i]
    double coag_z = 0, frag_z = 0;
    for (i = 0; i < para->Totalbin; i++) {
        for (j = 0; j < para->Totalbin; j++) {
            coag_z += Cz[i][j]; frag_z += Fz[i][j];
        }
        temp = coag_z + frag_z;
        if (temp > n[i]) {
            for (j = 0; j < para->Totalbin; j++) {
                Cz[i][j] = Cz[i][j] * n[i] / temp;
                Fz[i][j] = Fz[i][j] * n[i] / temp;
            }
            para->countnumber[6]++;
        }
    }
    //*/
    
    // consider the accretion from outer disk to background
    pMass_back += para->pAcc_rate / 2;
    // the calculation of delta_n_accr[i] is put into the initialize()
    
	pMass_tot = pMass + pMass_back + pMass_ghost;
	GPratio = pMass_tot / para->gMass_tot_i;
    
    pMass_con += para->pAcc_rate;
    
	return 0;
}

#endif


