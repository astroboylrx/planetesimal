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
#include "cog_frag_ingrav.h"

// declare some sub functions
int initialize();
int update_mass();
int planetesimal();

int main(int argc, char * argv[])
/* main function */
{
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    
    master = 0;
    int breakornot = 0;
    
    fop = (C_FileOp *)malloc(sizeof(C_FileOp));
    para = (C_Parameter *)malloc(sizeof(C_Parameter));
    
    /*************************************/
    /***** read the output file path *****/
    /*************************************/
    if (myid == master) {
        if (argc < 2) {
            printf("USAGE: %s <data_output_path> [-c]\n", argv[0]);
            exit(1);
        } else {
            switch (argc)
            {
                case 2:
                    strcpy(fop->outputpath, argv[1]);
                    break;
                case 3:
                    if (strcmp(argv[1], "-c") == 0 || strcmp(argv[1], "c") == 0) {
                        fop->conti = 1;
                        strcpy(fop->outputpath, argv[2]);
                        break;
                    }
                    strcpy(fop->outputpath, argv[1]);
                    strcpy(fop->cont, argv[2]);
                    if (strcmp(fop->cont, "-c") == 0 || strcmp(fop->cont, "c") == 0) {
                        fop->conti = 1;
                    } else {
                        fop->conti = 0;
                    }
                    //printf("fop->cont = %s, fop->conti = %d, argv[2] = %s, argv[1] = %s\n", fop->cont, fop->conti, argv[2], argv[1]);
                    break;
                    /********** waiting for implementation **********/
                    
                default:
                    printf("Only one arguments: USAGE: %s, <data_output_path> [-c]\n", argv[0]);
                    exit(1);
            }
        }
    }
        /********** end of block **********/
    
    /***********************************/
    /****** execute the main body ******/
    /***********************************/
    if (myid == master) {
        printf("Begin the whole program.\n");
        fop->openinput();
    }
    
    while (1) {
        MPI_Barrier(MPI_COMM_WORLD);
        // breakornot is necessary, you need to tell there is no more work
        if (myid == master) {
            if (fop->readinput(para) != 0) {
                breakornot = 1;
                MPI_Bcast(&breakornot, 1, MPI_INT, master, MPI_COMM_WORLD);
                breakornot = 0;
                break;
            } else {
                breakornot = 0;
                MPI_Bcast(&breakornot, 1, MPI_INT, master, MPI_COMM_WORLD);
            }
        } else {
            MPI_Bcast(&breakornot, 1, MPI_INT, master, MPI_COMM_WORLD);
            if (breakornot == 1) {
                breakornot = 0;
                break;
            }
        }
        
        MPI_Bcast(&para->Order_number, 1, MPI_INT, master, MPI_COMM_WORLD);
        MPI_Bcast(&para->Totalbin, 1, MPI_LONG, master, MPI_COMM_WORLD);
        MPI_Bcast(&para->R_min, 1, MPI_DOUBLE, master, MPI_COMM_WORLD);
        MPI_Bcast(&para->R_max, 1, MPI_DOUBLE, master, MPI_COMM_WORLD);
        MPI_Bcast(&para->Totaltime_inKyr, 1, MPI_DOUBLE, master, MPI_COMM_WORLD);
        MPI_Bcast(&para->Timestep_inyr, 1, MPI_DOUBLE, master, MPI_COMM_WORLD);
        MPI_Bcast(&para->Record_interval_inyr, 1, MPI_DOUBLE, master, MPI_COMM_WORLD);
        MPI_Bcast(&para->v_frag_th, 1, MPI_DOUBLE, master, MPI_COMM_WORLD);
        MPI_Bcast(&para->v_disr_th, 1, MPI_DOUBLE, master, MPI_COMM_WORLD);
        MPI_Bcast(&para->pMass_tot_i, 1, MPI_DOUBLE, master, MPI_COMM_WORLD);
        MPI_Bcast(&para->pMass_i, 1, MPI_DOUBLE, master, MPI_COMM_WORLD);
        MPI_Bcast(&para->Rho_solid, 1, MPI_DOUBLE, master, MPI_COMM_WORLD);
        MPI_Bcast(&para->Acc_rate, 1, MPI_DOUBLE, master, MPI_COMM_WORLD);
        MPI_Bcast(&para->GPratio_i, 1, MPI_DOUBLE, master, MPI_COMM_WORLD);
        MPI_Bcast(&para->R_inau, 1, MPI_DOUBLE, master, MPI_COMM_WORLD);
        MPI_Bcast(&para->R_trun, 1, MPI_DOUBLE, master, MPI_COMM_WORLD);
        MPI_Bcast(&para->Index_init, 1, MPI_DOUBLE, master, MPI_COMM_WORLD);
        MPI_Bcast(&para->Q_lr, 1, MPI_DOUBLE, master, MPI_COMM_WORLD);
        MPI_Bcast(&para->Alpha, 1, MPI_DOUBLE, master, MPI_COMM_WORLD);
        MPI_Bcast(&para->Fa, 1, MPI_DOUBLE, master, MPI_COMM_WORLD);
        MPI_Bcast(&para->Z_rel, 1, MPI_DOUBLE, master, MPI_COMM_WORLD);
        
        if (myid == master) {
            printf("\nNow processing task %d.\n\n", para->Order_number);
            fop->filename_gen(para);
            fop->openoutput();
        }
        MPI_Barrier(MPI_COMM_WORLD);
        para->DeriveParameter();
        planetesimal();
        
        if (myid == master) {
            fop->closeoutput();

        }
        printf("Processor %d: Task %d done!\n", myid, para->Order_number);

    }
    if (myid == master) {
        printf("The whole program ends.\n");
    }
    printf("Processor %d: Ready to Finalize.\n", myid);
    MPI_Finalize();
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
	double *errortest, errormass = 0;
    long count_iter = 0, count_iter2 = 0;
    
    errortest = (double *)malloc(sizeof(double)*para->Totalbin);
    delta = (C_Delta *)malloc(sizeof(C_Delta));
    storage = (C_Storage *)malloc(sizeof(C_Storage));
    delta->AllocateSpace(para->Totalbin);
    delta->InitializeDelta(para->Totalbin);
    storage->AllocateSpace(para->Record_times, para->Totalbin);
    storage->InitializeStorage(para->Record_times, para->Totalbin);
    tempdelta_coag = (double *)malloc(sizeof(double)*para->Totalbin);
    tempdelta_frag = (double *)malloc(sizeof(double)*para->Totalbin);
    
    /********** end of block **********/
    
    /********************************************/
    /***** initialize and open output files *****/
    /********************************************/
	printf("Task %d: Processor %d initialize...\n", para->Order_number, myid);

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
    if (myid == master) {
        fop->print_paras();
        fprintf(fop->fm, "%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", mytime_inkyr, pMass/M_earth, pMass_back/M_earth, tempghost/M_earth,
            pMass_tot/M_earth, para->gMass_tot_i/M_earth, GPratio, pMass_tot-pMass_con, errormass);
    }
	
	/********** end of block **********/
    
    /********************************************/
    /***** begin the iteration of main body *****/
    /********************************************/
	if (myid == master) {
        printf("Task %d: begin iteration...\n", para->Order_number);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    tempback = para->pMass_back_i;
    tempghost = 0;
    printf("Processor %d: my regime is %ld to %ld\n", myid, long(round(para->Totalbin * (1 - sqrt(1-float(myid)/float(numprocs))))), long(round(para->Totalbin * (1 - sqrt(1-float(myid+1)/float(numprocs)))))-1);
    do {
		count_iter2++;  /*count_iter, count_iter2 are two counters in charge of when to dump data in the files*/
        delta->InitializeDelta(para->Totalbin);

        /** condensation from background **/
        if (myid == master) {
            condensation();
        }
                
        //myregime = para->Totalbin / numprocs * myid;
        myregime = long(round(para->Totalbin * (1 - sqrt(1-float(myid)/float(numprocs)))));
        // don't need to -1 because there is <
        for (i = myregime; i < long(round(para->Totalbin * (1 - sqrt(1-float(myid+1)/float(numprocs))))); i++) {
            
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
                //printf("Processor %d: i = %ld, j = %ld, back = %e, ghost = %e\n", myid, i, j, pMass_back, pMass_ghost);
			}
 
            
		}
        
        MPI_Barrier(MPI_COMM_WORLD);
        //printf("Processor %d: i = %ld, j = %ld, back = %e, ghost = %e\n", myid, i, j, pMass_back, pMass_ghost);
        MPI_Allreduce(delta->delta_n_coag, tempdelta_coag, int(para->Totalbin), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(delta->delta_n_frag, tempdelta_frag, int(para->Totalbin), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        if (myid == master) {
            for (i = 0; i < para->Totalbin; i++) {
                delta->delta_n[i] = tempdelta_coag[i] + tempdelta_frag[i] + delta->delta_n_accr[i] + delta->delta_n_cond[i];//
                n[i] += delta->delta_n[i];
                
                if (isnan(n[i])) {
                    printf("n[%ld] = %e, c = %e, f = %e, a = %e, co = %e, n = %e \n", i, n[i], tempdelta_frag[i], tempdelta_frag[i], delta->delta_n_accr[i], delta->delta_n_cond[i], tempdelta_coag[i] + tempdelta_frag[i] + delta->delta_n_accr[i] + delta->delta_n_cond[i]);
                    printf("1\n");
                }
                
                //printf("n[%ld] = %e, coag[%ld] = %e, frag[%ld] = %e, accr[%ld] = %e, cond[%ld] = %e\n", i, n[i], i, tempdelta_coag[i], i, tempdelta_frag[i], i, delta->delta_n_accr[i], i, delta->delta_n_cond[i]);
                
                if (n[i] < 0) {
                    //
                    if (n[i] < -0.1) {
                        fprintf(fop->fo, "Task %d: n[%ld] = %e, M_error = %e, c[%ld] = %e, f[%ld] = %e\n", para->Order_number, i, n[i], pMass_tot-pMass_con, i,tempdelta_coag[i], i, tempdelta_frag[i]);
                    }
                    errormass += fabs(n[i] * m[i]);
                    n[i] = 0;
                    para->countnumber[8]++;
                    continue; //!@#
                    //printf("%ld\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", i, n[i], delta->delta_n_coag[i], delta->delta_n_frag[i], delta->delta_n_accr[i], delta->delta_n_cond[i], delta->delta_n_disr[i], delta->delta_n_hitr[i]);
                    printf("%ld\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", i, n[i], delta->delta_n_coag[i], delta->delta_n_frag[i], delta->delta_n_accr[i], delta->delta_n_cond[i], pMass_tot-pMass_con, errormass);
                    fprintf(fop->fo, "Error: n[%ld] became negative.\n", i);
                    end_t = clock();
                    printf("Task %d: finalize...\n", para->Order_number);
                    printf("Task %d: Time_cost: %fs\n", para->Order_number, (double)(end_t - start_t)/CLOCKS_PER_SEC);
                    fprintf(fop->fo, "Real Evolution time: %fKyr\n", mytime/Kyr);
                    fprintf(fop->fo, "Time_cost: %fs\n", (double)(end_t - start_t)/CLOCKS_PER_SEC);
                    delete delta;
                    delete storage;
                    exit(1);
                    //return 0;
                }
            }
            
        }
        
        MPI_Bcast(n, para->Totalbin, MPI_DOUBLE, master, MPI_COMM_WORLD);
        //MPI_Barrier(MPI_COMM_WORLD);
        update_mass();
        if (myid == master) {
            if (fabs(pMass_tot - pMass_con) > 10 * M_earth) {
                printf("Total mass deviates from original too much.\n");
            //exit(5);
            }
        }
        
        /*
        if (myid == master) {
            if (fabs((pMass_tot-pMass_con)/lasterror) > 5 && count_iter2 > 1) {
                double maxdelta = 0;
                long maxindex = 0;
                printf("Got you. lasterror = %e, pMass_tot-pMass_con = %e\n", lasterror, pMass_tot-pMass_con);
                for (i = 0; i < para->Totalbin; i++) {
                    //printf("n[%ld] = %e, delta->delta->n[%ld] = %e\n", i, n[i], i, delta->delta_n[i]);
                    if (maxdelta < fabs(delta->delta_n[i])) {
                        maxdelta = fabs(delta->delta_n[i]);
                        maxindex = i;
                    }
                }
                printf("maxdelta is delta_n[%ld] = %e\n", maxindex, maxdelta);
                printf("1\n");
                
            }
            lasterror = fabs(pMass_tot - pMass_con);
        }
         */
        //printf("%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", mytime_inkyr, pMass/M_earth, pMass_back/M_earth, pMass_ghost/M_earth, pMass_tot/M_earth, para->gMass_tot_i/M_earth, GPratio, pMass_tot-pMass_con);

        /** the ending of one loop, recording **/
		mytime += para->Timestep;
        mytime_inkyr = mytime/Kyr;
		count_iter++;
        
        if(count_iter * para->Timestep_inyr >= para->Record_interval_inyr)	{
            //fprintf(fop->fm, "%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", mytime_inkyr, pMass/M_earth, pMass_back/M_earth, pMass_ghost/M_earth,                   pMass_tot/M_earth, para->gMass_tot_i/M_earth, GPratio, pMass_tot-pMass_con);
			for(k = 0; k < para->Totalbin; k++) {
				storage->result[count_column][k] = n[k];
                //storage->coag[count_column][k] = delta->delta_n_coag[k];
                storage->coag[count_column][k] = tempdelta_coag[k];
                //storage->frag[count_column][k] = delta->delta_n_frag[k];
                storage->frag[count_column][k] = tempdelta_frag[k];
                storage->accr[count_column][k] = delta->delta_n_accr[k];
                storage->cond[count_column][k] = delta->delta_n_cond[k];
                /*
                 storage->disr[count_column][k] = delta->delta_n_disr[k];
                 storage->hitr[count_column][k] = delta->delta_n_hitr[k];
                 */
			}
            storage->GPratio[count_column] = GPratio;
            storage->Evotime[count_column] = mytime_inkyr;
			count_column++;
			count_iter = 0;
		}
        
        if (myid == master) {
            if ((long)(10000000*count_iter2*para->Timestep_inyr) % (long)(10000000*para->Record_interval_inyr) == 0) {
                printf("Task %d: Evolves to %7.2f yr, M_error = %6.2e (%6.2e), -: %ld, tot =  %6.2e (%6.2e), pM = %6.2e (%6.2e, %6.2e)\n", para->Order_number, mytime/yr, pMass_tot-pMass_con, errormass, para->countnumber[8], pMass_tot, pMass_con, pMass, pMass_back, tempghost);
            }
            fprintf(fop->fm, "%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", mytime_inkyr, pMass/M_earth, pMass_back/M_earth, tempghost/M_earth,
                    pMass_tot/M_earth, para->gMass_tot_i/M_earth, GPratio, pMass_tot - pMass_con, errormass);
            
        }
    } while(mytime <= para->Totaltime);
    	
	/********** end of block **********/
    
    /*****************************/
    /***** finalize the code *****/
    /*****************************/
    
	end_t = clock();
    for (int ii = 0 ; ii < 99; ii++) {
        MPI_Allreduce(&para->countnumber[ii], &para->countnumber[99], 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
        para->countnumber[ii] = para->countnumber[99];
    }
    if (myid == master) {
        printf("Task %d: finalize...\n", para->Order_number);
        printf("Task %d: Time_cost: %fs\n", para->Order_number, (double)(end_t - start_t)/CLOCKS_PER_SEC);
        fprintf(fop->fo, "Real Evolution time: %fKyr\n", mytime_inkyr);
        fprintf(fop->fo, "Time_cost: %fs\n", (double)(end_t - start_t)/CLOCKS_PER_SEC);
        for (int ii = 0; ii < 10; ii++) {
            printf("countnumber[%d] = %ld\n", ii, para->countnumber[ii]);
        }
        fop->dataoutput(count_column);
    }
	
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
    double v_turbulence, v_settle, v_ep_begin_2;
    double ConstA, ConstB, r_cutoff = 1.0, v_settle_max, Re_number; // seems r_cutoff = 1.75 will lead to the most plausible v_d distribution
    double St, Scale_h_d, loser;
    
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
    
    v_ep_begin_2 = 2.0*Grav*(4.0*Pi/3.0*para->Rho_solid*pow(para->R_in_ep_regime, 3.0))/para->R_in_ep_regime;
    /*
    if (myid == 0) {
        printf("v_ep_begin_2 = %f\n", v_ep_begin_2);
        for (i = 0; i < para->Totalbin; i++) {
            if (r[i] > para->R_in_ep_regime) {
                printf("v[%ld] = %f\n", i, sqrt(v_ep_begin_2*pow(para->R_in_ep_regime, 3.0)*pow(r[i], -3.0)));
            } else {
                printf("v[%ld] = %f\n", i, sqrt(2*Grav*m[i]/r[i]));
            }
                
        }
    }
     */
//#define ISO_FEEDING
#ifdef ISO_FEEDING
    for (i = 0; i < para->Totalbin; i++) {
        v_d_tur[i] = para->Omega_Kepler * 10.0 * pow(m[i]/3.0/M_solar, 1/3.0) * para->R_preset;
        //printf("v_d_tur[%ld] is %e\n", i, v_d_tur[i]);
        v_d_set[i] = 0;
    }
    for (i = 0; i < para->Totalbin; i++) {
        for (j = 0; j < para->Totalbin; j++) {
            v_d_relative[i][j] = sqrt(v_d_tur[i]*v_d_tur[i]+v_d_tur[j]*v_d_tur[j]);
            //printf("i = %ld, j = %ld, v_d_relative[i][j] = %e\n", i, j, v_d_relative[i][j]);
        }
    }
#else
    for (i = 0; i < para->Totalbin; i++) {
        //printf("%ld %f %f %f %f %f\n", i, r[i], n[i], v_d_set[i], v_d_tur[i], pow(2*Grav*m[i]/r[i], 0.5));
        for (j = 0; j < para->Totalbin; j++) {
            if (r[i] > para->R_in_ep_regime) {
                v_d_tur[i] = v_ep_begin_2*pow(para->R_in_ep_regime, 3.0)*pow(r[i], -3.0); v_d_set[i] = 0.0;
                if (r[j] > para->R_in_ep_regime) {
                    v_d_tur[j] = v_ep_begin_2*pow(para->R_in_ep_regime, 3.0)*pow(r[j], -3.0); v_d_set[j] = 0.0;
                    v_d_relative[i][j] = pow(v_ep_begin_2*pow(para->R_in_ep_regime, 3.0)*(pow(r[i], -3.0)+pow(r[j], -3.0)), 0.5);
                } else {
                    if (r[j] > para->R_in_grav_regime) {
                        v_d_tur[j] = sqrt(2*Grav*m[j]/r[j]);  v_d_set[j] = 0.0;
                        v_d_relative[i][j] = pow(v_ep_begin_2*pow(para->R_in_ep_regime, 3.0)*pow(r[i], -3.0)+2*Grav*m[j]/r[j], 0.5);
                    } else {
                        v_d_relative[i][j] = pow(v_ep_begin_2*pow(para->R_in_ep_regime, 3.0)*pow(r[i], -3.0)+v_d_set[j]*v_d_set[j]+v_d_tur[j]*v_d_tur[j], 0.5);
                    }
                }
            } else if (r[i] > para->R_in_grav_regime) {
                v_d_tur[i] = sqrt(2*Grav*m[i]/r[i]);  v_d_set[i] = 0.0;
                if (r[j] > para->R_in_ep_regime) {
                    v_d_tur[j] = v_ep_begin_2*pow(para->R_in_ep_regime, 3.0)*pow(r[j], -3.0); v_d_set[j] = 0.0;
                    v_d_relative[i][j] = pow(2*Grav*m[i]/r[i]+v_ep_begin_2*pow(para->R_in_ep_regime, 3.0)*pow(r[j], -3.0), 0.5);
                } else {
                    if (r[j] > para->R_in_grav_regime) {
                        v_d_tur[j] = sqrt(2*Grav*m[j]/r[j]);  v_d_set[j] = 0.0;
                        v_d_relative[i][j] = pow(2*Grav*(m[i]/r[i]+m[j]/r[j]), 0.5);
                    } else {
                        v_d_relative[i][j] = pow(2*Grav*m[i]/r[i]+v_d_set[j]*v_d_set[j]+v_d_tur[j]*v_d_tur[j], 0.5);
                    }
                }
            } else {
                if (r[j] > para->R_in_ep_regime) {
                    v_d_tur[j] = v_ep_begin_2*pow(para->R_in_ep_regime, 3.0)*pow(r[j], -3.0); v_d_set[j] = 0.0;
                    v_d_relative[i][j] = pow(v_d_set[i]*v_d_set[i]+v_d_tur[i]*v_d_tur[i]+v_ep_begin_2*pow(para->R_in_ep_regime, 3.0)*pow(r[j], -3.0), 0.5);
                } else {
                    if (r[j] > para->R_in_grav_regime) {
                        v_d_tur[j] = sqrt(2*Grav*m[j]/r[j]);  v_d_set[j] = 0.0;
                        v_d_relative[i][j] = pow(v_d_set[i]*v_d_set[i]+v_d_tur[i]*v_d_tur[i]+2*Grav*m[j]/r[j], 0.5);
                    } else {
                        v_d_relative[i][j] = pow((v_d_set[i]-v_d_set[j])*(v_d_set[i]-v_d_set[j])+(v_d_tur[i]*v_d_tur[i])+(v_d_tur[j]*v_d_tur[j]), 0.5);
                    }
                }
            }
            
            //printf("i = %ld, j = %ld, v_d_relative[i][j] = %e\n", i, j, v_d_relative[i][j]);
        }
        
    }
#endif
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
    pMass_con = pMass_tot;
    
	m_inearth[0] = (double)m[0]/M_earth;    /* might be useless */
    
    temp += pow(m[0], 1.0-para->Index_init);      //the initial mass distribution is a power-law, this is integrating
	if (myid == master) {
        for(i = 1; i < para->Totalbin; i++) {
            /*compute masses for each bin*/
            
            m[i] = m[0] * pow(para->Mass_index, i) * (1.0/para->Mass_index +(((para->Mass_index-1.0/para->Mass_index)/10000.0)*(rand()%5000+2500.0)));
            if (m[i]<=m[i-1]) printf("mass is wrong at mass[%ld].\n", i);
            m_inearth[i] = (double)m[i]/M_earth;
            //printf("m[%ld] = %f\n", i, m[i]);
            if(i < para->Truncate_bin){
                temp += pow(m[i],1.0-para->Index_init); /* compute the power-law coeffs. */
            }
        }
        
    }
    
    MPI_Bcast(m, para->Totalbin+1, MPI_DOUBLE, master, MPI_COMM_WORLD);
    MPI_Bcast(m_inearth, para->Totalbin, MPI_DOUBLE, master, MPI_COMM_WORLD);
    MPI_Bcast(&temp, 1, MPI_DOUBLE, master, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    
    dens_gen();
    r_gen();
    A = para->pMass_i / temp; /* pp_Mass = A * Integrate[m^(-Index_init)]. */
    B = para->pAcc_rate / 2.0 / temp;
    temp = 0.0;
    for(i = 0; i < para->Truncate_bin; i++) {
        
        n[i] = A * pow(m[i], -para->Index_init);
        //temp += n[i]*m[i];
        if (r[i] < para->R_in_grav_regime) {
            delta->delta_n_accr[i] = B * pow(m[i], -para->Index_init);
        } else {
            delta->delta_n_accr[i] = 0;
        }
        temp += delta->delta_n_accr[i] * m[i];
        if (myid == master) {
            printf("Processor 0: Initially r[%ld] = %e, m[%ld] = %e, n[%ld] = %e\n",i, r[i], i, m[i],  i, n[i]);
        }
    }
    para->pAcc_intoback = fabs(para->Acc_rate - temp);
    
    if (myid == master) {
        printf("Processor 0: Supply rate = %e, pMass_back_i = %e\n", m[0]*1.85e+07/(para->gMass_tot_i*0.8/para->pAcc_rate), para->pMass_back_i);
        printf("Processor 0: pAcc_rate = %e, para->intoback = %e, temp = %e\n", para->pAcc_rate, para->pAcc_intoback, temp);
    }
    //printf("pMass is %f, tempsum is %f.\n", pMass, temp);
    
    for(i = para->Truncate_bin; i < para->Totalbin; i++) {
        n[i] = 0;
        delta->delta_n_accr[i] = 0;
    }
	
    /** initialize the physical parameters: **/
    v_d_gen();
    
    GCS = (double **)malloc(sizeof(double *) * para->Totalbin);
    MEV = (double **)malloc(sizeof(double *) * para->Totalbin);
    for (i = 0; i < para->Totalbin; i++) {
        GCS[i] = (double *)malloc(sizeof(double) * para->Totalbin);
        MEV[i] = (double *)malloc(sizeof(double) * para->Totalbin);
    }
    for (i = 0; i < para->Totalbin; i++) {
        for (j = 0; j < para->Totalbin; j++) {
            // using gravitational focusing from A's book
            //temp = pow(2 * Grav * (m[i] + m[j]) / (r[i] + r[j]) / v_d_relative[i][j], 2);
            /*
            if (myid == 0) {
                printf("r[%ld] = %f, r[%ld] = %f, Grav_focusing_factor = %f\nanother case: %f,  new case: %lf\n", i, r[i], j, r[j], temp,
                       2 * Grav * (m[i] + m[j]) / (r[i] + r[j]) / pow(v_d_relative[i][j], 2),
                       2 * Grav * m[i] * m[j] / (r[i] + r[j]) / (m[i] * v_d_tur[i] * v_d_tur[i] + m[j] * v_d_tur[j] * v_d_tur[j]));
            }
            
            temp = 2 * Grav * (m[i] + m[j]) / (r[i] + r[j]) / pow(v_d_relative[i][j], 2);
            CC[i][j] = Pi * (r[i]+r[j]) * (r[i]+r[j]) * v_d_relative[i][j] * (1 + temp);
            
            temp = para->v_frag_th - v_d_relative[i][j];
            temp = (v_d_relative[i][j]/para->v_frag_th)*Heaviside(temp) + Heaviside(0-temp);
            Cij[i][j] = CC[i][j] * (1 - temp);
            Fij[i][j] = CC[i][j] * temp;
            //printf("CC[%ld][%ld] = %e, Cij[%ld][%ld] = %e, Fij[%ld][%ld] = %e\n", i, j, CC[i][j], i, j, Cij[i][j], i, j, Fij[i][j]);
            */
            
            GCS[i][j] = Pi * (r[i]+r[j]) * (r[i]+r[j]);
            MEV[i][j] = 2 * Grav * (m[i] + m[j]) / (r[i] + r[j]);
            if (isnan(GCS[i][j]) || isnan(MEV[i][j])) {
                printf("1");
            }
        }
    }
    
	return 0;
}


int update_mass()
{
	long i;
    pMass = 0;
    double temptempback;
    //pArea_tot = 0;
    
	for(i = 0; i < para->Totalbin; i++) {
		pMass += n[i] * m[i];
    }

    // consider the accretion from outer disk to background
    // the calculation of delta_n_accr[i] is put into the initialize()
    /*
    if (myid == master) {
        printf("Processor %d:(1) pMass = %e, back = %e, ghost = %e, tot = %e, ttb = %e, tb = %e\n", myid, pMass, pMass_back, tempghost, pMass_tot, temptempback, tempback);
    }
     */
    temptempback = tempback;
    MPI_Allreduce(&pMass_ghost, &tempghost, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&pMass_back, &tempback, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    tempback = tempback - temptempback * (numprocs - 1);

    pMass_back = tempback;
	pMass_tot = pMass + tempback + tempghost;
    /*
    if (myid == master) {
        printf("Processor %d:(2) pMass = %e, back = %e, ghost = %e, tot = %e, ttb = %e, tb = %e\n", myid, pMass, pMass_back, tempghost, pMass_tot, temptempback, tempback);
    }
     */
	GPratio = pMass_tot / para->gMass_tot_i;
    
    pMass_con += para->pAcc_rate;
    
	return 0;
}

#endif


