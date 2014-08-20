//
//  cog_frag_ingrav.cpp
//  planetesimal
//
//  Created by Rixin Li on 11/19/13.
//  Copyright (c) 2013 李 日新. All rights reserved.
//

#include "cog_frag_ingrav.h"
#define Rayleigh_Velocity
#define Test_Safronov_Number


int cog_frag_ingrav(long &i, long &j)
/** calculate the both coagulation and fragmentation in gravity regime between two mass bins **/
{
    if (n[i] == 0 || n[j] == 0) {
        return 0;
    }
    // The basic idea of this subroutine is from LS2012
    //
    // we simply assume that for the whole collisional cross section of target, the projectile will hit any point
    // with the same probability, so we divide B into several segments
    // l + B = R + r

    int count;
    long lr_bin, second_bin, third_bin, k;
    double *B, temp, b, l, b_crit, r_tot, m_tot, gamma, n_alpha, m_alpha, miu, V_esc, V_hr;
    double prob, xi, c1=2.43, c2=-0.0408, c3=1.86, c4=1.08;
    double Q_R, Q_RD, rho1 = 1.0e3, M_lr, z, tempz, check, M_second, M_third;
    
    /**************************************************/
    // here, because we assume the velocity dispersions of planetesimals in each mass bin
    // are subject to a Gaussian distribution, then the relative velocity, which is
    // sqrt(va^2+vb^2) must be subject to a Rayleigh distribution. In order to take this into
    // consideration, we produce a series of random numbers (random_factor) that are subject to a Rayleigh
    // distribution (sigma = 1), then impact_velocity = random_factor * v_d_relative[i][j]
    // the method to produce random numbers as a Rayleigh distribution is sqrt(-2*sigma^2*ln(1-x))
    double impact_velocity;
    //srand(time(NULL));
#ifdef Rayleigh_Velocity
    impact_velocity = v_d_relative[i][j] * sqrt(-2*log(1-(rand()%100000000+0.5/100000000)/100000000.0));
#else
    impact_velocity = v_d_relative[i][j];
#endif

    //printf("impact_velocity / v_d_relative[i][j] = %f\n", sqrt(-2*log(1-(rand()%10000)/10000.0)));
    /**************************************************/
    //while (impact_velocity/v_d_relative[i][j] > 3.71692 || impact_velocity/v_d_relative[i][j] < 0.0447325 ) {
    while (impact_velocity/v_d_relative[i][j] > 1.66511 || impact_velocity/v_d_relative[i][j] < 0.758528 ) {
        //printf("Processor %d: calculate impact velocity again\n", myid);
        para->countnumber[9]++;
        impact_velocity = v_d_relative[i][j] * sqrt(-2*log(1-(rand()%100000000+0.5/100000000)/100000000.0));
    }
    //printf("Processor %d: impact_velocity/v_d_relative[i][j] = %e\n", myid, impact_velocity/v_d_relative[i][j]);
    B = (double *)malloc(sizeof(double)*para->B_totseg);
    r_tot = r[i] + r[j];
    m_tot = m[i] + m[j];
    gamma = m[i] / m[j];
    temp = r_tot / para->B_totseg;
    B[0] =  temp / 2;
    //printf("r_tot is %lf, B[0] is %lf\n", r_tot, B[0]);
    for (count = 1; count < para->B_totseg; count++) {
        B[count] = B[count-1] + temp;
        //printf("B[%d] is %lf \n", count, B[count]);
    }
    b_crit = r[j] / r_tot;
    xi = (m[j] - m[i]) / (m_tot);
    n_alpha = m[i] * m[j] / m_tot;
    
    //temp = 2 * Grav * (m[i] + m[j]) / (r[i] + r[j]) / pow(impact_velocity, 2);
    //CC[i][j] = Pi * (r[i]+r[j]) * (r[i]+r[j]) * impact_velocity * (1 + temp);

    //z = n[i] * n[j] * CC[i][j] * (double)(para->Timestep) / para->pVolume;
    
#ifdef Test_Safronov_Number
    double Focusing_number = para->v_frag_th;
    z = n[i] * n[j] * GCS[i][j] * impact_velocity * (1 + Focusing_number) * para->Timestep / para->pVolume;
#else
    z = n[i] * n[j] * GCS[i][j] * impact_velocity * (1 + MEV[i][j] / impact_velocity / impact_velocity) * para->Timestep / para->pVolume;
#endif
    
    if (i == j) z = z * 0.5;
    //printf("i = %ld, j = %ld, Focusing number is %e\n", i, j, MEV[i][j] / impact_velocity / impact_velocity);
    /*
    if (z != 0) {// && (i == 36 || j == 36)) {
        printf("n[%ld] = %e, n[%ld] = %e, V_d = %e, V_im = %e, z = %e, other = %e\n", i, n[i], j, n[j], v_d_relative[i][j], impact_velocity, z, z/para->Timestep);
    }
    //*/
    /*
    if (z > n[i] || z > n[j]) {
        if (n[i] < n[j]) {
            z = n[i];
        } else {
            z = n[j];
        }
        para->countnumber[7]++;
    }
    //*/
    delta->delta_n_coag[i] -= z;
    delta->delta_n_coag[j] -= z;
    
    check = 0;
    for (count = 0; count < para->B_totseg; count++) {
        
        //prob = ((B[count]+temp/2) * (B[count]+temp/2) - (B[count]-temp/2) * (B[count]-temp/2)) / r_tot / r_tot;
        prob = 2 * B[count] / r_tot / para->B_totseg;
        //printf("prob is %f \n", prob);
        //check += prob;
        tempz = z * prob;
        // first, calculate b, l, m_alpha and miu
        b = B[count] / r_tot;
        l = r_tot - B[count];
        if (r[j] > b*r_tot+r[i]) {
            m_alpha = 1;
        } else{
            m_alpha = (3*r[i]*l*l-l*l*l) / (4*r[i]*r[i]*r[i]);
            //printf("m_alpha = %f\n", m_alpha);
        }
        
        miu = m_alpha * m[i] * m[j] / (m_alpha * m[i] + m[j]);
        
        // and then, calculate mutual escape velocity
        V_esc = sqrt(2*Grav*(m_alpha * m[i] + m[j])/r_tot);
        if (impact_velocity < V_esc) {
            // perfect merge happens
            para->countnumber[2]++;
            k = findbin(m_tot);
            if (k >= para->Totalbin) {
                pMass_ghost += tempz * m_tot;
            } else {
                delta->delta_n_coag[k] += tempz * m_tot / m[k];
            }

            continue;
        }
        
        if (b > b_crit) {
            // calculate V_hr
            V_hr = V_esc * (c1 * xi * xi * pow(1-b, 5./2) + c2 * xi * xi + c3 * pow(1-b, 5./2) + c4);
            //printf("V_hr/V_esc = %f\n", V_hr/V_esc);
            if (impact_velocity < V_hr) {
                // perfect merge happens
                para->countnumber[6]++;
                k = findbin(m_tot);
                if (k >= para->Totalbin) {
                    pMass_ghost += tempz * m_tot;
                } else {
                    delta->delta_n_coag[k] += tempz * m_tot / m[k];
                }
            
            } else {
                delta->delta_n_coag[i] += tempz;
                delta->delta_n_coag[j] += tempz;
                // hit and run happens, here we only consider ideal hit and run
                // actually, if v_d_relative is too high, the projectile will disrupt
                para->countnumber[3]++;
            }
            //continue;

        } else {
            // fragmentation happens
            // calcualte Q_R and Q_RD
            Q_R = miu * impact_velocity / 2 / m_tot;
            Q_RD = para->c_star * 4/5 * Pi * Grav * rho1 * pow(m_tot / rho1, 2./3);
            Q_RD = Q_RD * pow(0.25 * (gamma+1) * (gamma+1) / gamma, 2/3/para->miu_bar-1) * pow(n_alpha / miu, 2-3*para->miu_bar/2);
            temp = Q_R / Q_RD;
            if (temp < 1.8) {
                M_lr = m_tot * 0.5 * (2 - temp);
            } else {
                M_lr = m_tot * 0.1 / pow(1.8, -1.5) * pow(temp, -1.5);
            }
            if (M_lr < m[j]) {
                para->countnumber[4]++;
            } else {
                //printf("%e\n", (m_tot-M_lr)/m_tot);
                para->countnumber[5]++;
            }

            
#define Common_Frag
//#define Frag_Upper_limit
//#define Frag_Lower_limit
            
#ifdef Common_Frag
            if (M_lr < 0.5 * m_tot) {
                M_second = 0.5 * M_lr;
                M_third = 0.5 * M_second;
            } else {
                M_second = 0.5 * (m_tot - M_lr);
                M_third = 0.5 * M_second;
            }

            lr_bin = findbin(M_lr);
            second_bin = findbin(M_second);
            third_bin = findbin(M_third);
            
            if (lr_bin < 0) {
                pMass_back += tempz * M_lr;
            } else if (lr_bin >= para->Totalbin) {
                pMass_ghost += tempz * M_lr;
            } else {
                delta->delta_n_frag[lr_bin] += tempz * M_lr / m[lr_bin];
            }
            if (second_bin < 0) {
                pMass_back += tempz * M_second;
            } else if (second_bin >= para->Totalbin) {
                pMass_ghost += tempz * M_second;
            } else {
                delta->delta_n_frag[second_bin] += tempz * M_second / m[second_bin];
            }
            if (third_bin < 0) {
                pMass_back += tempz * M_third;
            } else if (third_bin >= para->Totalbin) {
                pMass_ghost += tempz * M_third;
            } else {
                delta->delta_n_frag[third_bin] += tempz * M_third / m[third_bin];
            }
            pMass_back += tempz * (m_tot - M_lr - M_second -M_third);
            
#endif
#ifdef Frag_Upper_limit
            if (M_lr > 0.5 * m_tot) {
                M_second = m_tot - M_lr;
            } else {
                M_second = 0.5 * M_lr;
            }
            lr_bin = findbin(M_lr);
            second_bin = findbin(M_second);
            if (lr_bin < 0) {
                pMass_back += tempz * M_lr;
            } else if (lr_bin >= para->Totalbin) {
                pMass_ghost += tempz * M_lr;
            } else {
                delta->delta_n_frag[lr_bin] += tempz * M_lr / m[lr_bin];
            }
            if (second_bin < 0) {
                pMass_back += tempz * M_second;
            } else if (second_bin >= para->Totalbin) {
                pMass_ghost += tempz * M_second;
            } else {
                delta->delta_n_frag[second_bin] += tempz * M_second / m[second_bin];
            }
            pMass_back += tempz * (m_tot - M_lr - M_second);
            
#endif
#ifdef Frag_Lower_limit
            lr_bin = findbin(M_lr);
            if (lr_bin < 0) {
                pMass_back += tempz * M_lr;
            } else if (lr_bin >= para->Totalbin) {
                pMass_ghost += tempz * M_lr;
            } else {
                delta->delta_n_frag[lr_bin] += tempz * M_lr / m[lr_bin];
            }
            pMass_back += tempz * (m_tot - M_lr);
            
#endif
            
            
            
        }
        
    }
    
    //printf("total prob is %f.\n", check);
    free(B);
    
    return 0;
}