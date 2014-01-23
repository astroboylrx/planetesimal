//
//  cog_frag_ingrav.cpp
//  planetesimal
//
//  Created by Rixin Li on 11/19/13.
//  Copyright (c) 2013 李 日新. All rights reserved.
//

#include "cog_frag_ingrav.h"

int cog_frag_ingrav(long &i, long &j)
/** calculate the both coagulation and fragmentation in gravity regime between two mass bins **/
{
    // The basic idea of this subroutine is from LS2012
    //
    // we simply assume that for the whole collisional cross section of target, the projectile will hit any point
    // with the same probability, so we divide B into several segments
    // l + B = R + r
    
    int count;
    long lr_bin, k;
    double *B, temp, b, l, b_crit, r_tot, m_tot, gamma, n_alpha, m_alpha, miu, V_esc, V_hr;
    double prob, xi, c1=2.43, c2=-0.0408, c3=1.86, c4=1.08;
    double Q_R, Q_RD, rho1 = 1.0e3, M_lr, z, tempz, check;
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
    
    z = partn[i][j] * partn[j][i] * CC[i][j] * (double)(para->Timestep) / para->pVolume;
    if (i == j) z = z * 0.5 * 0.5;
    if (z > partn[i][j] || z > partn[j][i]) {
        if (partn[i][j] < partn[j][i]) {
            z = partn[i][j];
        } else {
            z = partn[j][i];
        }
    }
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
            m_alpha = para->Rho_solid * Pi * l * l * (r[i] - l / 3.0);
        }
        
        miu = m_alpha * m[i] * m[j] / (m_alpha * m[i] + m[j]);
        
        // and then, calculate mutual escape velocity
        V_esc = sqrt(2*Grav*(m_alpha * m[i] + m[j])/r_tot);
        if (v_d_relative[i][j] < V_esc) {
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
            if (v_d_relative[i][j] < V_esc) {
                // perfect merge happens
                para->countnumber[2]++;
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
            Q_R = miu * v_d_relative[i][j] / 2 / m_tot;
            Q_RD = para->c_star * 4/5 * Pi * Grav * rho1 * pow(m_tot / rho1, 2./3);
            Q_RD = Q_RD * pow(0.25 * (gamma+1) * (gamma+1) / gamma, 2/3/para->miu_bar-1) * pow(n_alpha / m_alpha, 2-3*para->miu_bar/2);
            temp = Q_R / Q_RD;
            if (temp < 1.8) {
                M_lr = m_tot * 0.5 * (2 - temp);
            } else {
                M_lr = m_tot * 0.1 / pow(1.8, -1.5) * pow(temp, -1.5);
            }
            if (M_lr < m[j]) {
                para->countnumber[4]++;
            } else {
                para->countnumber[5]++;
            }
            
            lr_bin = findbin(M_lr);
            
            if (lr_bin < 0) {
                pMass_back += tempz * m_tot;
            } else if (lr_bin >= para->Totalbin) {
                pMass_ghost += tempz * M_lr;
                pMass_back += tempz * (m_tot - M_lr);
            } else {
                delta->delta_n_frag[lr_bin] += tempz;
                pMass_back += tempz * (m_tot - m[lr_bin]);
            }
            
        }
        
    }
    
    //printf("total prob is %f.\n", check);
    free(B);
    
    return 0;
}