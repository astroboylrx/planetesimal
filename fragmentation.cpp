//
//  fragmentation.cpp
//  planetesimal
//
//  Created by Rixin Li on 8/20/13.
//  Copyright (c) 2013 李 日新. All rights reserved.
//

#include "fragmentation.h"


int fragmentation(long &i, long &j)
/** calculate the fragmentation between two mass bins **/
{
    double z, mass_involved = 0, impact_velocity, temp;
    long i_result, j_result;
    
    if (n[i] == 0 || n[j] == 0) {
        return 0;
    }
    
    // we'll use the same collision kernel and set a threshold velocity for fragmentation
    // v_d_relative = pow((v_d_set[i]-v_d_set[j])*(v_d_set[i]-v_d_set[j])+v_d_tur[i]*v_d_tur[i]+v_d_tur[j]*v_d_tur[j], 0.5);
    
    // here we apply the probability (from BD2008) times collision kernel to calculate the amount of fragmentation
    impact_velocity = v_d_relative[i][j] * sqrt(-2*log(1-(rand()%100000000+0.5/100000000)/100000000.0));
    while (impact_velocity/v_d_relative[i][j] > 1.66511 || impact_velocity/v_d_relative[i][j] < 0.758528 ) {
        //printf("Processor %d: calculate impact velocity again\n", myid);
        para->countnumber[9]++;
        impact_velocity = v_d_relative[i][j] * sqrt(-2*log(1-(rand()%100000000+0.5/100000000)/100000000.0));
    }
    
    temp = para->v_frag_th - impact_velocity;
    temp = (impact_velocity / para->v_frag_th) * Heaviside(temp) + Heaviside(0-temp);
    //z = n[i] * n[j] * Fij[i][j] * (double)(para->Timestep) / para->pVolume;
    z = n[i] * n[j] * GCS[i][j] * impact_velocity * (1 + MEV[i][j] / impact_velocity / impact_velocity) * para->Timestep / para->pVolume * temp;
    if (i == j) z = z * 0.5;
    //printf("i = %ld, j = %ld, z_frag = %e\n", i, j, z);   //if enable this statement, there would be tricky
    /*
    if (z > n[i] || z > n[j]) {
        if (n[i] < n[j]) {
            z = n[i];
        } else {
            z = n[j];
        }
        para->countnumber[1]++;
    }
    */
    delta->delta_n_frag[i] -= z;
    delta->delta_n_frag[j] -= z;
    
    mass_involved += z * (m[i] + m[j]);
    

    // if the relative velocity dispersion exceed another threshold velocity for disruption,
    // we think this collision will be a catastrophic disruption
    if (v_d_relative[i][j] > para->v_disr_th && r[i]/r[j] > 0.1) {
        
        // We all keep the largest post-collision remnant and
        // let the rest mass become dust of background
        // Now we will treat the M_lr as 0.3M_tot
        // Later we may using the result of L2012
        
        long k;
        double m_lr;
        m_lr = para->Q_lr * (m[i] + m[j]);
        k = findbin(m_lr);
        if (k >= para->Totalbin) {
            pMass_ghost += para->Q_lr * mass_involved;
        } else {
            delta->delta_n_frag[k] = para->Q_lr * mass_involved / m[k];
        }
        pMass_back += (1 - para->Q_lr) * mass_involved;
        para->countnumber[1]++;
        return 0;
        
    }
   
    // for simplicity, we assume that fragmentation result is that
    // each planetesimals involved in the collision break up into
    // two smaller equal-mass planetesimals

    // Besides, if their radius difference larger than 10 times, then
    // we regard it as merging or cratering
    
    if (r[i]/r[j] > 0.1) {
        temp = m[i]/2.0;
        i_result = findbin(temp);
        if (i_result < 0) {
            pMass_back += z * m[i];
        } else {
            delta->delta_n_frag[i_result] += z * m[i] / m[i_result];
        }
        temp = m[j]/2.0;
        j_result = findbin(temp);
        if (j_result < 0) {
            pMass_back += z * m[j];
        } else {
            delta->delta_n_frag[j_result] += z * m[j] / m[j_result];
        }
    } else {
        
        //temp = ((r[j]+r[i])*(r[j]+r[i]) - (r[j]-r[i])*(r[j]-r[i])) / (r[j]+r[i]) / (r[j]+r[i]);
        temp = (2 * r[j] * 2 * r[i]) / (r[j]+r[i]) / (r[j]+r[i]);
        // cratering happens on a possibility
        delta->delta_n_frag[j] += z * temp;
        // some of them becomes background dust
        pMass_back += z * temp * m[i];
        // using (m[i]-m[j]) is in order to avoid non-conservation of mass
        
        // merging happens also on a possibility
        j_result = findbin(m[i]+m[j]);
        if (j_result >= para->Totalbin) {
            pMass_ghost += z * (1 - temp) * (m[j] + m[i]);
        } else {
            delta->delta_n_frag[j_result] += z * (1 - temp) * (m[j] + m[i]) / m[j_result];
        }
        

    }
    
    
    return 0;
}