//
//  fragmentation.cpp
//  planetesimal
//
//  Created by Rixin Li on 8/20/13.
//  Copyright (c) 2013 李 日新. All rights reserved.
//

#include "fragmentation.h"
#include "disruption.h"

int fragmentation(long &i, long &j)
/** calculate the fragmentation between two mass bins **/
{
    double z, mass_involved = 0;
    long i_result, j_result, temp;
    
    if (n[i] == 0 || n[j] == 0) {
        return 0;
    }
    
    // we'll use the same collision kernel and set a threshold velocity for fragmentation
    // v_d_relative = pow((v_d_set[i]-v_d_set[j])*(v_d_set[i]-v_d_set[j])+v_d_tur[i]*v_d_tur[i]+v_d_tur[j]*v_d_tur[j], 0.5);
    
    // here we apply the probability (from BD2008) times collision kernel to calculate the amount of fragmentation

    z = partn[i][j] * partn[j][i] * Fij[i][j] * (double)(para->Timestep) / para->pVolume;
    if (i == j) z = z * 0.5 * 0.5;
    //printf("i = %ld, j = %ld, z_frag = %e\n", i, j, z);   //if enable this statement, there would be tricky
    
    // remove mass from the original bin
    // and add them to the new bin
    if (z > partn[i][j] || z > partn[j][i]) {
        if (partn[i][j] < partn[j][i]) {
            z = partn[i][j];
        } else {
            z = partn[j][i];
        }
        para->countnumber[1]++;
    }
    delta->delta_n_frag[i] -= z;
    delta->delta_n_frag[j] -= z;
    
    mass_involved += z * (m[i] + m[j]);
    

    // if the relative velocity dispersion exceed another threshold velocity for disruption,
    // we think this collision will be a catastrophic disruption
    if (v_d_relative[i][j] > para->v_disr_th && r[i]/r[j] > 0.1) {
        disruption(i, j, mass_involved);
        return 0;
    }
   
    // for simplicity, we assume that fragmentation result is that
    // each planetesimals involved in the collision break up into
    // two smaller equal-mass planetesimals
    
    temp = m[i]/2.0;
    i_result = findbin(temp);
    if (i_result < 0) {
        pMass_back += z * m[i];
    } else {
        delta->delta_n_frag[i_result] += z * m[i] / m[i_result];
    }

    // Besides, if their radius difference larger than 10 times, then
    // we regard it as merging or cratering
    
    if (r[i]/r[j] > 0.1) {
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