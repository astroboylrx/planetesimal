//
//  coagulation.cpp
//  planetesimal
//
//  Created by Rixin Li on 8/20/13.
//  Copyright (c) 2013 李 日新. All rights reserved.
//

#include "coagulation.h"

int coagulation(long &i, long &j)
/** calculate the coagulation between two mass bins **/
{
    long k;
	double z, temp;
    
    if (n[i] == 0 || n[j] == 0) {
        return 0;
    }
    // using the coagulation equation to calculate the result
    // Aij denotes the relative probability of collision, named collision kernel
    // Aij = Pi * (r[i]+r[j]) * (r[i]+r[j]) * v_d_relative;
    
    z = partn[i][j] * partn[j][i] * Cij[i][j] * (double)(para->Timestep) / para->pVolume;

    //printf("i = %ld, j = %ld, z_coag = %e\n", i, j, z);  //if enable this statement, there would be tricky
    //printf("ratio[i] = %f, ratio[j] = %f, CC[i][j] = %f, v_d[i][j] = %f\n", pArea_ratio[i], pArea_ratio[j], CC[i][j], v_d_relative[i][j]);
    if (i == j) z = z * 0.5 * 0.5;
    
    if (z > partn[i][j] || z > partn[j][i]) {
        if (partn[i][j] < partn[j][i]) {
            z = partn[i][j];
        } else {
            z = partn[j][i];
        }
        para->countnumber[0]++;
        //z = ((partn[i][j]<partn[j][i])?partn[i][j]:partn[j][i]);
        //printf("i = %ld, j = %ld, z_coag changes to %e\n", i, j, z);
    }
    delta->delta_n_coag[i] -= z;
    delta->delta_n_coag[j] -= z;
    
	temp = m[i] + m[j]; /*the result mass assuming no mass loss during the
                         coagulation*/
    /*if it ends up in our ghost zone*/
    k = findbin(temp);
    if (k >= para->Totalbin) {
        pMass_ghost += z * temp;
    } else {
        /*otherwise, it ends up in bin k*/
        /*Record the mass change in bin k due to this coagulation*/
        delta->delta_n_coag[k] += z * temp / m[k];
    }
        
	return 0;
}
