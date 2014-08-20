//
//  coagulation.cpp
//  planetesimal
//
//  Created by Rixin Li on 8/20/13.
//  Copyright (c) 2013 李 日新. All rights reserved.
//

#include "coagulation.h"


int coagulation(long &i, long &j)
/** calculate the coagulation between two small-mass bins **/
{
    long k;
	double z, temp, impact_velocity;
    
    if (n[i] == 0 || n[j] == 0) {
        return 0;
    }
    // using the coagulation equation to calculate the result
    // Aij denotes the relative probability of collision, named collision kernel
    // Aij = Pi * (r[i]+r[j]) * (r[i]+r[j]) * v_d_relative;

    impact_velocity = v_d_relative[i][j] * sqrt(-2*log(1-(rand()%100000000+0.5/100000000)/100000000.0));
    while (impact_velocity/v_d_relative[i][j] > 1.66511 || impact_velocity/v_d_relative[i][j] < 0.758528 ) {
        //printf("Processor %d: calculate impact velocity again\n", myid);
        para->countnumber[9]++;
        impact_velocity = v_d_relative[i][j] * sqrt(-2*log(1-(rand()%100000000+0.5/100000000)/100000000.0));
    }
    
    temp = para->v_frag_th - impact_velocity;
    temp = (impact_velocity / para->v_frag_th) * Heaviside(temp) + Heaviside(0-temp);
    //z = n[i] * n[j] * Cij[i][j] * para->Timestep / para->pVolume;
    z = n[i] * n[j] * GCS[i][j] * impact_velocity * (1 + MEV[i][j] / impact_velocity / impact_velocity) * para->Timestep / para->pVolume * (1 - temp);
    
    if (i == j) z = z * 0.5;
    //printf("i = %ld, j = %ld, z_coag = %e\n", i, j, z);  //if enable this statement, there would be tricky    

    /*
    if (z > n[i] || z > n[j]) {
        if (n[i] < n[j]) {
            z = n[i];
        } else {
            z = n[j];
        }
        para->countnumber[0]++;
    }   
     */
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
