//
//  condensation.cpp
//  planetesimal
//
//  Created by Rixin Li on 8/20/13.
//  Copyright (c) 2013 李 日新. All rights reserved.
//

#include "condensation.h"


int condensation()
/** calculate the influence of condensation **/
{
    // for simplicity, since it would take 10^3 yrs for micro-size dust
    // to settle and growth to centimeter-size particles, now we just
    // set a ratio for the background materials to condense
    // assuming that it takes the same duration for centimeter-size particles
    // to grow to meter-size boulders
    //
    // based on the former calculation, if the dust-to-gas ratio reach 0.8, then
    // due to gravitational instabilities, 80 M_Earth solid materials will produce
    // ~2.5e+06 planetesimals (radius ~ 195km), leaving about half mass in background.
    // for the simplisity of simulaiton, we convert the total mass of planetesimals
    // into 100km-planetesimals, so we got 1.85e+07 planetesimals
    //
    double cond_ratio;
    
    if (r[0] < para->R_in_grav_regime) {
        cond_ratio = (para->Timestep/yr/2000.); // (+ 0.01/r[0])?
        delta->delta_n_cond[0] = pMass_back * cond_ratio / m[0];
        pMass_back *= (1 - cond_ratio);
        pMass_back += para->pAcc_intoback;
        
    } else {
        cond_ratio = 1.85e+07/(para->gMass_tot_i*0.8/para->pAcc_rate);
        delta->delta_n_cond[0] = cond_ratio;
        //pMass_back -= cond_ratio * m[0] * para->Timestep/yr;
        pMass_back = pMass_back + para->pAcc_rate - cond_ratio * m[0];
        //printf("cond_ratio * m[0] = %e\n", cond_ratio * m[0]);
        
    }
    
    
    return 0;
}