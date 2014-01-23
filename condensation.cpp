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
    double cond_ratio;
    cond_ratio = (para->Timestep/yr/2000.); // (+ 0.01/r[0])?
    delta->delta_n_cond[0] = pMass_back * cond_ratio / m[0];
    pMass_back *= (1 - cond_ratio);
    
    return 0;
}