//
//  disruption.cpp
//  planetesimal
//
//  Created by Rixin Li on 8/20/13.
//  Copyright (c) 2013 李 日新. All rights reserved.
//

#include "disruption.h"


int disruption(long &i, long &j, double &mass_involved)
/** calculate the result of catastrophic disruption **/
{
    // We all keep the largest post-collision remnant and
    // let the rest mass become dust of background
    // Now we will treat the M_lr as 0.3M_tot
    // Later we may using the result of L2012
    
    long k;
    double m_lr;
    m_lr = para->Q_lr * (m[i] + m[j]);
    
    k = j;
    while (k > 0 && m[k] > m_lr) {
        k--;
    }
    delta->delta_n_disr[k] = para->Q_lr * mass_involved / m[k];
    pMass_back += (1 - para->Q_lr) * mass_involved;
    
    return 0;
}