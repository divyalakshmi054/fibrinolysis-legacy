// S-law model of fibrinolysis
// Author: dl653@cornell.edu

#static::start
    tPA,PAI1,TAFI,AT,FIIa
#static::end

#dynamic::start
    FI,FIa,FDP,Plgn,Plasmin
#dynamic::end

#structure::start
    // name::{reactants} --> {products}
    // 1. name: reaction name (used in rate section)
    // 2. {reactants}: list of dynamic species acting as a reactants
    // 3. {products}: list of dynamic species acting as product
    
    // reactions in the pathway -
    r1::{FI} --> {FIa}
    r2::{Plgn} --> {Plasmin}
    r3::{FIa} --> {FDP}
#structure::end

#rate::start
    // name::{factor_1,factor_2,....}
    // 1: name = reaction name (must be in the #pragma::structure section)
    // 2: list of factor symbols (can be both static and dynamic, order does not matter)
    // 3: the order of the records does not matter here, we establish the order from the structure section
    r1::{FI,FIIa}
    r2::{Plgn,tPA,PAI1}
    r3::{FIa,Plasmin,TAFI}
#rate::end