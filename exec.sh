#!/bin/bash

# restart simulation or not
export restart_sim=false
# if restart_sim=true, restart_dir is mandatory
# normally information exist under "output/save/stepxxxxxxxxxx"
# export restart_dir=output/save/stepxxxxxxxxxx

## flags to specify how to treat temperature field
# solve temperature field or not
export solve_temp=true
# buoyancy force is added to x-momentum or not
export add_buoyancy=true

## durations
# maximum duration (in free-fall time)
export timemax=2.0e+2
# maximum duration (in wall time [s])
export wtimemax=6.0e+2
# logging rate (in free-fall time)
export log_rate=2.0e-1
# logging after (in free-fall time)
export log_after=0.0e+0
# save rate (in free-fall time)
export save_rate=2.0e+1
# save after (in free-fall time)
export save_after=0.0e+0
# statistics collection rate (in free-fall time)
export stat_rate=1.0e-1
# statistics collection after (in free-fall time)
export stat_after=1.0e+2

## domain
# domain lengths
export lx=1.0e+0
export ly=2.0e+0
# number of cell centers
export glisize=128
export gljsize=256

## treatment of diffusive terms
export implicitx=false
export implicity=false

## safety factors to decide time step size
## for advective and diffusive terms
export coef_dt_adv=0.95
export coef_dt_dif=0.95

## physical parameters
export Ra=1.0e+8
export Pr=1.0e+1

mpirun -n 2 --oversubscribe ./a.out
