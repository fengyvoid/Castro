#!/bin/bash

# echo the commands
set -x

DIM=3
EXEC=./Castro${DIM}d.gnu.MPI.ex

RUNPARAMS="
castro.ppm_type=0
castro.time_integration_method=0
problem.init_as_1d=1
geometry.prob_hi=1 0.125 0.125"

mpiexec -n 8 ${EXEC} inputs.64 amr.n_cell=64 8 8 ${RUNPARAMS} amr.plot_file=acoustic_pulse_64_plm_plt &> 64.out
mpiexec -n 8 ${EXEC} inputs.128 amr.n_cell=128 16 16  ${RUNPARAMS} amr.plot_file=acoustic_pulse_128_plm_plt &> 128.out
mpiexec -n 8 ${EXEC} inputs.256 amr.n_cell=256 32 32 ${RUNPARAMS} amr.plot_file=acoustic_pulse_256_plm_plt &> 256.out

RichardsonConvergenceTest${DIM}d.gnu.ex coarFile=acoustic_pulse_64_plm_plt00081 mediFile=acoustic_pulse_128_plm_plt00161 fineFile=acoustic_pulse_256_plm_plt00321 > convergence.${DIM}d.lo.plm.x.out

#mpiexec -n 4 ${EXEC} inputs.512 ${RUNPARAMS} amr.plot_file=acoustic_pulse_512_plm_plt &> 512.out

#RichardsonConvergenceTest${DIM}d.gnu.ex coarFile=acoustic_pulse_128_plm_plt00161 mediFile=acoustic_pulse_256_plm_plt00321 fineFile=acoustic_pulse_512_plm_plt00641 > convergence.${DIM}d.hi.plm.out

