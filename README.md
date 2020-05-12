# Reference-interaction site-model (RISM)-solvers

## Publisher
Author: Tsogbayar Tsednee (PhD)
Contact: tsog215@gmail.com

## Introduction:

The rism1d_solver_homonuclear_diatomic_liquid_hnc.m solves the RISM equation for uncharged homonuclear diatomic liquid in
the hypernetted chain approximation using the Picard iterative method. 

## Requirement:
Any version of Matlab

## Implementation details and running

The code oz_lj_ts_run_opt_one_parm.m uses the fminsearch.m function from Matlab software to find a optimal value of a free parameter. A code oz_lj_ts.m solves the Ornstein-Zernike integral equation for the Lennard-Jones potential together with a one-parameter Verlet-modified closure. The code oz_lj_ts.m  uses a code lsint.m which calculates the sine transform of a function.
You may download them and run the oz_lj_ts_run_opt_one_parm.m directly. You can download and run the OZ_solver_one_component_plasma_hnc.m directly.

## Copyright / License

These codes are published as freeware. Basically, you have the right to use them and modify any of them.

Find the GNU General Public License at:
https://www.gnu.org/licenses/gpl-3.0.en.html

### Acknowledgments
Author is very thankful for advices of Professor T. Luchko at California State University, Northridge regarding to an usage of a fminsearch function of Matlab software.
