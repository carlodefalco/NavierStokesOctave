# NavierStokesOctave
Octave olver for 2D incompressible Navier-Stokes equations on a Cartesian grid using Inf-Sup stable Q2-P1 elements [[1]](#1)

## Stationary solver

Solve 

$$ - \nu\ (\partial_{xx} u_x + \partial_{yy} u_x) + u_x\ \partial_x u_x + u_y\ \partial_y u_x + \partial_x\ p = f_x $$

$$ - \nu\ (\partial_{xx} u_y + \partial_{yy} u_y) + u_x\ \partial_{x} u_y + u_y\ \partial_{y} u_y + \partial_{y}\ p  = f_y $$

constrained by

$$ \partial_{x} u_x + \partial_{y} u_y = 0. $$

Initialize by solving the Stokes problem

$$ - \nu\ (\partial_{xx} u_x + \partial_{yy} u_x) + \partial_x\ p    = f_x $$

$$ - \nu\ (\partial_{xx} u_y + \partial_{yy} u_y) + \partial_{y}\ p  = f_y $$

constrained by

$$ \partial_{x} u_x + \partial_{y} u_y = 0. $$

Then perform fiexed-point iterations

$$ - \nu\ (\partial_{xx} u_x^{(k+1)} + \partial_{yy} u_x^{(k+1)}) + u_x^{(k)}\ \partial_x u_x^{(k)}  + u_y^{(k)}\ \partial_y u_x^{(k)} + \partial_x\ p^{(k+1)} = f_x $$

$$ - \nu\ (\partial_{xx} u_y^{(k+1)} + \partial_{yy} u_y^{(k+1)}) + u_x^{(k)}\ \partial_{x} u_y^{(k)} + u_y^{(k)}\ \partial_{y} u_y^{(k)} + \partial_{y}\ p^{(k+1)}  = f_y $$

$$ \partial_{x} u_x^{(k+1)} + \partial_{y} u_y^{(k+1)} = 0. $$

## References
<a id="1">[1]</a> 
Daniele Boffi , Franco Brezzi, Michel Fortin.
Mixed Finite Element Methods and Applications.
Sprimger. Ch. 8 Sec. 6
