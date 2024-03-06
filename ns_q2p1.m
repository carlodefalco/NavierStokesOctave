close all
clear all

## ## Stationary solver
##
## Solve 
##
## $$ - \nu\ (\partial_{xx} u_x + \partial_{yy} u_x) + u_x\ \partial_x u_x + u_y\ \partial_y u_x + \partial_x\ p = f_x $$
##
## $$ - \nu\ (\partial_{xx} u_y + \partial_{yy} u_y) + u_x\ \partial_{x} u_y + u_y\ \partial_{y} u_y + \partial_{y}\ p  = f_y $$
##
## constrained by
##
## $$ \partial_{x} u_x + \partial_{y} u_y = 0. $$
##
## Initialize by solving the Stokes problem
##
## $$ - \nu\ (\partial_{xx} u_x + \partial_{yy} u_x) + \partial_x\ p    = f_x $$
##
## $$ - \nu\ (\partial_{xx} u_y + \partial_{yy} u_y) + \partial_{y}\ p  = f_y $$
##
## constrained by
##
## $$ \partial_{x} u_x + \partial_{y} u_y = 0. $$
##
## Then perform fiexed-point iterations
##
## $$ - \nu\ (\partial_{xx} u_x^{(k+1)} + \partial_{yy} u_x^{(k+1)}) + u_x^{(k)}\ \partial_x u_x^{(k)}  + u_y^{(k)}\ \partial_y u_x^{(k)} + \partial_x\ p^{(k+1)} = f_x $$
##
## $$ - \nu\ (\partial_{xx} u_y^{(k+1)} + \partial_{yy} u_y^{(k+1)}) + u_x^{(k)}\ \partial_{x} u_y^{(k)} + u_y^{(k)}\ \partial_{y} u_y^{(k)} + \partial_{y}\ p^{(k+1)}  = f_y $$
##
## $$ \partial_{x} u_x^{(k+1)} + \partial_{y} u_y^{(k+1)} = 0. $$





## Example problem, see https://doi.org/10.1002/fld.2337 for solutions (without convection)

nr = 25;
nc = 25;
hx = 1 / nc;
hy = 1 / nr;

A  = op_gradu_gradv (nr, nc);
M  = op_u_v (hx, hy, nr, nc);
Bx = op_dudx_p (hx, hy, nr, nc);
By = op_dudy_p (hx, hy, nr, nc);
 
[X, Y] = meshgrid (linspace (0, 1, 2*nc+1), linspace (0, 3, 2*nr+1));
ndof_ux = ndof_uy = numel (X);
dnodes = find ((X(:) <= eps) ...
	       | (X(:) >= (1)*(1-eps)) ...
	       | (Y(:) <= eps) ...
	       | (Y(:) >= (3)*(1-eps)));

ndof_p = nr*nc*3;
nu  = 1/32;
LHS = [nu*A       0*A      -Bx;
       0*A        nu*A     -By;
       -Bx.'      -By.'     sparse(columns(Bx), columns(Bx))];

dnodes_stokes = union (dnodes, dnodes+ndof_ux);
inodes_stokes = setdiff (1:(ndof_ux+ndof_uy+ndof_p), dnodes_stokes);

U = zeros (ndof_ux, 1);
U(Y(:) >= (3)*(1-eps)) = 1;
U(Y(:) <= eps)         = -1;
V = zeros (ndof_uy, 1);
p = zeros (ndof_p,  1);

SOL= [U;V;p];
RHS = zeros (ndof_ux+ndof_uy+ndof_p, 1);

SOL(inodes_stokes) = LHS(inodes_stokes, inodes_stokes) \ (RHS(inodes_stokes) - LHS(inodes_stokes, dnodes_stokes) * SOL(dnodes_stokes));
U = SOL (1 : ndof_ux);
V = SOL ((1 : ndof_uy) + ndof_ux);
p = SOL ((1 : ndof_p) + ndof_ux  + ndof_uy);

figure (1)
streamline (X, Y, reshape (U, size (X)), reshape (V, size (X)), linspace(0, 1, 10), linspace(0, 3, 10));
hold all
quiver (X, Y, reshape (U, size (X)), reshape (V, size (X)));
axis image

## Now enable convection and perform fixed point iterations

Cx = op_u0_gradu (hx, hy, nr, nc, U, V, U);
Cy = op_u0_gradu (hx, hy, nr, nc, U, V, V);

RHS(1:ndof_ux) = -Cx;
RHS((1:ndof_uy)+ndof_ux) = -Cy;

SOL(inodes_stokes) = LHS(inodes_stokes, inodes_stokes) \ (RHS(inodes_stokes) - LHS(inodes_stokes, dnodes_stokes) * SOL(dnodes_stokes));
U = SOL (1 : ndof_ux);
V = SOL ((1 : ndof_uy) + ndof_ux);
p = SOL ((1 : ndof_p) + ndof_ux  + ndof_uy);

figure (2)
streamline (X, Y, reshape (U, size (X)), reshape (V, size (X)), [linspace(0, 1, 10), linspace(1, 0, 10)], [linspace(0, 3, 10), linspace(0, 3, 10)]);
hold all
quiver (X, Y, reshape (U, size (X)), reshape (V, size (X)));
axis image

Cx = op_u0_gradu (hx, hy, nr, nc, U, V, U);
Cy = op_u0_gradu (hx, hy, nr, nc, U, V, V);

RHS(1:ndof_ux) = -Cx;
RHS((1:ndof_uy)+ndof_ux) = -Cy;

SOL(inodes_stokes) = LHS(inodes_stokes, inodes_stokes) \ (RHS(inodes_stokes) - LHS(inodes_stokes, dnodes_stokes) * SOL(dnodes_stokes));
U = SOL (1 : ndof_ux);
V = SOL ((1 : ndof_uy) + ndof_ux);
p = SOL ((1 : ndof_p) + ndof_ux  + ndof_uy);

figure (3)
streamline (X, Y, reshape (U, size (X)), reshape (V, size (X)), [linspace(0, 1, 10), linspace(1, 0, 10)], [linspace(0, 3, 10), linspace(0, 3, 10)]);
hold all
quiver (X, Y, reshape (U, size (X)), reshape (V, size (X)));
axis image
hold off
drawnow

Cx = op_u0_gradu (hx, hy, nr, nc, U, V, U);
Cy = op_u0_gradu (hx, hy, nr, nc, U, V, V);

RHS(1:ndof_ux) = -Cx;
RHS((1:ndof_uy)+ndof_ux) = -Cy;

SOL(inodes_stokes) = LHS(inodes_stokes, inodes_stokes) \ (RHS(inodes_stokes) - LHS(inodes_stokes, dnodes_stokes) * SOL(dnodes_stokes));
U = SOL (1 : ndof_ux);
V = SOL ((1 : ndof_uy) + ndof_ux);
p = SOL ((1 : ndof_p) + ndof_ux  + ndof_uy);

close (3)
figure (3)
streamline (X, Y, reshape (U, size (X)), reshape (V, size (X)), [linspace(0, 1, 10), linspace(1, 0, 10)], [linspace(0, 3, 10), linspace(0, 3, 10)]);
hold all
quiver (X, Y, reshape (U, size (X)), reshape (V, size (X)));
axis image
hold off
drawnow

Cx = op_u0_gradu (hx, hy, nr, nc, U, V, U);
Cy = op_u0_gradu (hx, hy, nr, nc, U, V, V);

RHS(1:ndof_ux) = -Cx;
RHS((1:ndof_uy)+ndof_ux) = -Cy;

SOL(inodes_stokes) = LHS(inodes_stokes, inodes_stokes) \ (RHS(inodes_stokes) - LHS(inodes_stokes, dnodes_stokes) * SOL(dnodes_stokes));
U = SOL (1 : ndof_ux);
V = SOL ((1 : ndof_uy) + ndof_ux);
p = SOL ((1 : ndof_p) + ndof_ux  + ndof_uy);

close (3)
figure (3)
streamline (X, Y, reshape (U, size (X)), reshape (V, size (X)), [linspace(0, 1, 10), linspace(1, 0, 10)], [linspace(0, 3, 10), linspace(0, 3, 10)]);
hold all
quiver (X, Y, reshape (U, size (X)), reshape (V, size (X)));
axis image
hold off
drawnow

Cx = op_u0_gradu (hx, hy, nr, nc, U, V, U);
Cy = op_u0_gradu (hx, hy, nr, nc, U, V, V);

RHS(1:ndof_ux) = -Cx;
RHS((1:ndof_uy)+ndof_ux) = -Cy;

SOL(inodes_stokes) = LHS(inodes_stokes, inodes_stokes) \ (RHS(inodes_stokes) - LHS(inodes_stokes, dnodes_stokes) * SOL(dnodes_stokes));
U = SOL (1 : ndof_ux);
V = SOL ((1 : ndof_uy) + ndof_ux);
p = SOL ((1 : ndof_p) + ndof_ux  + ndof_uy);

close (3)
figure (3)
streamline (X, Y, reshape (U, size (X)), reshape (V, size (X)), [linspace(0, 1, 10), linspace(1, 0, 10)], [linspace(0, 3, 10), linspace(0, 3, 10)]);
hold all
quiver (X, Y, reshape (U, size (X)), reshape (V, size (X)));
axis image
hold off
drawnow
