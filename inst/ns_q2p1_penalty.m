1;

## ## Stationary solver, rectangular domain with cylinder
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



## Example problem, flow past cylinder

nr = 25;
nc = 75;
hx = 1 / nc;
hy = 1 / nr;

L =  3;
H =  1;
r =  L/10;

A  = op_gradu_gradv (nr, nc);
M  = op_u_v (hx, hy, nr, nc);
Bx = op_dudx_p (hx, hy, nr, nc);
By = op_dudy_p (hx, hy, nr, nc);
 
[X, Y] = meshgrid (linspace (0, L, 2*nc+1), linspace (0, H, 2*nr+1));
ndof_ux = ndof_uy = numel (X);
dnodes = find ((X(:) <= eps) ...
	       | (Y(:) <= eps) ...
	       | (Y(:) >= (H)*(1-eps)) ...
	       | ( ((X(:)-L/3).^2 + Y(:).^2) <= r^2));

ndof_p = nr*nc*3;
nu  = 1/10.;
LHS = [nu*A                               sparse(rows(A), columns(A))      -Bx;
       sparse(rows(A), columns(A))        nu*A                             -By;
       -Bx.'                              -By.'                 sparse(columns(Bx), columns(Bx))];

dnodes_stokes = union (dnodes, dnodes+ndof_ux);
inodes_stokes = setdiff (1:(ndof_ux+ndof_uy+ndof_p), dnodes_stokes);

U = zeros (ndof_ux, 1);
U(X(:) <= eps)    = 1;
V = zeros (ndof_uy, 1);
p = zeros (ndof_p,  1);

SOL= [U;V;p];
RHS = zeros (ndof_ux+ndof_uy+ndof_p, 1);

SOL(inodes_stokes) = LHS(inodes_stokes, inodes_stokes) \ (RHS(inodes_stokes) - LHS(inodes_stokes, dnodes_stokes) * SOL(dnodes_stokes));
U = SOL (1 : ndof_ux);
V = SOL ((1 : ndof_uy) + ndof_ux);
p = SOL ((1 : ndof_p) + ndof_ux  + ndof_uy);

figure (1)
streamline (X, Y, reshape (U, size (X)), reshape (V, size (X)), [linspace(0, L, 10), linspace(L, 0, 10)], [linspace(0, H, 10), linspace(0, H, 10)]);
hold all
quiver (X, Y, reshape (U, size (X)), reshape (V, size (X)));
axis image

## Now enable convection and perform fixed point iterations
figure (2)
for kk = 1 : 100
  tic ()
  Cx = op_u0_gradu (hx, hy, nr, nc, U, V, U);
  Cy = op_u0_gradu (hx, hy, nr, nc, U, V, V);

  RHS(1:ndof_ux) = -Cx;
  RHS((1:ndof_uy)+ndof_ux) = -Cy;

  SOL(inodes_stokes) = .6*SOL(inodes_stokes) + .4* (LHS(inodes_stokes, inodes_stokes) \ (RHS(inodes_stokes) - LHS(inodes_stokes, dnodes_stokes) * SOL(dnodes_stokes)));
  U = SOL (1 : ndof_ux);
  V = SOL ((1 : ndof_uy) + ndof_ux);
  p = SOL ((1 : ndof_p) + ndof_ux  + ndof_uy);
  toc ()

  if (mod (kk, 10) == 0)
    %%close (2)
    figure (2)
    pcolor (X, Y, reshape (U, size (X)))
    colorbar ('ylabel', 'horizontal velocity')
    hold all
    hax1 = streamline (X, Y, reshape (U, size (X)), reshape (V, size (X)), ones(1, 20), linspace(0, H, 20), [ 1, 150]);
    hax2 = streamline (X, Y, reshape (U, size (X)), reshape (V, size (X)), ones(1, 20), linspace(0, H, 20), [-1, 50]);
    %%quiver (X, Y, reshape (U, size (X)), reshape (V, size (X)));
    axis image
    title (sprintf ("iteration %d", kk))
    drawnow
    hold off
  endif
  
endfor

