## ## Transient solver, rectangular domain with cylinder
##
## Solve 
##
## $$ \partial_{t} u_x - \nu\ (\partial_{xx} u_x + \partial_{yy} u_x) + u_x\ \partial_x u_x + u_y\ \partial_y u_x + \partial_x\ p = f_x $$
##
## $$ \partial_{t} u_y - \nu\ (\partial_{xx} u_y + \partial_{yy} u_y) + u_x\ \partial_{x} u_y + u_y\ \partial_{y} u_y + \partial_{y}\ p  = f_y $$
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
## $$  \partial_{t} u_x - \nu\ (\partial_{xx} u_x^{(k+1)} + \partial_{yy} u_x^{(k+1)}) + u_x^{(k)}\ \partial_{x} u_x^{(k)} + u_y^{(k)}\ \partial_{y} u_x^{(k)} + \partial_{x}\ p^{(k+1)} = f_x $$
##
## $$  \partial_{t} u_y - \nu\ (\partial_{xx} u_y^{(k+1)} + \partial_{yy} u_y^{(k+1)}) + u_x^{(k)}\ \partial_{x} u_y^{(k)} + u_y^{(k)}\ \partial_{y} u_y^{(k)} + \partial_{y}\ p^{(k+1)}  = f_y $$
##
## $$ \partial_{x} u_x^{(k+1)} + \partial_{y} u_y^{(k+1)} = 0. $$



## Example problem, flow past cylinder

nr = 32;
nc = 32;

L =  32;
H =  32;
hx = L / nc;
hy = H / nr;


A  = op_gradu_gradv (nr, nc);
M  = op_u_v (hx, hy, nr, nc);
Bx = op_dudx_p (hx, hy, nr, nc);
By = op_dudy_p (hx, hy, nr, nc);
 
[X, Y] = meshgrid (linspace (0, L, 2*nc+1), linspace (0, H, 2*nr+1));
ndof_ux = ndof_uy = numel (X);


## FIXME : periodic boundary conditions are currently broken ...
mnodesx = find ((X(:) <= eps));
snodesx = find (abs (X(:)-L) <= eps(L));

[~, jj] = sort(Y(:)(mnodesx));
mnodesx = mnodesx(jj);
[~, jj] = sort(Y(:)(snodesx));
snodesx = snodesx(jj);

mnodesy = find ((Y(:) <= eps));
snodesy = find (abs (Y(:)-H) <= eps(H));

[~, jj] = sort(X(:)(mnodesy));
mnodesy = mnodesy(jj);
[~, jj] = sort(X(:)(snodesy));
snodesy = snodesy(jj);

mnodesu = [mnodesx(:); mnodesy(:)];
mnodesv = mnodesu + ndof_ux;

snodesu = [snodesx(:); snodesy];
snodesv = snodesu + ndof_uy;

unodesu = setdiff (1:ndof_ux, snodesu);
unodesv = setdiff ((1:ndof_uy)+ndof_ux, snodesv);

unodes = union (unodesu, unodesv);

ndof_p = nr*nc*3;
nu     = 1;


U = .1 * cos (2*pi* .5 * Y/H)(:);
V = 0  * cos (2*pi* 2 * X/L)(:);
p = zeros (ndof_p,  1);

dt  = .01 * L/nc/norm (U(:), inf)
LHS = [(1/dt)*M+nu*A                      sparse(rows(A), columns(A))      -Bx;
       sparse(rows(A), columns(A))        (1/dt)*M+nu*A                    -By;
       -Bx.'                              -By.'                             sparse(columns(Bx), columns(Bx))];

LHS(mnodesu, :) += LHS(snodesu, :);
LHS(:, mnodesu) += LHS(:, snodesu);

LHS(mnodesv, :) += LHS(snodesv, :);
LHS(:, mnodesv) += LHS(:, snodesv);

RHS = zeros (ndof_ux+ndof_uy+ndof_p, 1);
SOL = zeros (ndof_ux+ndof_uy+ndof_p, 1);

##figure (2)
ii = 0;
for kk = 1 : ceil (L/dt)
  tic ()
  Cx = op_u0_gradu (hx, hy, nr, nc, U, V, U);
  Cy = op_u0_gradu (hx, hy, nr, nc, U, V, V);

  RHS(1:ndof_ux) = -Cx+(1/dt)*M*U;
  RHS((1:ndof_uy)+ndof_ux) = -Cy+(1/dt)*M*V;
  RHS(mnodesu) += RHS(snodesu);
  RHS(mnodesv) += RHS(snodesv);
  
  SOL(unodes)  = LHS(unodes, unodes) \ RHS(unodes);
  SOL(snodesu) = SOL(mnodesu);
  SOL(snodesv) = SOL(mnodesv);
  
  U = SOL (1 : ndof_ux);
  V = SOL ((1 : ndof_uy) + ndof_ux);
  p = SOL ((1 : ndof_p) + ndof_ux  + ndof_uy);
  toc ()

  if (mod (kk, 2) == 0)
    hold off
    figure (2)
    pcolor (X, Y, reshape (sqrt (U.^2 + V.^2), size (X)))
    colorbar ('ylabel', 'velocity magnitude')
    hold all
    hax1 = streamline (X, Y, reshape (U, size (X)), reshape (V, size (X)), X(1:5:end,1:5:end)(:), Y(1:5:end,1:5:end)(:), [  1, 50]);
    hax2 = streamline (X, Y, reshape (U, size (X)), reshape (V, size (X)), X(1:5:end,1:5:end)(:), Y(1:5:end,1:5:end)(:), [ -1, 50]);
    %%quiver (X, Y, reshape (U, size (X)), reshape (V, size (X)));
    axis image
    title (sprintf ("iteration %d", kk))
    drawnow
    print ('-dpng', sprintf ('frame_%4.4d.png', ++ii))
    hold off
  endif
  save ('-binary', '-z', sprintf("iter_%4.4d.octbin.gz", kk), 'U', 'V', 'p', 'X', 'Y', 'kk')
endfor

