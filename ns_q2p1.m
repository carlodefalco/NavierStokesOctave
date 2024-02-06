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



function ii = rc2i (r, c, nrows, ncols=c)
  ii = sub2ind([nrows, ncols], r, c);
endfunction

function [r, c] = i2rc (ii, nrows, ncols = 1)
  [r, c] = ind2sub ([nrows, ncols], ii);
endfunction

function p = p_ref (ii, x, y)
  switch (ii)
    case 1
      p = 1 - x - y;
    case 2
      p = x;
    case 3
      p = y;
    otherwise
      p = nan;
  endswitch
endfunction

function v = v_ref (ii, x, y)
  [r, c] = i2rc (ii, 3);
  v = ones (size (x));
  switch (c)
    case 1
      v .*= (x-1/2).*(x-1)   ./ ((0-1/2) .* (0-1));
    case 2
      v .*= (x-0)  .*(x-1)   ./ ((1/2-0) .* (1/2-1));
    case 3
      v .*= (x-0)  .*(x-1/2) ./ ((1-0)   .* (1-1/2));
    otherwise
      v = nan;
  endswitch

  switch (r)
    case 1
      v .*= (y-1/2).*(y-1)   ./ ((0-1/2) .* (0-1));
    case 2
      v .*= (y-0)  .*(y-1)   ./ ((1/2-0) .* (1/2-1));
    case 3
      v .*= (y-0)  .*(y-1/2) ./ ((1-0)   .* (1-1/2));
    otherwise
      v = nan;
  endswitch
endfunction

function v = v_dx_ref (ii, x, y)
  [r, c] = i2rc (ii, 3);
  v = ones (size (x));
  switch (c)
    case 1
      v .*= ((x-1/2) + (x-1))   ./ ((0-1/2) .* (0-1));
    case 2
      v .*= ((x-0)   + (x-1))   ./ ((1/2-0) .* (1/2-1));
    case 3
      v .*= ((x-0)   + (x-1/2)) ./ ((1-0)   .* (1-1/2));
    otherwise
      v = nan;
  endswitch

  switch (r)
    case 1
      v .*= (y-1/2).*(y-1)   ./ ((0-1/2) .* (0-1));
    case 2
      v .*= (y-0)  .*(y-1)   ./ ((1/2-0) .* (1/2-1));
    case 3
      v .*= (y-0)  .*(y-1/2) ./ ((1-0)   .* (1-1/2));
    otherwise
      v = nan;
  endswitch
endfunction

function v = v_dy_ref (ii, x, y)
  [r, c] = i2rc (ii, 3);
  v = ones (size (x));
  switch (c)
    case 1
      v .*= ((x-1/2) .* (x-1))   ./ ((0-1/2) .* (0-1));
    case 2
      v .*= ((x-0)   .* (x-1))   ./ ((1/2-0) .* (1/2-1));
    case 3
      v .*= ((x-0)   .* (x-1/2)) ./ ((1-0)   .* (1-1/2));
    otherwise
      v = nan;
  endswitch

  switch (r)
    case 1
      v .*= ((y-1/2) + (y-1))   ./ ((0-1/2) .* (0-1));
    case 2
      v .*= ((y-0)   + (y-1))   ./ ((1/2-0) .* (1/2-1));
    case 3
      v .*= ((y-0)   + (y-1/2)) ./ ((1-0)   .* (1-1/2));
    otherwise
      v = nan;
  endswitch
endfunction

%{
[x, y] = meshgrid (linspace (0,1,100));
for ii = 1 : 9
v = v_ref (ii, x, y);
figure
mesh (x, y, v)
endfor
%}

function [nx, ny, weights] = get_rule ()
  
  bp =  [-9.061798459386640e-01
	 -5.384693101056831e-01
	 0
	 5.384693101056831e-01
	 9.061798459386640e-01];
  bp = (bp + 1)/2;

  wf = [ 2.369268850561891e-01
	 4.786286704993664e-01
	 5.688888888888889e-01
	 4.786286704993664e-01
	 2.369268850561891e-01];
  wf = 1/2 * wf;

  [nx, ny] = meshgrid (bp);
  weights = wf * wf';

endfunction

function A = ref_gradu_gradv ()
  ##  28/45       -1/5      -1/30       -1/5     -16/45        1/9      -1/30        1/9      -1/45
  ##   -1/5      88/45       -1/5     -16/45     -16/15     -16/45        1/9          0        1/9
  ##  -1/30       -1/5      28/45        1/9     -16/45       -1/5      -1/45        1/9      -1/30
  ##   -1/5     -16/45        1/9      88/45     -16/15          0       -1/5     -16/45        1/9
  ## -16/45     -16/15     -16/45     -16/15     256/45     -16/15     -16/45     -16/15     -16/45
  ##    1/9     -16/45       -1/5          0     -16/15      88/45        1/9     -16/45       -1/5
  ##  -1/30        1/9      -1/45       -1/5     -16/45        1/9      28/45       -1/5      -1/30
  ##    1/9          0        1/9     -16/45     -16/15     -16/45       -1/5      88/45       -1/5
  ##  -1/45        1/9      -1/30        1/9     -16/45       -1/5      -1/30       -1/5      28/45
  [nx, ny, weights] = get_rule ();
  A = zeros (9);
  for is = 1 : 9
    dudx = v_dx_ref (is, nx(:), ny(:));
    dudy = v_dy_ref (is, nx(:), ny(:));
    for js = is : 9
      dvdx = v_dx_ref (js, nx(:), ny(:));
      dvdy = v_dy_ref (js, nx(:), ny(:));
      A(is, js) = weights(:)' * (dudx .* dvdx + dudy .* dvdy);
      A(js, is) = A(is, js);
    endfor
  endfor
endfunction

function A = ref_u_v ()
  ##  4/225      2/225     -1/225      2/225      1/225     -1/450     -1/225     -1/450      1/900
  ##  2/225     16/225      2/225      1/225      8/225      1/225     -1/450     -4/225     -1/450
  ## -1/225      2/225      4/225     -1/450      1/225      2/225      1/900     -1/450     -1/225
  ##  2/225      1/225     -1/450     16/225      8/225     -4/225      2/225      1/225     -1/450
  ##  1/225      8/225      1/225      8/225     64/225      8/225      1/225      8/225      1/225
  ## -1/450      1/225      2/225     -4/225      8/225     16/225     -1/450      1/225      2/225
  ## -1/225     -1/450      1/900      2/225      1/225     -1/450      4/225      2/225     -1/225
  ## -1/450     -4/225     -1/450      1/225      8/225      1/225      2/225     16/225      2/225
  ##  1/900     -1/450     -1/225     -1/450      1/225      2/225     -1/225      2/225      4/225
  [nx, ny, weights] = get_rule ();
  A = zeros (9);
  for is = 1 : 9
    u = v_ref (is, nx(:), ny(:));
    for js = is : 9
      v = v_ref (js, nx(:), ny(:));
      A(is, js) = weights(:)' * (u(:) .* v(:));
      A(js, is) = A(is, js);
    endfor
  endfor
endfunction

function A = ref_dudx_p ()
  ## -5/36      -1/36          *
  ##  -2/9       -1/9       -1/3
  ##  1/36      -1/36       -1/6
  ##   1/9       -1/9          *
  ##   4/9       -4/9          *
  ##   1/9       -1/9          *
  ##  1/36       5/36          *
  ##  -2/9        5/9        1/3
  ## -5/36       5/36        1/6
  [nx, ny, weights] = get_rule ();
  A = zeros (9, 3);
  for is = 1 : 9
    u = v_dx_ref (is, nx(:), ny(:));
    for js = 1 : 3
      p = p_ref (js, nx(:), ny(:));
      A(is, js) = weights(:)' * (u(:) .* p(:));
    endfor
  endfor
endfunction

function A = ref_dudy_p ()
  ## -5/36          *      -1/36
  ##   1/9          *       -1/9
  ##  1/36          *       5/36
  ##  -2/9       -1/3       -1/9
  ##   4/9          *       -4/9
  ##  -2/9        1/3        5/9
  ##  1/36       -1/6      -1/36
  ##   1/9          *       -1/9
  ## -5/36        1/6       5/36
  [nx, ny, weights] = get_rule ();
  A = zeros (9, 3);
  for is = 1 : 9
    u = v_dy_ref (is, nx(:), ny(:));
    for js = 1 : 3
      p = p_ref (js, nx(:), ny(:));
      A(is, js) = weights(:)' * (u(:) .* p(:));
    endfor
  endfor
endfunction


function iglob = q2_connectivity (iloc, iel, numrows, numcols=1)
  [rel, cel]   = i2rc (iel, numrows, numcols);
  [rloc, cloc] = i2rc (iloc, 3, 3);
  r = 2*(rel-1) + rloc;
  c = 2*(cel-1) + cloc;
  iglob = rc2i (r, c, 2*numrows+1, 2*numcols+1);
endfunction

function iglob = p1_connectivity (iloc, iel, numrows=1, numcols=1)
  iglob = iel-1 + iloc;
endfunction

function A = op_gradu_gradv (numrows, numcols)
  Aref = ref_gradu_gradv ();
  A = sparse ((2*numrows+1)*(2*numcols+1), (2*numrows+1)*(2*numcols+1));
  for iel = 1 : numcols*numrows
    for jnode = 1 : 9
      jglob = q2_connectivity (jnode, iel, numrows, numcols);
      for inode = 1 : 9
	iglob = q2_connectivity (inode, iel, numrows, numcols);
	A(iglob, jglob) += Aref(inode, jnode);
      endfor
    endfor
  endfor
endfunction

function A = op_u_v (hx, hy, numrows, numcols)
  Aref = ref_u_v ();
  A = sparse ((2*numrows+1)*(2*numcols+1), (2*numrows+1)*(2*numcols+1));
  for iel = 1 : numcols*numrows
    for jnode = 1 : 9
      jglob = q2_connectivity (jnode, iel, numrows, numcols);
      for inode = 1 : 9
	iglob = q2_connectivity (inode, iel, numrows, numcols);
	A(iglob, jglob) += (hx*hy) * Aref(inode, jnode);
      endfor
    endfor
  endfor
endfunction

function A = op_dudx_p (hx, hy, numrows, numcols)
  Aref = ref_dudx_p ();
  A = sparse ((2*numrows+1)*(2*numcols+1), (3*numrows*numcols));
  for iel = 1 : numcols*numrows
    for inode = 1 : 9
      iglob = q2_connectivity (inode, iel, numrows, numcols);
      for jnode = 1 : 3
	jglob = p1_connectivity (jnode, iel, numrows, numcols);
	A(iglob, jglob) += hy * Aref(inode, jnode);
      endfor
    endfor
  endfor
endfunction

function A = op_dudy_p (hx, hy, numrows, numcols)
  Aref = ref_dudy_p ();
  A = sparse ((2*numrows+1)*(2*numcols+1), (3*numrows*numcols));
  for iel = 1 : numcols*numrows
    for inode = 1 : 9
      iglob = q2_connectivity (inode, iel, numrows, numcols);
      for jnode = 1 : 3
	jglob = p1_connectivity (jnode, iel, numrows, numcols);
	A(iglob, jglob) += hx * Aref(inode, jnode);
      endfor
    endfor
  endfor
endfunction

function C = op_u0_gradu (hx, hy, numrows, numcols, ux0, uy0, u)
  [nx, ny, weights] = get_rule ();
  C = zeros ((2*numrows+1)*(2*numcols+1), 1);
  for iel = 1 : numcols*numrows
    for inode = 1 : 9
      iglob = q2_connectivity (inode, iel, numrows, numcols);
      tmp = v_ref (inode, nx(:), ny(:));
      x = ux0(iglob) * tmp;
      y = uy0(iglob) * tmp;
      for jnode = 1 : 9
	jglob = q2_connectivity (jnode, iel, numrows, numcols);
	dx = u(jglob) * v_dx_ref (jnode, nx(:), ny(:));
	dy = u(jglob) * v_dy_ref (jnode, nx(:), ny(:));
	C(iglob) += hy * weights(:)' * (dx(:) .* x(:)) + hx * weights(:)' * (dy(:) .* y(:));
      endfor
    endfor
  endfor
endfunction


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

%{
U = zeros (ndof_ux, 1);
U(Y(:) >= (1)*(1-eps)) = 1;
F = ones (ndof_ux, 1);
inodes = setdiff (1:ndof_ux, dnodes);

U(inodes) = A(inodes, inodes) \ ((M*F)(inodes) -A(inodes, dnodes) * U(dnodes));

figure (1)
surf (X, Y, reshape (U, size (X)))
%}


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
