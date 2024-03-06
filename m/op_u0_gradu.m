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
