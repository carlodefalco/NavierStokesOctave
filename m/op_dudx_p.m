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
