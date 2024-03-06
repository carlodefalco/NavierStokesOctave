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
