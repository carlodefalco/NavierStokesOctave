function iglob = q2_connectivity (iloc, iel, numrows, numcols=1)
  [rel, cel]   = i2rc (iel, numrows, numcols);
  [rloc, cloc] = i2rc (iloc, 3, 3);
  r = 2*(rel-1) + rloc;
  c = 2*(cel-1) + cloc;
  iglob = rc2i (r, c, 2*numrows+1, 2*numcols+1);
endfunction
