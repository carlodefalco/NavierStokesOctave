function ii = rc2i (r, c, nrows, ncols=c)
  ii = sub2ind([nrows, ncols], r, c);
endfunction
