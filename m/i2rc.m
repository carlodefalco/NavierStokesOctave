function [r, c] = i2rc (ii, nrows, ncols = 1)
  [r, c] = ind2sub ([nrows, ncols], ii);
endfunction
