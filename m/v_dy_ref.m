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
