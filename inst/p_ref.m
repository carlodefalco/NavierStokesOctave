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
