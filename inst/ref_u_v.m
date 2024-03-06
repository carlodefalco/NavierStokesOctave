function A = ref_u_v ()
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
