function A = ref_gradu_gradv ()
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
