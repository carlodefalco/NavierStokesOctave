function A = ref_dudx_p ()
  [nx, ny, weights] = get_rule ();
  A = zeros (9, 3);
  for is = 1 : 9
    u = v_dx_ref (is, nx(:), ny(:));
    for js = 1 : 3
      p = p_ref (js, nx(:), ny(:));
      A(is, js) = weights(:)' * (u(:) .* p(:));
    endfor
  endfor
endfunction
