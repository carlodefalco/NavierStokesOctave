function A = ref_p_q ()
  [nx, ny, weights] = get_rule ();
  A = zeros (3);
  for is = 1 : 3
    p = p_ref (is, nx(:), ny(:));
    for js = is : 3
      q = q_ref (js, nx(:), ny(:));
      A(is, js) = weights(:)' * (p(:) .* q(:));
      A(js, is) = A(is, js);
    endfor
  endfor
endfunction
