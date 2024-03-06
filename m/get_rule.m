## function [nx, ny, weights] = get_rule ()
  
##   bp = [-9.061798459386640e-01
## 	   -5.384693101056831e-01
## 	    0
## 	    5.384693101056831e-01
## 	    9.061798459386640e-01];
##   bp = (bp + 1)/2;

##   wf = [ 2.369268850561891e-01
## 	    4.786286704993664e-01
## 	    5.688888888888889e-01
## 	    4.786286704993664e-01
## 	    2.369268850561891e-01];
##   wf = 1/2 * wf;

##   [nx, ny] = meshgrid (bp);
##   weights = wf * wf';

## endfunction

function [nx, ny, weights] = get_rule ()
  
  bp =  [0; .5; 1];
  wf = [1/6; 2/3; 1/6];

  [nx, ny] = meshgrid (bp);
  weights = wf * wf';

endfunction
