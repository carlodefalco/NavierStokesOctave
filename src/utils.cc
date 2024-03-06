/*

Copyright (C) 2024 Carlo de Falco

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

Author: Carlo de Falco <carlo@guglielmo.lan>
Created: 2024-02-06

*/

#include <octave/oct.h>
#include <iostream>

constexpr double nx[25]{
  4.691007703066802e-02, 4.691007703066802e-02, 4.691007703066802e-02, 4.691007703066802e-02, 4.691007703066802e-02,
    2.307653449471584e-01, 2.307653449471584e-01, 2.307653449471584e-01, 2.307653449471584e-01, 2.307653449471584e-01,
    5.000000000000000e-01, 5.000000000000000e-01, 5.000000000000000e-01, 5.000000000000000e-01, 5.000000000000000e-01,
    7.692346550528415e-01, 7.692346550528415e-01, 7.692346550528415e-01, 7.692346550528415e-01, 7.692346550528415e-01,
    9.530899229693319e-01, 9.530899229693319e-01, 9.530899229693319e-01, 9.530899229693319e-01, 9.530899229693319e-01};

constexpr double ny[25]{
  4.691007703066802e-02, 2.307653449471584e-01, 5.000000000000000e-01, 7.692346550528415e-01, 9.530899229693319e-01,
    4.691007703066802e-02, 2.307653449471584e-01, 5.000000000000000e-01, 7.692346550528415e-01, 9.530899229693319e-01,
    4.691007703066802e-02, 2.307653449471584e-01, 5.000000000000000e-01, 7.692346550528415e-01, 9.530899229693319e-01,
    4.691007703066802e-02, 2.307653449471584e-01, 5.000000000000000e-01, 7.692346550528415e-01, 9.530899229693319e-01,
    4.691007703066802e-02, 2.307653449471584e-01, 5.000000000000000e-01, 7.692346550528415e-01, 9.530899229693319e-01};

constexpr double weights[25]{
  1.403358721560716e-02, 2.835000000000000e-02, 3.369626809688023e-02, 2.835000000000000e-02, 1.403358721560716e-02,
    2.835000000000000e-02, 5.727135105599777e-02, 6.807163313768767e-02, 5.727135105599777e-02, 2.835000000000000e-02,
    3.369626809688023e-02, 6.807163313768767e-02, 8.090864197530864e-02, 6.807163313768767e-02, 3.369626809688023e-02,
    2.835000000000000e-02, 5.727135105599777e-02, 6.807163313768767e-02, 5.727135105599777e-02, 2.835000000000000e-02,
    1.403358721560716e-02, 2.835000000000000e-02, 3.369626809688023e-02, 2.835000000000000e-02, 1.403358721560716e-02};

Array<double>
p_ref (octave_idx_type ii, const Array<double> x, const Array<double> y)
{
  Array<double> ret(x.dims());
  switch (ii) {
  case 1 :
    for (auto kk = 0; kk < x.numel (); ++kk)
      ret.xelem (kk) = 1. - x(kk) - y(kk);
    break;

  case 2 :
    for (auto kk = 0; kk < x.numel (); ++kk)
      ret.xelem (kk) = x(kk);
    break;
  case 3 :
    for (auto kk = 0; kk < x.numel (); ++kk)
      ret.xelem (kk) = y(kk);
    break;
  default :
    throw ("p_ref: incorrect basis function number\n");
    break;
  }
  return ret;
};

// PKG_ADD: autoload ('p_ref', canonicalize_file_name ('utils.oct'))
DEFUN_DLD(p_ref, args, nargout,
	  "-*- texinfo -*-\n\
@deftypefn {} {@var{retval} =} utils (@var{input1}, @var{input2})\n\
@seealso{}\n\
@end deftypefn")
{
  octave_value_list retval;
  int nargin = args.length ();

  retval(0) = p_ref (args(0).idx_type_value (), args(1).array_value (),
		     args(2).array_value ());

  return retval;
}

Array<double>
v_ref (octave_idx_type ii, const Array<double> x, const Array<double> y) {
  Array<double> v(x.dims());
  Array<idx_vector> rc = ind2sub ({3, 0}, ii-1);
  int row = rc(0)(0)+1;
  int col = rc(1)(0)+1;
  switch (col) {
  case 1 :
    for (auto kk = 0; kk < x.numel (); ++ kk) {
      v (kk) = (x(kk)-1./2.) * (x(kk)-1.) / ((0.-1./2.) * (0.-1.));
    }
    break;
  case 2 :
    for (auto kk = 0; kk < x.numel (); ++ kk) {
      v (kk) = x(kk) * (x(kk)-1.) / ((1./2.-0.) * (1./2-1));
    }
    break;
  case 3 :
    for (auto kk = 0; kk < x.numel (); ++ kk) {
      v (kk) = x(kk) * (x(kk)-1./2.) / ((1.-0) * (1.-1./2.));
    }
    break;
  default :
    throw ("v_ref: incorrect basis function number\n");
    break;
  }

  switch (row) {
  case 1 :
    for (auto kk = 0; kk < y.numel (); ++ kk) {
      v (kk) *= (y(kk)-1/2) * (y(kk)-1) / ((0.-1./2.) * (0.-1.));
    }
    break;
  case 2 :
    for (auto kk = 0; kk < x.numel (); ++ kk) {
      v (kk) *= y(kk) * (y(kk)-1.) / ((1./2.-0.) * (1./2.-1.));
    }
    break;
  case 3 :
    for (auto kk = 0; kk < y.numel (); ++ kk) {
      v (kk) *= y(kk) * (y(kk)-1./2.) / ((1.-0) * (1.-1./2.));
    }
    break;
  default :
    throw ("v_ref: incorrect basis function number\n");
    break;
  }
  return v;
};

// PKG_ADD: autoload ('v_ref', canonicalize_file_name ('utils.oct'))
DEFUN_DLD(v_ref, args, nargout,
	  "-*- texinfo -*-\n\
@deftypefn {} {@var{retval} =} utils (@var{input1}, @var{input2})\n\
@seealso{}\n\
@end deftypefn")
{
  octave_value_list retval;
  int nargin = args.length ();

  retval(0) = v_ref (args(0).idx_type_value (), args(1).array_value (),
		     args(2).array_value ());

  return retval;
}

Array<double>
v_dx_ref (octave_idx_type ii, const Array<double> x, const Array<double> y) {
  Array<double> v(x.dims());
  Array<idx_vector> rc = ind2sub ({3, 0}, ii-1);
  int row = rc(0)(0)+1;
  int col = rc(1)(0)+1;
  switch (col) {
  case 1 :
    for (auto kk = 0; kk < x.numel (); ++ kk) {
      v(kk) = ((x(kk)-1./2.) + (x(kk)-1.)) / ((0.-1./2.)*(0.-1.));
    }
    break;
  case 2 :
    for (auto kk = 0; kk < x.numel (); ++ kk) {
      v(kk) = ((x(kk)-0.) + (x(kk)-1.)) / ((1./2.-0.) * (1./2.-1.));
    }
    break;
  case 3 :
    for (auto kk = 0; kk < x.numel (); ++ kk) {
      v(kk) = ((x(kk)-0.) + (x(kk)-1./2.)) / ((1.-0.) * (1.-1./2.));
    }
    break;
  default :
    throw ("v_dx_ref: incorrect basis function number\n");
  }

  switch (row) {
  case 1 :
    for (auto kk = 0; kk < x.numel (); ++ kk) {
      v(kk) *= (y(kk)-1./2.) * (y(kk)-1.) / ((0.-1./2.) * (0.-1.));
    }
    break;
  case 2 :
    for (auto kk = 0; kk < x.numel (); ++ kk) {
      v(kk) *= (y(kk)-0.) *(y(kk)-1) / ((1./2.-0) * (1./2.-1.));
    }
    break;
  case 3 :
    for (auto kk = 0; kk < x.numel (); ++ kk) {
      v(kk) *= (y(kk)-0.)*(y(kk)-1./2.) / ((1.-0.) * (1.-1./2.));
    }
    break;
  default :
    throw ("v_dx_ref: incorrect basis function number\n");
  }
  return v;
}


// PKG_ADD: autoload ('v_dx_ref', canonicalize_file_name ('utils.oct'))
DEFUN_DLD(v_dx_ref, args, nargout,
	  "-*- texinfo -*-\n\
@deftypefn {} {@var{retval} =} utils (@var{input1}, @var{input2})\n\
@seealso{}\n\
@end deftypefn")
{
  octave_value_list retval;
  int nargin = args.length ();

  retval(0) = v_dx_ref (args(0).idx_type_value (), args(1).array_value (),
			args(2).array_value ());

  return retval;
}

Array<double>
v_dy_ref (octave_idx_type ii, const Array<double> x, const Array<double> y) {
  Array<double> v(x.dims());
  Array<idx_vector> rc = ind2sub ({3, 0}, ii-1);
  int row = rc(0)(0)+1;
  int col = rc(1)(0)+1;
  switch (col) {
  case 1 :
    for (auto kk = 0; kk < x.numel (); ++ kk) {
      v(kk) = ((x(kk)-1./2.) * (x(kk)-1.)) / ((0.-1./2.) * (0.-1.));
    }
    break;
  case 2 :
    for (auto kk = 0; kk < x.numel (); ++ kk) {
      v(kk) = ((x(kk)-0.) * (x(kk)-1.)) / ((1./2.-0.) * (1./2.-1.));
    }
    break;
  case 3 :
    for (auto kk = 0; kk < x.numel (); ++ kk) {
      v(kk) = ((x(kk)-0.) * (x(kk)-1./2.)) / ((1.-0.) * (1.-1./2.));
    }
    break;
  default :
    throw ("v_dy_ref: incorrect basis function number\n");
  }

  switch (row) {
  case 1 :
    for (auto kk = 0; kk < x.numel (); ++ kk) {
      v(kk) *= ((y(kk)-1./2.) + (y(kk)-1.)) / ((0.-1./2.) * (0.-1.));
    }
    break;
  case 2 :
    for (auto kk = 0; kk < x.numel (); ++ kk) {
      v(kk) *= ((y(kk)-0.) + (y(kk)-1.)) / ((1./2.-0.) * (1./2.-1.));
    }
    break;
  case 3 :
    for (auto kk = 0; kk < x.numel (); ++ kk) {
      v(kk) *= ((y(kk)-0.) + (y(kk)-1./2.)) / ((1.-0.) * (1.-1./2.));
    }
    break;
  default :
    throw ("v_dy_ref: incorrect basis function number\n");
  }
  return v;
}

// PKG_ADD: autoload ('v_dy_ref', canonicalize_file_name ('utils.oct'))
DEFUN_DLD(v_dy_ref, args, nargout,
	  "-*- texinfo -*-\n\
// @deftypefn {} {@var{retval} =} utils (@var{input1}, @var{input2})\n\
// @seealso{}\n\
// @end deftypefn")
{
  octave_value_list retval;
  int nargin = args.length ();

  retval(0) = v_dy_ref (args(0).idx_type_value (), args(1).array_value (),
			args(2).array_value ());

  return retval;
}

// PKG_ADD: autoload ('op_u0_gradu', canonicalize_file_name ('utils.oct'))
DEFUN_DLD(op_u0_gradu, args, nargout,
	  "-*- texinfo -*-\n\
// @deftypefn {} {@var{retval} =} utils (@var{input1}, @var{input2})\n\
// @seealso{}\n\
// @end deftypefn")
{

  double hx = args(0).double_value ();
  double hy = args(1).double_value ();
  octave_idx_type numrows = args(2).idx_type_value ();
  octave_idx_type numcols = args(3).idx_type_value ();
  ColumnVector ux0 = args(4).column_vector_value ();
  ColumnVector uy0 = args(5).column_vector_value ();
  ColumnVector u = args(6).column_vector_value ();

  ColumnVector C((2*numrows+1)*(2*numcols+1), 0.);

  
  for iel = 1 : numcols*numrows
    for inode = 1 : 9
      iglob = q2_connectivity (inode, iel, numrows, numcols);
      tmp = v_ref (inode, nx(:), ny(:));
      x = ux0(iglob) * tmp;
      y = uy0(iglob) * tmp;
      for jnode = 1 : 9
	jglob = q2_connectivity (jnode, iel, numrows, numcols);
	dx = u(jglob) * v_dx_ref (jnode, nx(:), ny(:));
	dy = u(jglob) * v_dy_ref (jnode, nx(:), ny(:));
	C(iglob) += hy * weights(:)' * (dx(:) .* x(:)) + hx * weights(:)' * (dy(:) .* y(:));
      endfor
    endfor
  endfor
	}
