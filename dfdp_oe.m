## Copyright (C) 2002 E. Joshua Rigler
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

## dfdp_oe.m
##
## Usage:
##
##   prt = dfdp_oe (x, f, p)
##
## x           = independent variables
## f           = predicted values
## p           = parameters; this should be a column vector,
##               with the filters organized like [B;F].  The length
##               of each filter must be designated in the global structure
##               leasqr_misc.[b|a|f|c|d]
##
##
## This function calculates the partial derivatives df/dp for 
## an OE model for use with leasqr.m.  The input parameters are
## therefore severly restricted.  This routine is more specialized
## than the default dfdp routine, which numerically calculates the
## partials, so it should provide faster convergence.
##

function prt = dfdp_oe (x, f, p)

  ## Set the leasqr_misc structure to global
  global leasqr_misc;

  ## Determine the lengths of the various polynomials
  ## (note, you MUST set the leasqr_misc structure
  ##  before calling this function)
  nb = leasqr_misc.nb;
  nf = leasqr_misc.nf;

  ## Make sure x, f, and p are column vectors
  if (columns(x) > rows(x) )
    x = x';
  endif

  if (columns(f) > rows(f) )
    f = f';
  endif

  if (columns(p) > rows(p))
    p = p';
  endif


  ## Pull B and F out of P
  B = p (1:nb);
  F = p (nb+1:length(p));


  ## Initialize Jacobian matrix
  prt = zeros (length(f), length(p) );

  ## Calculate partials
  for j=1:nb

    ## dy/db
    prt (:,j) = filter ([1], [1;F], \
    			[zeros(j-1,1); x(1:length(x)-j+1)]);
  endfor

  for j=1:nf

    ## dy/dF
    prt (:,j+nb) = -(filter ([1], [1;F], \
    			     [zeros(j,1); f(1:length(f)-j)]) );
  endfor



endfunction
