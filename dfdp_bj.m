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

## dfdp_bj.m
##
## Usage:
##
##   prt = dfdp_bj (x, f, p)
##
## x           = input and output time series; should be a matrix comprised of
##               two equal length column vectors, with column 1 being the input
##               and column 2 being the output time series
## f           = predicted values
## p           = parameters; this should be a column vector,
##               with the filters organized like [B;F;C;D].  The length
##               of each filter must be designated in the global structure
##               leasqr_misc.[b|a|f|c|d]
##
##
## This function calculates the partial derivatives df/dp for 
## an Box-Jenkins model for use with leasqr.m.  The input parameters are
## therefore severly restricted.  This routine is more specialized
## than the default dfdp routine, which numerically calculates the
## partials, so it should provide faster convergence.
##

function prt = dfdp_bj (x, f, p)

  ## Set the leasqr_misc structure to global
  global leasqr_misc;

  ## Determine the lengths of the various polynomials
  ## (note, you MUST set the leasqr_misc structure
  ##  before calling this function)
  nb = leasqr_misc.nb;
  nf = leasqr_misc.nf;
  nc = leasqr_misc.nc;
  nd = leasqr_misc.nd;

  ## Make sure x, f, and p are column vectors (for now, we assume
  ## that the time series is >> two steps, and use this following
  ## simple heuristic to arrange 'x' appropriately
  if (columns(x) > rows(x) )
    x = x';
  endif

  if (columns(f) > rows(f) )
    f = f';
  endif

  if (columns(p) > rows(p))
    p = p';
  endif


  ## Pull A, B and C out of P
  B = p (1:nb);
  F = p (nb+1:nf+nb);
  C = p (nb+nf+1:nb+nf+nc);
  D = p (nb+nf+nc+1:length(p));


  ## Pull ITS and OTS out of x
  its = x (:,1);
  ots = x (:,2);


  ## Initialize Jacobian matrix
  prt = zeros (length(f), length(p) );

  ## Calculate partials
  resids = ots-f;
  for j=1:nb
    ## dy/db
    prt (:,j) = \
	filter ([1], [1;F], \
		filter ( [1;D],[1;C], \
			[zeros(j-1,1); its(1:length(its)-j+1)] ) );
  endfor

  for j=1:nf
    ## dy/df
    prt (:,j+nb) = \
	-filter ([1], [1;F], \
		 filter ([1;D], [1;C], \
			 [zeros(j,1); ots(1:length(ots)-j)]) + \
		 [zeros(j,1); f(1:length(f)-j)]);
  endfor

  for j=1:nc
    ## dy/dc
    prt (:,j+nb+nf) = \
	filter ([1], [1;C], [zeros(j,1); resids(1:length(resids)-j)]);
  endfor

  for j=1:nd
    ## dy/dd
    prt (:,j+nb+nf+nc) = \
	filter ([1], [1;F], \
		filter([B], [1;C], \
		       [zeros(j,1); its(1:length(its)-j)]) ) - \
	filter ([1], [1;C], [zeros(j,1); ots(1:length(ots)-j)]);
  endfor

endfunction
