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

## dfdp_arx.m
##
## Usage:
##
##   prt = dfdp_arx (x, f, p)
##
## x           = input/output matrix, with column 1 being the input, and
##               column 2 being the output
## f           = predicted values
## p           = parameters; this should be a column vector,
##               with the filters organized like [B;A].  The length
##               of each filter must be designated in the global structure
##               leasqr_misc.[b|a|f|c|d]
##
##
## This function calculates the partial derivatives df/dp for 
## an ARX model for use with leasqr.m.  The input parameters are
## therefore severly restricted.  This routine is more specialized
## than the default dfdp routine, which numerically calculates the
## partials, so it _should_ provide faster convergence.
##

function prt = dfdp_arx (x, f, p)

  ## Set the leasqr_misc structure to global
  global leasqr_misc;

  ## Determine the lengths of the various polynomials
  ## (note, you MUST set the leasqr_misc structure
  ##  before calling this function)
  nb = leasqr_misc.nb;
  na = leasqr_misc.na;

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


  ## How many inputs?
  nits = columns (x) - 1;


  ## Pull A and B out of P
  for i=1:nits
    B(:,i) = p ( (nb*(i-1))+1:nb*i);
  endfor
  A = p (nb*i+1:length(p));


  ## Pull ITS and OTS out of x
  its = x (:,1:nits);
  ots = x (:,nits+1);


  ## Initialize Jacobian matrix
  prt = zeros (length(f), length(p) );

  ## Calculate partials
  for i=1:nits
    for j=1:nb
      
      ## dy/db
      prt (:,j + ((i-1)*nb) ) = filter ([1], [1], \
    					[zeros(j-1,1); its(1:length(its)-j+1,i)]);
    endfor
  endfor

  for j=1:na
    ## dy/dA
    prt (:,j + nb*i) = -filter ([1], [1], \
    				  [zeros(j,1); ots(1:length(ots)-j)]);
    
  endfor



endfunction
