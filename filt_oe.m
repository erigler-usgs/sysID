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

## filt_oe.m
##
## Usage:
##
##   y = filt_oe (in_ts, params)
##
## in_ts       = input time series, or independent variables
## params      = paramters to solve for; this should be a column vector,
##               with the filters organized like [B;F].  The length
##               of each filter must be designated in the global structure
##               leasqr_misc.[b|a|f|c|d]
##
##
## This simple function is primarily designed to work with the 
## "leasqr.m" non-linear least-squares optimization function.
## Therefore the input parameters are severely restricted, but 
## othewise it should work nicely as a simulator for OE models.
## 
## Still to do:
##
## - Make it handle multiple inputs/outputs (difficult or impossible)
## - Allow B to be of arbitrary length (difficult)
##

function y = filt_oe (its, p)

  ## Set the leasqr_misc structure to global
  global leasqr_misc;

  ## Determine the lengths of the various polynomials
  ## (note, you MUST set the leasqr_misc structure
  ##  before calling this function)
  nb = leasqr_misc.nb;
  nf = leasqr_misc.nf;

  ## Make sure that the filter and input vectors are
  ## both vectors of scalars
  if (min(size(its)) > 1 ||
      min(size(p))   > 1 )
    error (["Can only handle vectors of scalars for now\n"]);
  endif

  ## Make sure that the filter and input vectors are
  ## both column vectors
  if (columns(p) > rows(p) )
    p = p';
  endif

  if (columns(its) > rows(its))
    its=its';
  endif


  ## Pull F and B out of P
  B = p (1:nb);
  F = p (nb+1:length(p));

  ## for now we won't worry about non-zero initial conditions
  y = filter ( [B], [1;F], its);


endfunction
