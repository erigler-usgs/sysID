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

## filt_bj.m
##
## Usage:
##
##   y = filt_bj ([in_out_ts], params)
##
## in_out_ts   = input/output time series; this should be a matrix with
##               two columns.  The first column is the input time
##               time series, the second is the output time series.
##               Each row corresponds to an observation.
##
## params      = paramters to solve for; this should be a column vector,
##               with the filters organized like [B;F;C;D].  The length
##               of each filter must be designated in the global structure
##               leasqr_misc.[b|a|f|c|d]
##
##
## This simple function is primarily designed to work with the
## "leasqr.m" non-linear least-squares optimization function.
## Therefore the input parameters are severely restricted, but 
## othewise it should work nicely as a simulator for Box-Jenkins 
## models.
## 
##
## Still to do:
##
## - Make it handle multiple inputs/outputs (difficult or impossible)
## - Allow B to be of arbitrary length (difficult)
##

function y = filt_bj (iots, p)

  ## Set the leasqr_misc structure to global
  global leasqr_misc;

  ## Determine the lengths of the various polynomials
  ## (note, you MUST set the leasqr_misc structure
  ##  before calling this function)
  nb = leasqr_misc.nb;
  nf = leasqr_misc.nf;
  nc = leasqr_misc.nc;
  nd = leasqr_misc.nd;

  ## Make sure that the filter is a vector of scalars
  if (min(size(p)) > 1 )
    error (["The coefficients must be in a single vector of scalars\n"]);
  endif

  ## Make sure that the iots is a two column matrix
  if (min(size(iots)) > 2 )
    error (["The input/output time series must be a two-column matrix\n",
	    "with rows corresponding to different observation times\n"]);
  endif


  ## Make sure that the filter and input vectors are
  ## both made of column vectors
  if (columns(p) > rows(p) )
    p = p';
  endif

  ## This is probably safe, since a time series by definition
  ## must be more than a single observation, and will very likely
  ## be much more than two observations
  if (columns(iots) > rows(iots))
    iots=iots';
  endif


  ## Pull A, B and C out of P
  B = p (1:nb);
  F = p (nb+1:nf+nb);
  C = p (nb+nf+1:nb+nf+nc);
  D = p (nb+nf+nc+1:length(p));


  ## Pull ITS and OTS out of IOTS
  its = iots (:,1);
  ots = iots (:,2);


  ## for now we won't worry about non-zero initial conditions
  y = filter (B, [1;F], filter ([1;D], [1;C], its) ) + \
      filter ([[1;C]-[1;D]], [1;C], ots);

endfunction
