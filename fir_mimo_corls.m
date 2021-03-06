### Copyright (C) 2005 E. Joshua Rigler
###
### This program is free software; you can redistribute it and/or modify
### it under the terms of the GNU General Public License as published by
### the Free Software Foundation; either version 2 of the License, or
### (at your option) any later version.
###
### This program is distributed in the hope that it will be useful,
### but WITHOUT ANY WARRANTY; without even the implied warranty of
### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
### GNU General Public License for more details.
###
### You should have received a copy of the GNU General Public License
### along with this program; if not, write to the Free Software
### Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

###
### Usage:
###
### [theta, lags, lhs, rhs] = fir_mimo_corls (its, ots[, lagsize])
###
### Produces a multi-input/multi-output Finite Impulse Response (FIR) linear 
### filter based on given input and output time series.  This filter is the 
### solution to system of linear equations derived from the auto- and cross-
### correlation vectors, rather than directly from an observation-based 
### regression matrix.  This, in theory, provides a "consistent" estimate
### of the estimated parameters (i.e., given infinite training data, the
###  solution is guaranteed to converge to an optimal solution), and helps
### attenuate noise that might otherwise bias parameter estimates.
###
### INPUTS:
###
### its       - Input time series matrix.  Columns must represent different 
###             input veriables, rows must represent observations at equally
###             spaced and monotonically increasing moments in time.
###
### ots       - Output time series matrix.  Columns must represent different 
###             output veriables, rows must represent observations at equally
###             spaced and monotonically increasing moments in time.
###
### lagsize   - Either a scalar, or two element vector.  If a scalar, this
###             indicates the mimimum and maximum lag about zero lag.  If a
###             two element vector, lagsize(1) indicates the minimum lag
###             from zero, and lagsize(2) indicates the maximium lag from zero.
###
###
### OUTPUTS:
###
### theta     - MIMO fir filter matrix; rows correspond to time-lags
###             divided by the number of different inputs; columns 
###             correspond to the different outputs; (for example,
###             a 2 input (a,b), 3 output (x,y,z), 3 lag filter 
###             with lags 0, 1, and 2, the matrix would look like:
###
###                [ a(0,x) a(0,y) a(0,z) ]
###                [ b(0,x) b(0,y) b(0,z) ]
###                [ a(1,x) a(1,y) a(1,z) ]
###                [ b(1,x) b(1,y) b(1,z) ]
###                [ a(2,x) a(2,y) a(2,z) ]
###                [ b(2,x) b(2,y) b(2,z) ]
###
###             This matrix can be used in conjunction with the function 
###             filt_mimo.m and the appropriate set of imput column vectors
###             to generate predicted MIMO output.
###
### lags      - vector of corresponding lag values (units are whatever the
###             sampling interval is).  This may also be used in conjunction
###             witht he filt_mimo.m function.
###
### rhs       - full cross-correlation matrix for output/input time series.
###             It is the right hand side of the regression equation.
###
### lhs       - full auto-correlation matrix for the input time series.
###             It is the regression matrix, and left hand side of our
###             equation.  This, multiplied by theta, should give the
###             rhs matrix.  If not, something is very wrong!
###

### Revision History:
###
### 2005-07-20  First version with GNU/GPL license statement included
###             for public distribution (and possible modifications).  
###             No known bugs at this time.
###
### $Log$

function [theta,lag,PHI,CC] = fir_mimo_corls(its, ots, lagsize)

  if ( nargin < 2 || nargin > 3)
    usage (["\n\n[filter, lag] = ", \
	    "fir_mimo_corls (input_ts, output_ts [, lagsize])\n"]);
  endif

  ## Set a default lagsize of +/-10
  if nargin == 2
    lagsize = 10;
  endif

  ## If "lagsize" is a two-element vector, then assigne the max value of
  ## that vector to lagsize, and determine the min and max lag for use
  ## in calculating the correlation functions
  if (is_vector(lagsize) && max(size(lagsize)) == 2)
    min_lag = min(lagsize);
    max_lag = max(lagsize);
    maxlagsize = max(abs(lagsize));
    ## if for some reason a vector was passed of identical values, assume that
    ## the user wants one negative and one positive
    ## -- This prohibits generating 1st order models --
    ##if (min_lag == max_lag)
    ##  min_lag = -min_lag;
    ##endif
  elseif (is_scalar(lagsize) ) ## do nothing
    min_lag = -lagsize;
    max_lag = lagsize;
    maxlagsize = lagsize;
  else
    error(["\n\"lagsize\" must be either a scalar, or a two-element vector\n",\
	   "indicating the min and max lag (relative to zero) for the correlation\n",\
	   "functions."]);
  endif

  ## Make sure the input and output time series matrices are composed of
  ## column vectors.  This is a bit of a kludge, since it assumes that the
  ## length(theta)<<length(its).
  if (rows(its) < columns(its))
    its = its';
  endif
  if (rows(ots) < columns(ots))
    ots = ots';
  endif


  ## Calculate auto-correlation functions for input series
  autoin = flipud(xcorr (its, max_lag-min_lag));


  ## Arrange the auto-correlation regression matrix (Left-Hand side)
  dim_phis = [columns(its),columns(its)];
  PHI = zeros((max_lag-min_lag+1) * dim_phis(1));

  for i=1:dim_phis(1):(max_lag-min_lag+1)*dim_phis(1)
    for j=1:dim_phis(2):(max_lag-min_lag+1)*dim_phis(2)

      ## These are indexes for the auto-correlation matrix
      ii = ceil (i/dim_phis(1));
      jj = ceil (j/dim_phis(2));

      ## If we use the "reshape" command, we must keep in mind that it reshapes
      ## according to FORTRAN standard indexing (meaning that it progresses down
      ## the first column, then to the top of the second and down, then to the
      ## top of the third and down, etc.).  This is why we take the transpose here.
      PHI ( i:i+dim_phis(1)-1 , j:j+dim_phis(2)-1 ) = \
	  reshape(autoin ( abs(jj - ((max_lag-min_lag+1)+1)) + (ii-1),:), \
		   dim_phis(2), dim_phis(1) )';

    endfor
  endfor

  
  ## Calculate the cross-correlation functions for the output and input
  ## time series.  This is equivalent to the assignment to "autoin" above,
  ## except that the xcorr.m function does not accept a matrix for both
  ## input and output (X & Y).
  k = 1;

  for i=1:columns (its)
    for j=1:columns (ots)
      ## "flipud" is to fix problem that occurs since change in xcorr.m
      ## in octaveforge that was made in 2004 (EJR, 2/28/2005)
      cross_inout(:,k) = flipud(xcorr(its(:,i), ots(:,j), maxlagsize));
      #cross_inout (:,k) = xcov (its (:,i), ots(:,j), maxlagsize);
      k++;
    endfor
  endfor

  ## Shift for appropriate lags
  cross_inout = cross_inout( (maxlagsize+1) + min_lag: (maxlagsize+1) + max_lag , : );


  ## Arrange the cross-correlation matrix (Right-Hand side)
  dim_ccs = [columns(its),columns(ots)];
  CC = reshape(cross_inout(1,:), dim_ccs(2), dim_ccs(1) )';

  for i=2:(max_lag-min_lag+1)
    
    ## If we use the "reshape" command, we must keep in mind that it reshapes
    ## according to FORTRAN standard indexing (meaning that it progresses down
    ## the first column, then to the top of the second and down, then to the
    ## top of the third and down, etc.).  This is why we take the transpose here.
    CC = [CC ; reshape(cross_inout (i,:), dim_ccs(2), dim_ccs(1) )'];

  endfor


  ## Now use Single Value Decomposition to solve the system of equations
  [U,S,V] = svd (PHI);
  ## take out very small values of S
  S = diag (S);
  small = find(S/max(S) <= 1e-9);

  if small
    S(small) = 0;
  endif
  S = diag (S);

  ## Solve for the filter paramters (I don't remember exactly where I found
  ## this equation for solving this using the output from the SVD function,
  ## but it is essentially a form of reverse substitution)
  theta = V*(S \ (U' * CC) );
  lag = [min_lag:max_lag];

  ## ...and that's all folks!

endfunction
