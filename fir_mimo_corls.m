###
### fir_mimo_corls.m
###
### This routine will produce a simple Finite Impulse Response (FIR) linear 
### filter based on given input and output time series.  This filter is the 
### solution to system of linear equations derived from the auto- and cross-
### correlation vectors, rather than the actual data-based regression matrix.
### This should, in theory, result in a "consistent", if not necessarily un-
### biased estimate of the optimal parameters (that's what was said in 
### "Nonlinear System Identification", Nelles (2001) about a nearly identical 
### algorithm designed to determine ARX model parameters...since this model 
### has no recursion, I think it is also guaranteed to be "unbiased").
###
### Usage:
###
### [theta] = fir_mimo_corls (its, ots[, lagsize])
###
### filt     - vector of optimal filter coefficients
### lag      - vector of corresponding lag values (units are whatever the
###            sampling interval is)
### crosscov - full cross-correlation vector for input/output time series
### autoin   - full auto-correlation vector for the input time series
### autoout  - full auto-correlation vector for the output time seriess
###

function [theta,lag] = fir_mimo_corls (its, ots, lagsize)

  if ( nargin < 2 || nargin > 3)
    usage (["\n\n[filter, lag] = ", \
	    "fir_mimo_corls (input_ts, output_ts [, lagsize])\n"]);
  endif

  ## Set a default lagsize of 10
  if nargin == 2
    lagsize = 10;
  endif

  ## Make sure the input and output time series matrices are composed of
  ## column vectors.  This is a bit of a kludge, since it assumes that the
  ## length(theta)<<length(its).
  if (rows(its)<columns(its))
    its = its';
  endif
  if (rows(ots) < columns(ots))
    ots = ots';
  endif

  ## Calculate auto-correlation functions for input series
  ## (we make it twice the lagsize, so that we can get a square
  ##  regression matrix when calculating symetric filters about
  ##  zero lag)
  autoin = xcorr (its,lagsize*2);

  ## Arrange the auto-correlation regression matrix (Left-Hand side)
  dim_phis = [columns(its),columns(its)];
  PHI = zeros ((2*lagsize+1) * dim_phis(1));

  for i=1:dim_phis(1):(2*lagsize+1)*dim_phis(1)
    for j=1:dim_phis(2):(2*lagsize+1)*dim_phis(2)

      ## These are indexes for the auto-correlation matrix
      ii = ceil (i/dim_phis(1));
      jj = ceil (j/dim_phis(2));

      ## If we use the "reshape" command, we must keep in mind that it reshapes
      ## according to FORTRAN standard indexing (meaning that it progresses down
      ## the first column, then to the top of the second and down, then to the
      ## top of the third and down, etc.).  This is why we reverse the indices, 
      ## then take the transpose here.
      PHI ( i:i+dim_phis(1)-1 , j:j+dim_phis(2)-1 ) = \
	  reshape (autoin ( abs(jj - ((2*lagsize+1)+1)) + (ii-1),:), \
		   dim_phis(2), dim_phis(1) )';
    endfor
  endfor

  ## Calculate the cross-correlation functions for the output and input
  ## time series.  This is equivalent to the assignment to "autoin" above,
  ## except that the xcov.m function does not accept a matrix for both
  ## input and output (X & Y).  Also, according to (Robinson, 1983)
  k = 1;

  for i=1:columns (ots)
    for j=1:columns (its)

      cross_inout (:,k) = flipud (xcorr (ots (:,i),its(:,j), lagsize) );
      k++;

    endfor
  endfor

  ## Arrange the cross-correlation matrix (Right-Hand side)
  dim_ccs = [columns(its),columns(ots)];
  CC = reshape (cross_inout (1,:), dim_ccs(2), dim_ccs(1) )';

  for i=2:(2*lagsize+1)
    
    ## If we use the "reshape" command, we must keep in mind that it reshapes
    ## according to FORTRAN standard indexing (meaning that it progresses down
    ## the first column, then to the top of the second and down, then to the
    ## top of the third and down, etc.).  This is why we reverse the indices, 
    ## then take the transpose here.
    CC = [CC ; reshape (cross_inout (i,:), dim_ccs(2), dim_ccs(1) )'];

    ## (That was a heck of a lot easier than setting up the 
    ##  auto-correlation matrix!!!)

  endfor

  ## Now use Single Value Decomposition to solve the system of equations
  [U,S,V] = svd (PHI);
  ## take out very small values of S
  S = diag (S);
  small = find (S/max(S) <= 1e-9);

  if small
    S(small) = 0;
  endif
  S = diag (S);

  ## Solve for the filter paramters (I don't remember exactly where I found
  ## this equation for solving this using the output from SVD, but it is
  ## essentially a form of reverse substitution)
  theta = V*(S \ (U' * CC) );
  lag = [-lagsize:lagsize];

  ## ...and that's all folks!

endfunction
