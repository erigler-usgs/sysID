###
### fir_mimo_regress
###
###
### This routine will produce a coupled Finite Impulse Response (FIR) linear 
### filter based on given input and output time series.  This filter is the 
### solution to system of linear equations arranged in a "regression matrix."
###
###
### Usage:
###
### [theta, lags, lhs, rhs] = fir_mimo_regress (its, ots[, lagsize])
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
### theta     - vector of optimal filter coefficients
### lags      - vector of corresponding lag values (units are whatever the
###             sampling interval is)
### PHI       - regression matrix that, when it is multiplied
###             with the theta matrix, one gets the in-sample prediction
###             that minimized the least-squares criterion.
###

function [theta ,lags , PHI] = fir_mimo_regress (its, ots, lagsize)

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
  if (rows(its)<columns(its))
    its = its';
  endif
  if (rows(ots) < columns(ots))
    ots = ots';
  endif


  ## Determine the number of inputs and outputs
  nits = columns (its);
  nots = columns (ots);



  ## Determine the length and width of the regression matrix.
  ##
  ## The problem here is that we waste computations if our lags
  ## are not symmetric about zero because we could technically
  ## have a smaller regression matrix.  For now, however, we
  ## will simply calculate symmetric, acausal filters, then
  ## extract the lags we want at the end.  ** NOTE:  not only
  ## is this more computationally expensive, it is not the 
  ## optimal set of parameters!  We REALLY should fix this to
  ## set up the proper regression matrix for lag times that
  ## are not symmetric about zero.
  lphi = rows(its);
  wphi = nits * (2*maxlagsize+1);



  ## Pad the input time series with zeros so that when we
  ## generate PHI below, the zero-lags line up properly
  ## with the output time series (this was originally
  ## done in the EJR function filt_mimo.m...there is 
  ## probably a better way to do this).
  its = [zeros((floor((wphi/nits)/2)),nits);\
	 its;\
	 zeros((floor((wphi/nits)/2)),nits)];


  ## This is so warnings about "empty lists" don't occur when we
  ## make a call to the 'shift' function with a zero-length shift
  empty_list_elements_ok = true;

  ## Initialize the regression matrix
  PHI = zeros (lphi,wphi);

  ## Fill the regression matrix
  
  for j=1:nits:wphi
    
    PHI (:,wphi-j+1-nits+1 : wphi-j+1) = \
	shift (its,-(j-1)/nits) (1:lphi,:);
    
  endfor

  ## Now, calculate the pseudo-inverse of the regression matrix
  ## using SVD.  I'm not entirely sure why I use SVD here, but
  ## it still works.

  PHI2 = PHI'*PHI;
  [U,S,V] = svd (PHI2);

  ## Remove very small values of S (1e-9 is fairly arbitrary)
  S = diag (S);
  small = find (S/max(S) <= 1e-9);
  if small
    S(small) = 0;
  endif
  S = diag (S);
  
  ## This is the actual pseudo-inverse
  PHI2_pi = V * inv(S) * U';

  ## Finally, calculate the optimal parameters.
  theta = PHI2_pi * PHI' * ots;
  

  ## Even though this isn't technically correct for non-symmetric
  ## lag-times, extract the coefficients for the lags we want to
  ## return.

  lags = [min_lag:max_lag];

  idx1 =  (maxlagsize*nits+1) + (min_lag*nits);
  idx2 =  (maxlagsize*nits+1) + (max_lag*nits + nits - 1);

  theta = theta (idx1:idx2);
  
  
endfunction

