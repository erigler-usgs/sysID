## eureka.m
##
## This routine will determine the linear prediction filter coefficients
## for a single input / single output (SISO) system given the auto- and
## cross-correlations for the respective input and output time series.
## It uses a recursive algorithm described in "Multichannel Time Series
## Analysis" (Robinson, 1983) by a similar name.
##
## Inputs
##
## 

function [f, a] = eureka (lf, rin, gin)

  if (nargin < 3 || nargin > 3)
    usage (["\n\n[filt, err] = eureka ( length, auto, cross)\n",\
	    "\n (where \"length\" is a scalar describing a symetric\n",\
	    "  filter of maxlag=length, OR it is a 2-element vector\n",\
	    "  describing the minimum and maximum lag of the desired\n",\
	    "  output filter)"]);
  endif

  ## Check first input parameter
  if (is_scalar(lf))
    minlag = -lf;
    maxlag = lf;
  elseif (is_vector(lf) && (length(lf)==2) )
    minlag = lf(1);
    maxlag = lf(2);
  else
    error ("First input parameter must be a scalar, or two-element vector\n");
  endif

  ## Check second and third input parameters
  if (is_vector(rin) && is_vector(gin))
    if (length(rin) != length(gin))
      error ("Second and Third input parameters must be equal length vectors\n");
    endif
  else
    error ("Second and Third input paramters must be equal length vectors\n");
  endif

  ## Since the variable rin is supposed to be an autocorrelation function, it
  ## should be symetric, with a peak at the zero-lag, and of an odd length.
  ## This is consistent of the output from the function xcov.  We can only
  ## assume that the variable gin was created in a similar fashion, but there
  ## is no guarantee that its peak is at zero-lag.
  if( (length(rin)/2) == (floor(length(rin)/2)) ) ||
    ( max(rin) != rin(length(rin)/2) )
    error(["Autocorrelation vector should be symetric, and of odd length,\n",\
	   "with its peak at zero-lag."]);
  endif

  ## Check to see if the bounds of the desired filter lie outside the range
  ## of the correlation functions passed in
  if (max (abs([minlag,maxlag]) ) > length (rin)/2 -1 )
    error ("Requested filter is outside bounds of correlation functions\n");
  endif


  ## That's about all the checking we can do, move on...


  ## Initialize the output paramters
  f = zeros (maxlag-minlag+1,1);
  a = f;

  ## Extract appropriate sections of the correlation functions, and
  ## rearrange r so that it would form a toplitz matrix.
  
  GAVE UP HERE, I REALIZED THAT THIS ALGORITHM CANNOT WORK THE WAY I WANT IT TO!!!
  


  ## Starting here is a nearly line-for-line translation of the FORTRAN
  ## subroutine "eureka" found in Robinson (1983).  It may be difficult
  ## to follow
  
  v = r(1);
  d = r(2);

  a(1) = 1;
  f(1)


endfunction
