### fir_siso_rls.m
###
### This routine calculates an adaptive FIR filter given "input" and 
### "output" time series, using a Recursive Least Squares algorithm.
### The name of this function implies that only FIRs are returned,
### but adaptive AR models can be returned as well, if the input
### and output vectors are the same, and the range of lag values
### does NOT include zero (if it included zero-lag, all of the response
### would be in the zero lag, and everything else would be zero!).
### This function is NOT designed to calculate IIR filters with both
### a numberator AND denominator.
### 
### Usage:
###
### [thetas, lags, errs] = fir_siso_rls (its, ots, lagsize,\
###                                      [theta0, lambda, delta])
###
### INPUTS:
###
### its       - Input time series.  For now this only handles single 
###             channel inputs.  Should be a column vector, but we'll
###             make an effort to fix it if it's not.
###
### ots       - Output time series.  For now this only handles single
###             channel outputs.  It actually shouldn't be all that
###             difficult to implement a multi-channel output (SIMO)
###             version.  Use column vectors.
###
### lagsize   - Either a scalar, or two element vector.  If a scalar, this
###             indicates the mimimum and maximum lag about zero lag.  If a
###             two element vector, lagsize(1) indicates the minimum lag
###             from zero, and lagsize(2) indicates the maximium lag from 
###             zero.  If the user is attempting to calculate an adaptive
###             AR filter, they MUST use the second form, and be certain
###             that all lags are greater than zero.
###
### theta0    - Allows the user to provide the initial value for the SISO
###             model.  If it's set to a scalar value, all values of theta0
###             will be set to that value.  If it is not set, an initial 
###             theta vector will be set to a vector of zeros.
###
### lambda    - Determines a "forgetting factor".  It defaults to .9, which
###             corresponds to an effective memory of ~10 samples.  .99 
###             corresponds to an effective memory of ~100 samples.
###
### delta     - The initial value for R is a delta*I, and delta is a small,
###             positive scalar.  This means that the initial Rinv value(s)
###             will be quite large.
### 
### OUTPUTS:
###
### thetas    - This 2-D matrix contains the model coefficients in rows
###             that correspond to each time step.  The length of thetas
###             will be equal to the length of the input/output time
###             series, and each filter will correspond to the point in
###             time coinciding with the zero-lag in the filter.
###
### lags      - A vector of integers corresponding to the time-lagged
###             model coefficients in each row of thetas.
###
### errs      - This column vector contains the error between ots, and
###             the modelled output, for time steps coinciding with the
###             zero lag of the filter(s).  Subtracting errs from the
###             appropriate ots time steps will result in
###             a vector of modelled output.  If this function is
###             used to determine an adaptive AR filter, (I believe) the
###             errs vector corresponds to an adaptively "pre-whitened"
###             input time series that can then be used to calculate an 
###             adaptive FIR model.
###

function [thetas, lags, errs, pred_err] = fir_siso_rls (input, output, lagsize, theta0, lambda, delta)

  ## We need to have at least three input parameters
  if (nargin < 3)
    error ("\nNot input enough parameters\n");
  endif

  ## For now this is only for SISO filters
  if (min (size (input))) > 1 || (min (size (output))) > 1
    error ("\nOnly Single Input / Single Output time series supported\n");
  endif

  ## Make certain input and output are column vectors
  if (columns (input) > 1)
    input = input';
  endif

  if (columns (output) > 1)
    output = output';
  endif

  ## Check of lambda and/or delta were set, or do we use defaults
  if (nargin > 6)
    error ("\nToo many input parameters\n");
  elseif (nargin == 3)
    theta0 = 0;
    lambda = .9;
    delta = .1;
  elseif (nargin == 4)
    lambda = .9;
    delta = .1;
  elseif (nargin == 5)
    delta = .1;
  endif

  ## alpha is simply 1-lambda...lambda and alpha define the "forgetting"
  ## factor
  alpha = 1-lambda;

  ## If lagsize is a scalar, assume that the FIR is length (2*lagsize+1)
  ## and centered at zero lag.  If lagsize is a 2-element vector, assume
  ## that the lag vector starts at lagsize(1), and ends at lagsize(2).
  ## This routine cannot currently handle lag vectors that are not
  ## simple indexes (i.e. sequential, increasing integers), so exit
  ## if lagsize is longer than 2 elements (i.e. don't trust the user
  ## to create his/her own lag vector properly).
  if (length (lagsize) > 2)
    error (["\nLag vector must be a scalar, or a 2-element vector\n", \
	    "that indicates the minimum and maximum lag indices.\n"]);
  endif

  if (length (lagsize) == 1)
    ## Assume a negative number is a mistake, and correct it
    lags = [-abs(lagsize):abs(lagsize)]';
  else
    if (lagsize (1) >= lagsize (2))
      error (["\nLag vector must be comprised of sequentially\n", \
	      "INCREASING integers.\n"]);
    endif
    lags = [lagsize(1):lagsize(2)];
  endif

  ## OK, now the input and output time series need to be of equal
  ## length, AND they must be at least as long as the "lags" vector
  if (length (input) != length (output) || \
      (length (lags) >= length (input) || \
       length (lags) >= length (output) ) )
    error (["\nInput and Output time series don't match,\n", \
	    "or they are not long enough to determine a \n", \
	    "filter recursively.\n"]);
  endif

  if ( is_scalar (theta0) )
    theta = ones (length (lags),1) * theta0;
  elseif ( is_vector (theta0) && length(theta0) == length(lags))
    theta = theta0;
  else
    error(["\nTheta0 must be either a scalar, or a vector of\n",\
	   "equal length to the lag vector.\n "]);
  endif

  ## Initialize Rinv
  Rinv = [eye (length (lags)) .* delta]^-1;


  ## OK, maybe NOW we can do some real work...

  ## The following is really unwieldy, and the reader needs to think
  ## hard if he/she wishes to understand.  This routine (currently)
  ## determines filter coefficients to predict the output at ZERO-
  ## lag.  This is easy to handle if we assume that the filter is 
  ## centered about zero-lag, or that it is always the same offset
  ## from zero-lag.  We want to be able to pass more-or-less arbitrary
  ## lag vectors to this routine, so we need a couple of variables to
  ## keep track of the shift from input's minimum index, to the zero
  ## lag in the output.  That's what the following does (there's 
  ## probably a MUCH better way to do this):
  ## 
  ## "mii" == "minimum input index"
  ##          (minimum index of input vector...note that these
  ##           indexs are not actually octave-style indexes, but
  ##           can be negative and zero...we'll shift things
  ##           appropriately when necessary)
  ## "zoi" == "zero output index"
  ##          (index of output vector that corresponds to zero lag)
  mii = max ([0, min (-lags)]);
  zoi = mii - min (-lags) ;

  ## Loop as many times as possible given 1) the length of the lag
  ## vector, and 2) the relative position of the output's zero lag
  ## to the min and max of the lag vector.
  errs = zeros ((length (input) - mii) - length (lags) + 1 - \
		max ([0,(zoi - (mii + length(lags)) + 1 ) ] ),1);
  thetas = zeros ((length (input) - mii) - length (lags) + 1 - \
		  max ([0,(zoi - (mii + length(lags)) + 1 ) ] ),length(lags) );

  
  ## The calculation of the following static values has been taken out
  ## of the loop in a very straight-forward attempt to speed things up
  ## slightly.
  l_lags = length(lags);

  for i = 1 : (length (input) - mii) - l_lags + 1 - \
    max ([0,(zoi - (mii + l_lags) + 1 ) ] )


    ## I think this was an attempt at a variable forgetting factor.
    ## I don't know if it worked well or not
    ##alpha = 1/sqrt((i+1))
    ##lambda = 1-alpha;

    Phi = input ( (mii+i) + l_lags - 1 :-1: (mii+i));
    ## just reverse index above, rather than a call to flipud
    ##Phi = flipud (Phi);
    yhat = theta' * Phi;
 
    err = output(zoi+i) - yhat;

    ## Probably remove this in favor of calculating the real prediction error
    errs(i) = err;

    Rinv = (1/lambda) * (Rinv - ((Rinv*Phi*Phi'*Rinv) / 
				 ((lambda/alpha) + Phi'*Rinv*Phi) ) );

    theta = theta + alpha*Rinv*Phi*err;
    thetas(i,:) = theta';

    pred_err(i) = output(zoi+i) - (theta' * Phi);


    printf ("\r %d of %d",i,(length (input) - mii) - l_lags + 1 - \
	    max ([0,(zoi - (mii + l_lags) + 1 ) ] ));
    fflush (stdout);
    #usleep (10000);
    #plot (lags,[theta0';theta']);
  endfor



  printf("\n");
  fflush (stdout);

endfunction
