## fir_siso_rls
##
## Calculates adaptive FIR filter given "input" and "output" time
## series, using a Recursive Least Squares algorithm.  


function [fir_mtrx] = fir_siso_rls (input, output, lags, lambda, delta)

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

  ## Check of lambda and/or delta were set, or do we use defalts
  if (nargin > 5)
    error ("\nToo many input parameters\n");
  elseif (nargin == 3)
    lambda = .9;
    delta = .1;
  elseif (nargin == 4)
    delta = .1;
  endif
  
  ## alpha is simply 1-lambda...lambda and alpha define the "forgetting"
  ## factor
  alpha = 1-lambda;

  ## If lags is a scalar, assume that the FIR is length (2*lags+1)
  ## and centered at zero lag.  If lags is a 2-element vector, assume
  ## that the lag vector starts at lags(1), and ends at lags (2).
  ## This routine cannot currently handle lag vectors that are not
  ## simple indexes (i.e. sequential, increasing integers), so exit
  ## if lags is longer than 2 elements (i.e. don't trust the user
  ## to create his/her own lag vector properly).
  if (length (lags) > 2)
    error (["\nLag vector must be a scalar, or a 2-element vector\n", \
	    "that indicates the minimum and maximum lag indices.\n"]);
  endif

  if (length (lags) == 1)
    ## Assume a negative number is a mistake, and correct it
    lag = [-abs(lags):abs(lags)]';
  else
    if (lags (1) >= lags (2))
      error (["\nLag vector must be comprised of sequentially\n", \
	      "INCREASING integers.\n"]);
    endif
    lag = [lags(1):lags(2)];
  endif

  ## OK, now the input and output time series need to be of equal
  ## length, AND they must be at least as long as the "lag" vector
  if (length (input) != length (output) || \
      (length (lag) >= length (input) || \
       length (lag) >= length (output) ) )
    error (["\nInput and Output time series don't match,\n", \
	    "or they are not long enough to determine a \n", \
	    "filter recursively.\n"]);
  endif



  ## OK, maybe NOW we can do some real work...


  ## For now we'll just hard code whether our initial filter is
  ## a zero vector, or seed it with the static filter
  theta = zeros (length (lag),1);

  ## Initialize Rinv
  Rinv = eye (length (lag)) .* delta;


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
  mii = max ([0, min (-lag)]);
  zoi = mii - min (-lag) ;

  ## Loop as many times as possible given 1) the length of the lag
  ## vector, and 2) the relative position of the output's zero lag
  ## to the min and max of the lag vector.
  for i = 1 : (length (input) - mii) - length (lag) + 1 - \
    max ([0,(zoi - (mii + length(lag)) + 1 ) ] )

    ##alpha = 1/sqrt((i+1))
    ##lambda = 1-alpha;

    Phi = input ( (mii+i) : (mii+i) + length (lag) - 1 );
    Phi = flipud (Phi);
    yhat = theta' * Phi;
 
    err = output(zoi+i) - yhat;

    Rinv = (1/lambda) * (Rinv - ((Rinv*Phi*Phi'*Rinv) / 
				 ((lambda/alpha) + Phi'*Rinv*Phi) ) );

    theta = theta + alpha*Rinv*Phi*err;

    printf ("\r %d",i);
    fflush (stdout);
    usleep (100000);
    plot (lag,theta);
  endfor

  keyboard

endfunction
