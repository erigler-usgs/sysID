### fir_mimo_rls2.m
###
### This routine calculates an adaptive FIR filter given "input" and 
### "output" time series, using a Recursive Least Squares algorithm.
### This routine uses the algorithm described in O. Nelles (2001).
### The name of this function implies that only FIRs are returned,
### but adaptive AR models can be returned as well, if the input
### and output vectors are the same, and the range of lag values
### does NOT include zero (if it included zero-lag, all of the response
### would be in the zero lag, and everything else would be zero!).
###
### 
### Usage:
###
### [thetas, lags, errs, Pbar] = fir_mimo_rls2 (its, ots, lagsize,\
###                                             [theta0, lambda, delta])
###
### INPUTS:
###
### its       - Input time series.  Columns must represent different
###             input variables; rows must represent observations at
###             identical, monotinically increasing moments in time.
###
### ots       - Output time series.   Columns must represent different
###             output variables; rows must represent observations at
###             identical, monotinically increasing moments in time
###             corresponding to the input times.
###
### lagsize   - Either a scalar, or two element vector.  If a scalar, this
###             indicates the mimimum and maximum lag about zero lag.  If a
###             two element vector, lagsize(1) indicates the minimum lag
###             from zero, and lagsize(2) indicates the maximium lag from 
###             zero.  
###
### theta0    - Allows the user to provide the initial value for the MIMO
###             model.  If it's set to a scalar value, all values of theta0
###             will be set to that value.  If it is not set, an initial 
###             theta vector will be set to a vector of zeros.
###
### Pbar0     - Allows the user to provide the initial value for Pbar, or the
###             inverse of the estimated Hessian matrix.  This is so that the
###             routine can be run "online".  If it is a scalar, it is then
###             its inverse is multiplied by the identity matrix.  If it's
###             not provided, the default value is .1.
###
### lambda    - Determines a "forgetting factor".  It defaults to .9, which
###             corresponds to an effective memory of ~10 samples.  .99 
###             corresponds to an effective memory of ~100 samples.
###
###
### 
### OUTPUTS:
###
### thetas    - Structure containing matrix filters for each time step.
###             "theta_00000" is the "initial guess".  "theta_xxxxx"
###             is the filter matrix for each time step.  Until
###             decent N-D matrix support is available in Octave,
###             this format will limit the length the number of
###             time steps to be modeled to 99,999 (probably just as
###              well, otherwise we really push memory limits).
###
### lags      - A vector of integers corresponding to the time-lags of
###             filter coefficient matrices in each element of thetas.
###
### errs      - This column vector contains the error between ots, and
###             the modelled output(s), for time steps coinciding with the
###             zero lag of the filter(s).  If one subtracts errs from the 
###             output data series, the result is the "predicted" output
###             at zero-lag.
###
### Pbar      - The Hessian matrix at the last valid time-step.  This is
###             primarily to use as an initial value when making subsequent
###             calls to this routine.
###

function [thetas_mtrx, lags, errs, Pbar] = fir_mimo_rls2 (its, \
							  ots, \
							  lagsize, \
							  theta0, \
							  Pbar0, \
							  lambda)

  ## We need to have at least three input parameters
  if (nargin < 3)
    error ("\nNot input enough parameters\n");
  endif

  ## TESTING TRUE MIMO MODELS
  ##
  ## For now this is only for MISO filters
  #if (min (size (ots))) > 1
  #  error ("\nOnly Single Output time series supported\n");
  #endif

  ## Make sure the input and output time series matrices are composed of
  ## column vectors.  This is a bit of a kludge, since it assumes that the
  ## number of inputs/outputs is << length (its).
  if (rows(its) < columns(its))
    its = its';
  endif
  if (rows(ots) < columns(ots))
    ots = ots';
  endif

  ## Do we use defaults?
  if (nargin > 7)
    error ("\nToo many input parameters\n");
  elseif (nargin == 3)
    theta0 = 0;
    Rinv0 = .1;
    lambda = .9;
  elseif (nargin == 4)
    Rinv = .1;
    lambda = .9;
  elseif (nargin == 5)
    lambda = .9;
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
    lags = [lagsize(1):lagsize(2)]';
  endif

  ## The input and output time series need to be of equal
  ## length, AND they must be at least as long as the "lags" vector
  if (rows (its) != rows (ots) || \
      (rows (lags) >= rows (its) || \
       rows (lags) >= rows (ots) ) )
    error (["\nInput and Output time series don't match,\n", \
	    "or they are not long enough to determine a \n", \
	    "filter recursively.\n"]);
  endif


  ## A bunch of potentially useful static variables
  llags  = rows (lags);    ## number of lags
  nits   = columns (its);  ## number of different inputs
  nots   = columns (ots);  ## number of different outputs
  loits  = rows (its);     ## number of input/output data points
  ltheta = llags*nits;     ## length of theta 
  wtheta = nots;           ## width of theta

  ## Vectors and/or matrices should be initialized if possible
  yhat = zeros (1,nots);
  errs = zeros (loits,nots);
  PHI  = zeros (ltheta,1);


  ## If an initial value theta0 was provided, 1) is it the correct size?
  ## 2) convert scalars to properly sized matrix of the same value; and
  ## 3) copy the value to the first element of thetas (i.e. theta_00000)
  if (is_scalar (theta0))
    thetas.theta_00000 = theta0 * ones (columns(its)*rows(lags), nots);
  else 
    if (rows (theta0) / nits != llags ||
	columns (theta0) != nots )
      error (["\nFilter's initial value has wrong dimensions!\n"]);
    endif
    thetas.theta_00000 = theta0;
  endif
  ## This is the variable for the loop -- initializing it here is really 
  ## kind of pointless.
  theta = thetas.theta_00000;


###
###
### THIS IS DIFFICULT...HOW DOES ONE DEFINE Pbar FOR MIMO FILTER COEFFICIENTS?
###
### Answer:  You don't!  (for now)  We will do a MISO filter.  I _believe_ that
###          the solution may be simply to reshape the MIMO filter into a
###          single, long vector (didn't we do this in Stat-OD), but I don't have
###          time to check this out.
###
###          Actually, I think that Pbar will be OK as is.  I am going to check
###          this out...I was correct.  Eventually, I should remove this comment
###          block, but leave it for now so I don't forget.
###


  ## Initialize Pbar
  ## This isn't really necessary since discovering that saving a structure
  ## of relatively large square matrices can REALLY eat up memory, and that
  ## this structure really isn't all that useful anyway.  I'm just leaving
  ## it for now, and removing the update from the loop.  Pbar will just be
  ## a single matrix returned at the end.
  if (is_scalar (Pbar0))
    Pbars.pbar_00000 = [eye (ltheta) .* Pbar0]^-1;
  else
    if (size (Pbar0) != size (eye(ltheta)) )
      error (["Wrong size Hessian matrix!\n"]);
    endif
    Pbars.pbar_00000 = Pbar0;
  endif
  ## This is the variable for the loop -- initializing it here is really 
  ## kind of pointless.
  Pbar = Pbars.pbar_00000;
  
  ## We'll need this later
  EYE = eye (ltheta);


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
  moi = max ([0, max (-lags)]);
  zoi = mii - min (-lags) ;

  ## Now, we need to pad the input (its) and output (ots) time series
  ## with zeros.  This means that parameters at the beginning and end
  ## will not be very reliable.  I should attempt to determine the
  ## proper initial conditions, but I'm not sure how to do this with
  ## MISO/MIMO models, or with models with negative lags.  Zeros it 
  ## is then!
  empty_list_elements_ok = true;
  its = [zeros (zoi, nits); its; zeros (max([0,llags-zoi-1]), nits)];
  ots = [zeros (zoi, nots); ots; zeros (max([0,llags-zoi-1]), nots)];
  empty_list_elements_ok = false;


  ## OK, maybe NOW we can do some real work...  

  ## Begin the main loop
  for i=1:loits


    ## Determine PHI for this time step
    for j=1:nits:ltheta
      
      ##printf ("\r i=%d  j=%d",i,j);
      ##fflush (stdout);
      
      jj = ceil (j/nits);
      
      PHI (j:j+nits-1) = its (abs(jj - (ltheta/nits + 1) ) + (i-1), :)';
    endfor


    yhat = PHI' * theta;
 
    err = ots (zoi+i,:) - yhat;


    ## Only update these if the zero-padding is not being used
    ## (this should be some kind of option)
    if (i > (max(abs([zoi,mii]) ) ) && \
	i <= (loits - max(abs ([zoi,moi]) ) ) )


      gamma = (1 / (PHI' * Pbar * PHI + lambda) ) * Pbar * PHI;
      
      Pbar = (1 / lambda) * (EYE - gamma * PHI') * Pbar;
      
      theta = theta + gamma * err;
      
    endif
    

    ## For now, with MISO filts, we just make a 2-D matrix, then
    ## save to a named file...
    ##eval (sprintf ("thetas.theta_%05d = theta;",i));
    ##eval (sprintf ("Rinvs.rinv_%05d = Rinv;",i));


    ## Need to reshape this matrix for storage when we don't have N-D matrix
    ## support.
    ##thetas_mtrx (i,:) = theta';
    thetas_mtrx = theta';
    errs (i,:) = err;

    printf ("\r %d of %d",i,loits);


    fflush (stdout);
    ##usleep (10000);

    ##grid ("on");
    ##for k=1:nits
    ## 
    ##  figure (k);
    ##  plot (lags,[thetas.theta_00000(k:nits:rows(theta)),theta(k:nits:rows(theta))]);
    ## 
    ##endfor
    ##pause

  endfor

  printf("\n");
  fflush (stdout);

endfunction
