## imp_gen_rls.m
##
## THIS ROUTINE NEEDS TO BE BROKEN UP INTO SEPARATE FUNCTIONS.
## IT SHOULD REALLY BE NO MORE THAN A WRAPPER SCRIPT.
##
## Adaptively pre-whitens imput and output time-series, then generates 
## a matrix of time-varying discrete impulse response functions. 
##

function [thetas_struct, tsteps, pred_errs, filt_lags] = imp_gen_rls (input, \
								      output, \
								      lagsize, \
								      lambda)
  
  if (nargin < 3 || nargin > 4)
    usage (["\n\n [thetas, tsteps, errs, lags] = ",\
	    "imp_gen_rls (input, output, lagsize [,lambda]\n"]);
  endif

  ## Only a single input vector is allowed for now
  if (! is_vector (input))
    error ("Only works with 1-D input vectors\n");
  endif

  ## Default
  loops = 1;
  ll = length (input);

  ## If subds is set, determine the length of the subdivisions of the
  ## input and output time series.
  if (nargin != 4 )
    lambda = .99; # assumes ~100 time step effective memory
  endif

  ## Make sure input and output vectors are same length
  if (rows (output) != length (input))
    error ("Input and Output dimenions don't match\n");
  endif

  ## If "lagsize" is a two-element vector, then assigne the max value of
  ## that vector to lagsize, and determine the min and max lag for use
  ## in calculating the correlation functions
  lags = zeros(1,2);
  if (is_vector(lagsize) && max(size(lagsize)) == 2)
    lags(1) = min(lagsize);
    lags(2) = max(lagsize);
    maxlagsize = max(abs(lagsize));
    ## if for some reason a vector was passed of identical values, assume that
    ## the user wants one negative and one positive
    if (lags(1) == lags(2))
      lags(1) = -lags(1);
    endif
  elseif (is_scalar(lagsize) ) ## do nothing
    min_lag = -lagsize;
    max_lag = lagsize;
    maxlagsize = lagsize;
  else
    error(["\n\"lagsize\" must be either a scalar, or a two-element vector\n",\
	   "indicating the min and max lag (relative to zero) for the correlation\n",\
	   "functions."]);
  endif


  ##
  ## Adaptively pre-whiten input time series with an adaptive AR filter.
  ## Save the vector of filters to do the initial pre-whitening of the
  ## output time series below.
  ##

  ## Calculate a static AR filter, then seed the adaptive algorithm
  ## The "pre-whitened" input vector is simply the error vector from
  ## fir_siso_rls.
  [AR_its_corls, lag_ar_its] = fir_mimo_corls (input-mean(input), \
					       input-mean(input), \
					       [1 ,abs(lags(2))]);
  [ARs_its_rls, lag_ars_its, its_prw] = fir_siso_rls (input-mean(input), \
						      input-mean(input), \
						      [1,abs(lags(2))], \
						      AR_its_corls, \
						      lambda);


  ## Loop over each of the output time series.
  for i=1:columns(output)

    ## First, pre-whiten the output time series with the input AR filter
    ## (note that the output here has shorted by the length of the AR filter,
    ##  it will do so again with the next pre-whitening)
    ots_rls = conv_siso_rls (output(:,i)-mean(output(:,i)), ARs_its_rls, lag_ars_its);
    ots_prw1 = [output(:,i)](abs(lags(2))+1:length(output(:,i))) - ots_rls;

    ## Next, calculate another AR filter based on the pre-whitened output,
    ## and pre-whiten the output a second time.  Theoretically this is not
    ## pre-whitening, since we don't refilter the input.  Think of it rather
    ## as the mapping of the diurnal variability in flux to the daily mean.
    ## (this thinking is probably very flawed, and is likely the reason
    ##  that the resulting impulse response functions seem to be significantly
    ##  diminished in absolute power)
    
    ## Calculate a static AR filter, then seed the adaptive algorithm.
    ## The "pre-whitened" data here is, again, simply the error vector
    ## returned by fir_siso_rls
    [AR_ots_prw_corls] = fir_mimo_corls (ots_prw1 - mean(ots_prw1), \
    					 ots_prw1 - mean(ots_prw1), \
    					 [1 ,abs(lags(2))]);
    [ARs_ots_prw_rls, lag_ars_ots, ots_prw2] = fir_siso_rls (ots_prw1-mean(ots_prw1), \
    							     ots_prw1-mean(ots_prw1), \
    							     [1,abs(lags(2))], \
    							     AR_ots_prw_corls, \
    							     lambda);

    ## If NOT commented, temporarily taking out second prewhitening...
    #ots_prw2 = ots_prw1 (abs(lags(2))+1:length(ots_prw1));

    
    ## At this point, we have its_prw, and ots_prw2 with which to calculate the 
    ## adaptive FIR filters.  They are different lengths, each shorted from its
    ## parent time series vector by the length of the AR filter.  We must be sure
    ## to pass the appropriately lagged vectors to the filter generator.  Now, 
    ## generate a static filter as a seed to the adaptive algorithm, then
    ## calculate the time-varying vilters.
    [filt_its_ots_corls] = fir_mimo_corls (its_prw (abs(lags(2))+1:length(its_prw)), \
					   ots_prw2, [lags(1),lags(2)]);
    [filt_its_ots_rls, lag_filt_iots, errs] = \
	fir_siso_rls (its_prw (abs(lags(2))+1:length(its_prw)), \
		      ots_prw2, [lags(1),lags(2)], filt_its_ots_corls, lambda);

    pred_errs (:,i) = errs;

    
    ## Save this 2-D matrix to a structure, since we can't do N-D arrays in octave.
    ## (structures require strings to identify their subelements, thus the need
    ##  for the "eval" below)

    ## Normally, this would all be good, but for 5 years of data, and 70 channels,
    ## and 166 filter coefficients each, and double precision numbers, I will run 
    ## out of memory (even with swap!).
    #if (nargout > 1)
    #
    #  if (columns(output) > 999)
    #	error ("Too many  matrices requested\n");
    #  endif
    #  
    #  eval (sprintf ("thetas_struct.filt_%04d = filt_its_ots_rls;",i) );
    #
    #endif

    ## This is not really a structure now...
    thetas_struct = filt_its_ots_rls;

    ## Print some kind of status meter
    printf("Finished %d'th filter.\n",i);
    fflush(stdout);


    tsteps = [length(AR_its_corls) + length(AR_ots_prw_corls) - lags(1) + 1: \
	      length(input) - lags(2)]';
    
    filt_lags = lag_filt_iots;
    
    
    eval (sprintf ("save -mat-binary imp_gen_rls_%04d.mat filt_its_ots_rls errs tsteps filt_lags;",i))

  endfor




  keyboard

endfunction
