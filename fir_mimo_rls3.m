### fir_mimo_rls3.m
###
### This routine calculates an "adaptive" FIR filter given input and 
### output time series, using a recursive least-squares (RLS) algorithm.
### The name of this function implies that only Finite Impulse Response
### functions are returned, but adaptive AR models can be returned as well
### if one passes appropriately lagged input and output time series. This 
### routine uses the algorithm described in Johansson (1993, pp. 263-267).
###
### 
### Usage:
###
### [theta, Pbar, lags, theta_mtx, Pbar_mtx, errs] = fir_mimo_rls3 \
###               (i_ts, o_ts, lagsize [, theta_0, Pbar_0, lambda, W]);
###
### INPUTS:
###
### i_ts      - Input time series.  Columns must represent different
###             input variables; rows represent observations made at
###             monotonically increasing moments in time.  
###
### o_ts      - Output time series.   Columns must represent different
###             output variables; rows represent observations made at
###             monotonically increasing moments in time.
###
### lagsize   - Either a scalar, or two element vector.  If a scalar, this
###             indicates the mimimum and maximum lag about zero lag.  If a
###             two element vector, lagsize(1) indicates the minimum lag
###             index, and lagsize(2) indicates the maximium lag index.  
###
###
###       Note: Positive lags correspond to observations made in the past;
###             negative lags correspond to observations made in the future.
###
###             Since this routine is designed to be used in an online, 
###             and possibly realtime configuration, it does not do 
###             any zero-padding, or pre-filtering to set nice initial 
###             conditions.  Therefore, the i_ts time series will always be 
###             longer than the o_ts time squotationeries by a number of samples
###             equal to the number of lags (defined by lagsize) minus 1. 
###
###             Furthermore, the first value of i_ts must lag in time the
###             first value of o_ts by a number of samples equal to the 
###             maximum lag, and the last value of i_ts must lag in time the
###             last value of o_ts by a number of samples equal to the 
###             minimum lag.  Hopefully, the following examples make this 
###             more clear:
###
###
###             lags = [-4  2]
###             o_ts  =                 [0 1 2 3 4 5]
###             i_ts  =           [-2 -1 0 1 2 3 4 5 6 7 8 9]
###
###             lags = [ 2  4]
###             o_ts  =             [0 1 2 3 4 5]
###             i_ts  = [-4 -2 -3 -1 0 1 2 3]
###
###             lags = [-4 -2]
###             o_ts  =                 [0 1 2 3 4 5]
###             i_ts  =                     [2 3 4 5 6 7 8 9]
###
###
###
### theta_0   - Allows the user to provide the initial value for the MIMO
###             model.  If it's set to a scalar value, all values of theta0
###             will be set to that value.  If it is a matrix of dimension
###             (lagsize * #inputs, #outputs), it is assumed that the model
###             is paramaterized according to:
###
###             theta = | A1_1 A2_1 ... AN_1 |
###                     | B1_1 B2_1 ... BN_1 |
###                     |  :    :    :   :   |
###                     | A1_2 A2_2 ... AN_2 |
###                     | B1_2 B2_2 ... BN_2 |
###                     |  :    :    :   :   |
###                     | A1_N A2_N ... AN_N |
###                     | B1_N B2_N ... BN_N |
###                     |  :    :    :   :   |
###
###
###             If it is a column vector of length (lagsize * #inputs * #outputs)
###             it is assumed that the model is parameterized according to:
###
###             theta = | A1_1 |        This is the parameterization used internally
###                     | A2_1 |        for this routine.  It can be easily derived
###                     |  :   |        from the previous, more natural parameter-
###                     | AN_1 |        ization by doing the following in Octave:
###                     |  :   |        
###                     | B1_1 |        octave>
###                     | B2_1 |        octave> theta_internal = theta'(:)
###                     |  :   |        octave>
###                     | BN_1 |
###                     |  :   |
###                     |  :   |
###                     | A1_2 |
###                     | A2_2 |
###                     |  :   |
###                     | AN_2 |
###                     |  :   |
###                     | B1_2 |
###                     | B2_2 |
###                     |  :   |
###                     | BN_2 |
###                     |  :   |
###                     |  :   |
###                     | A1_N |
###                     | A2_N |
###                     |  :   |
###                     | AN_N |
###                     |  :   |
###                     | B1_N |
###                     | B2_N |
###                     |  :   |
###                     | BN_N |
###                     |  :   |
###
###             (default value = 0)
###
### Pbar_0    - State error covariance.  Allows the user to provide the initial 
###             value for Pbar, or the inverse of the estimated Hessian matrix.  
###             This is so that the routine can be run "online". This input
###             parameter must be either 1) a scalar, which will be multiplied
###             by an identity matrix of a dimension equal to the total number
###             of state parameters (i.e. #lags * #inputs * #outputs),
###             2) a vector corresponding to the diagonal elements of the
###             covariance matrix, 3) a vector of length (#lags * #inputs * #outputs)^2
###             that can be reshaped into 4) a square, symmetric, true a priori
###             covariance matrix calculated using the user's favorite method.
###             (default value = 1)
###
### lambda    - Exponential forgetting factor.  This is a scalar fraction between
###             zero and one that is designed to prevent Pbar from getting
###             so small that new measurements no longer affect the state param-
###             eters.  The effective memory is 1/(1-lambda).  If this is equal to 
###             zero, the state is really nothing more than the ratio of inputs
###             to outputs at that time step.  If this is equal to 1, the algorithm
###             converges to the global least-squares solution.  A value of .99
###             has a memory of approximately 100 samles.
###             (default value = 1)
###
### W         - Observation weighting factor  If this parameter is a scalar or
###             a vector of length equal to the number of outputs, these values
###             will constitute the diagonal elements of a square matrix of
###             dimension (#outputs,#outputs).  If this parameter is a matrix
###             of this dimension, it is passed on "as is".  If this parameter is a 
###             matrix of dimension (length(o_ts), #outputs), it is assumed that
###             W is a set of time-depedent weigting vectors, each of which consti-
###             tute the diagonal elements of a square matrix of dimension #outputs.
###             (default value = 1)
###
###
### 
### OUTPUTS:
###
### theta     - This is the final set of model coefficients, parameterized according
###             the first form described above.  This is compatible with the form
###             returned by the least-squares solution given in fir_mimo_corls.m
###             and fir_mimo_regress.m.  It can also be passed back into this routine
###             "as is" if it is implemented online.
###
### Pbar      - This is the final covariance matrix.  It can be passed back into
###             this routine "as is" if it is implemented online.
###
### lags      - A vector of integers corresponding to the time-lags of
###             filter coefficient matrices in each element of thetas.
###
### errs      - This column vector contains the error between o_ts, and
###             the modelled output(s), for time steps coinciding with the
###             zero lag of the filter(s).  If one subtracts errs from the 
###             output data series, the result is the "predicted" output
###             at zero-lag.
###
### theta_mtx - This is a matrix of time-varying parameter vectors that are
###             parameterized according to the transpose of the second form described 
###             above.  This simply allows for compact storage until Octave truly 
###             implements N-D matrices. In order to convert a single row of this
###             matrix back into form 1, do the following:
###
###             octave>
###             octave> reshape (theta_mtx(i,:), size (theta'))';
###             octave>
###
### Pbar_mtx  - This is a matrix of reshaped, time-varying state-covariance
###             matrices.  Each row constitutes an entire covariance matrix.
###             In order to convert a single row back into a square covariance
###             matrix, do the following:
###
###             octave>
###             octave> reshape (Pbar_mtx(i,:), size(Pbar'))';
###             octave>
###


function [theta, Pbar, lags, errs, theta_mtx, Pbar_mtx] = fir_mimo_rls3 \
  (i_ts, o_ts, lagsize, theta0, Pbar0, lambda, W)

  ## We need to have at least three input parameters
  if (nargin < 3)
    error ("\nNot input enough parameters\n");
  endif


  ## Make sure the input and output time series matrices are composed of
  ## column vectors.  This is a bit of a kludge, since it assumes that the
  ## number of inputs/outputs is << length (i_ts).
  if (rows(i_ts) < columns(i_ts))
    i_ts = i_ts';
  endif
  if (rows(o_ts) < columns(o_ts))
    o_ts = o_ts';
  endif


  ## Do we use defaults?
  if (nargin > 7)
    error ("\nToo many input parameters\n");
  elseif (nargin == 3)
    theta0 = 0;
    Pbar0 = 1;
    lambda = 1;
    W = 1
  elseif (nargin == 4)
    Pbar0 = 1;
    lambda = 1;
    W = 1;
  elseif (nargin == 5)
    lambda = 1;
    W = 1;
  elseif (nargin == 6)
    W = 1;
  endif


  ## If lagsize is a scalar, assume that the FIR is length (2*lagsize+1)
  ## and centered at zero lag.  If lagsize is a 2-element vector, assume
  ## that the lag vector starts at lagsize(1), and ends at lagsize(2).
  ## This routine cannot currently handle lag vectors that are not
  ## simple indexes (i.e. sequential, increasing integers), so exit
  ## if lagsize is longer than 2 elements (i.e. don't trust the user
  ## to create his/her own lag vector properly).
  if (length (lagsize) > 2)
    error (["\nLag vector must be a scalar, or a 2element vector\n", \
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

  

  ## Make certain that the input and output time series are at least
  ## of the proper length, even if we can't guarantee that they are
  ## properly lagged.
  if ((rows (i_ts) - rows (o_ts)) != length(lags)-1 ||
      (length (lags) >= rows (i_ts)) )
    error (["\nInput and Output time series don\'t match,\n ", \
	    "or they are not long enough to determine a \n", \
	    "filter recursively.\n"]);
  endif
    


  ## A bunch of potentially useful static variables
  llags   = rows (lags);          ## number of lags
  ni_ts   = columns (i_ts);       ## number of different inputs
  no_ts   = columns (o_ts);       ## number of different outputs
  lo_ts   = rows (o_ts);          ## number of output data points
  li_ts   = rows (i_ts);          ## number of input data points

  ## Vectors and/or matrices should be initialized if possible
  yhat = zeros (no_ts,1);
  PHI = zeros (llags*ni_ts,1);

  ## These are matrices designed to store the error, state, and co-
  ## variance matrix at each time step.  Only initialize them if the
  ## user wants them returned because they suck up a LOT OF MEMORY!
  if (nargout > 3)
    errs = zeros (lo_ts,no_ts);
  endif
  if (nargout > 4)
    theta_mtx = zeros (lo_ts,llags*ni_ts*no_ts );
  endif
  if (nargout > 5)
    Pbar_mtx = zeros (lo_ts, (llags*ni_ts*no_ts)^2);
  endif


  ## If an initial value theta0 was provided, 1) is it the correct size?
  ## 2) convert scalars to properly sized matrix of the same value; and
  ## 3) copy the value to the first element of theta.
  if (is_scalar (theta0))
    theta = theta0 * ones (columns(i_ts)*rows(lags), no_ts);
  elseif (is_vector (theta0))
    if (length (theta0) != (ni_ts*llags*no_ts))
      error (["\nFilter's initial value has wrong dimensions!\n"]);
    endif
    ## This is horrible, because we will simply reshape it again later,
    ## but for now parameterize theta according to form #1 in the Help section.
    theta = reshape (theta0, [no_ts, ni_ts*llags])';
  elseif (is_matrix (theta0))
    if (rows (theta0) / ni_ts != llags || \
	columns (theta0) != no_ts )
      error (["\nFilter's initial value has wrong dimensions!\n"]);
    endif
    theta = theta0;
  else
    error (["\nFilter is wrong data-type, it must be a scalar or properly\n",\
	    "sized vector or matirx.\n"]);
  endif


  ## Initialize Pbar
  if (is_scalar (Pbar0))
    Pbar = [eye (ni_ts*llags*no_ts) .* Pbar0];
  elseif (is_vector (Pbar0))
    if (length (Pbar0) == (ni_ts*llags*no_ts))
      Pbar = diag (Pbar0);
    elseif (length (Pbar0) == (ni_ts*llags*no_ts)^2)
      Pbar = reshape (Pbar0, (ni_ts*llags*no_ts), (ni_ts*llags*no_ts))';
    else
      error(["\nCovariance's initial value has wrong dimension!\n"]);
    endif

  elseif (is_matrix (Pbar0))
    if (size (Pbar0) != size (eye(ni_ts*llags*no_ts)) )
      error (["\nWrong size Covariance matrix!\n"]);
    endif
    Pbar = Pbar0;
  else
    error (["\nCovariance matrix is wrong data-type, it must be a scalar or ",\
	    "properly sized vector or matirx.\n"]);
  endif


  ## Make sure lambda is appropriate.
  if (!is_scalar (lambda) ||
      lambda > 1 ||
      lambda < 0)
    error (["\"lambda\" must be a scalar betwen zero and one.\n"]);
  endif


  ## Make sure the weighting matrix W is set up properly.
  tv_W = 0;
  if (is_scalar (W))
    W = eye (no_ts) * W;
  elseif (is_vector (W))
    if (length(W) == no_ts)
      ## This pertinent for multiple outputs, with constant weighting
      W = diag (W);
    elseif (size(W) == size (o_ts))
      ## This is pertinent for single outputs, with time-varying weighting.
      ## Change nothing in W, but set a flag so that the algorithm knows
      ## how to handle it later on.
      tv_W = 1;
    else
      error (["\n\"W\" is not an appropriate dimension, read the help file.\n"]);
    endif
  elseif (size (W) == size (o_ts))
    ## this is pertinent for multiple outputs, with time-varying weighting.
    ## Change nothing in W, but set a flag so that the algorithm knows how
    ## to handle it later on.
    tv_W = 1;
  elseif (size (W) == size (eye(no_ts)))
    ## I can find nothing in the literature that says it is OK to treat this
    ## as anything but a set of weights.  This means that anything on the
    ## off-diagonals is not appropriate
    if (W != diag (diag (W)))
      error (["\n\"W\" can have no off-diagonal values.\n"]);
    endif
  else
    error (["\"W\" is not an appropriate dimension, read the help file.\n"]);
  endif
      




  ## OK, maybe NOW we can do some real work...  

  ## Reshape theta.  This will need to be switched back at the end.
  theta_tmp = theta'(:);

  ## Begin the main loop
  for i=1:lo_ts


    ## Determine PHI for this time step
    for j=1:ni_ts:llags*ni_ts
      
      ##printf ("\r i=%d  j=%d",i,j);
      ##fflush (stdout);
      
      jj = ceil (j/ni_ts);
      PHI (j:j+ni_ts-1) = i_ts (abs(jj - (llags + 1) ) + (i-1), :)';

    endfor


    ## This converts PHI into a form compatible with the parameterization
    ## of theta described by form #2 in the Help section.
    PHI_tmp = kron (PHI,eye(no_ts));

    yhat = PHI_tmp' * theta_tmp; 
    err = o_ts (i,:)' - yhat;


    ## Check whether W is time-varying or not.
    if (tv_W == 1)
      W_tmp = diag (W(i,:));
    else
      W_tmp = W;
    endif


    ## This is mostly to suppress warnings about singular matrices,
    ## but maybe someday I'll actually check the "rcond".
   [denom,rcond] = inverse (PHI_tmp' * Pbar * PHI_tmp + lambda * (W_tmp^-1) );

    gamma = Pbar * PHI_tmp * denom;
    Pbar = (1 / lambda) * (Pbar - gamma * PHI_tmp' * Pbar);
    theta_tmp = theta_tmp + gamma * err;


    ## Store errs, state, and covariance matrices for each time step.
    ## Only store these if the user wants them returned, otherwise they
    ## suck up a WHOLE LOT OF MEMORY (especially the covariances)!
    if (nargout > 3)
      errs (i,:) = err';
    endif
    if (nargout > 4)
      theta_mtx(i,:) = theta_tmp';
    endif
    if (nargout > 5)
      Pbar_mtx(i,:) = reshape (Pbar', 1, rows (Pbar)*columns(Pbar));
    endif


    ## This probably slows things down slightly, but that's OK for now.
    printf ("\r %d of %d",i,lo_ts);
    fflush (stdout);


  endfor

  ## reshape the last state vector into the proper form
  theta = reshape (theta_tmp, size (theta'))';
  

  ## Just giving the user a new line for the next interactive prompt.
  printf("\n");
  fflush (stdout);


endfunction
