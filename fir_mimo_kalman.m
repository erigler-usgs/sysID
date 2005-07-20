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
### [theta, Pbar, lags, errs, theta_mtx, Pbar_mtx] = fir_mimo_kalman \
###                 (i_ts, o_ts, lagsize [, theta_0, Pbar_0, R1, R2]);
###
###
### Calculates an "adaptive" multi-input/multi-output FIR filter given 
### input and output time series, using a form of Kalman filter algorithm.
### The routine uses an algorithm described in Johansson (System Modeling
### and Identification, 1993; pp. 120;269), which is very similar to the 
### algorithm described in Nelles (Nonlinear System Identifications, 2001;
### pp. 66) but we designate an observation error covariance matrix rather 
### than observation weights.
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
###             longer than the o_ts time series by a number of samples
###             equal to the number of lags defined by lagsize minus 1. 
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
###             model.  If it's set to a scalar value, all values of theta_0
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
###                     | AN_1 |        ization in Octave by doing the following:
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
###             (default = zeros(#lags*#inputs,#outputs))
###
### Pbar_0    - State error covariance.  Allows the user to provide the initial 
###             value for Pbar, or the inverse of the estimated Hessian matrix.  
###             This is so that the routine can be run "online". This input
###             parameter must be either 1) a scalar, which will be multiplied
###             by an identity matrix of a dimension equal to the total number
###             of state parameters (i.e. #lags * #inputs * #outputs),
###             2) a vector corresponding to the diagonal elements of the
###             covariance matrix, or 3) a square, symmetric, true a priori
###             covariance matrix calculated using the user's favorite method.
###             (default = 1)
###
### R1        - State or process noise covariance.  This is a vector or matrix 
###             describing the strength of the time variance of the parameters
###             in theta, and can be tuned to affect stability and convergence
###             properties. This input parameter must be either 1) a scalar, which
###             will be multiplied by an identity matrix of a dimension equal 
###             to the number of state parameters (i.e. #lags * #inputs * #outputs),
###             2) a vector corresponding to the diagonal elements of the
###             covariance matrix, 3) a square, symmetric, true a priori
###             covariance matrix calculated using the user's favorite method.
###             (Note that this matrix is designed to operate on a state/parameter
###              in the second form described in the help section), or 4) a matrix
###             with a number of columns equal to the number of state parameters
###             and a number of rows equal to the number of output observations.
###             The last is a form that can be used to have time-dependent diagonal
###             process noise covariance matrix; eventually this routine will allow 
###             for full time-dependent process noise covariance matrices.
###             (default = 0)
###
### R2        - Measurement noise covariance.  If this parameter is a scalar or
###             a vector of length equal to the number of outputs, these values
###             will constitute the diagonal elements of a square matrix of
###             dimension (#outputs,#outputs).  If this parameter is a matrix
###             of this dimension, it is passed on.  If this parameter is a 
###             matrix of dimension (length(o_ts), #outputs), it is assumed that
###             R2 is a set of time-depedent weigting vectors, where each element 
###             is considered as the inverse of the weight for that moment in time.
###              (in other words, R2 = 1/W from the fir_mimo_rls3 function)
###             Eventually this routine will allow for full time-dependent
###             measurement noise covariance matrices.
###             (default = 1)
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
###             octave> reshape (Pbar_mtx(i,:), size(Pbar))';
###             octave>
###

### Revision History:
###
### 2005-07-20  First version with GNU/GPL license statement included
###             for public distribution (and possible modifications).  
###             - the "reshaped" storage matrices should be modified
###               in order to take advantage of the fact that Octave
###               now has good N-D matrix support.  Other than that,
###               there are no known bugs at this time.

function [theta, Pbar, lags, errs, theta_mtx, Pbar_mtx] = fir_mimo_kalman\
      (i_ts, o_ts, lagsize, theta_0, Pbar_0, R1, R2)

  ## We need to have at least three input parameters
  if (nargin < 3)
    error ("\nNot input enough parameters\n");
  endif


  ## Make sure the input and output time series matrices are composed of
  ## column vectors.  This is a bit of a kludge, since it assumes that the
  ## length (theta) << length (i_ts).
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
    theta_0 = 0;
    Pbar_0 = 1;
    R1 = 0;
    R2 = 1;
  elseif (nargin == 4)
    Pbar_0 = 1;
    R1 = 0;
    R2 = 1;
  elseif (nargin == 5)
    R1 = 0;
    R2 = 1;
  elseif (nargin == 6)
    R2 = 1;
  endif


  ## If lagsize is a scalar, assume that the FIR is length (2*lagsize+1)
  ## and centered at zero lag.  If lagsize is a 2-element vector, assume
  ## that the lag vector starts at lagsize(1), and ends at lagsize(2).
  ## This routine cannot currently handle lag vectors that are not
  ## simple indexes (i.e. sequential, increasing integers), so exit
  ## if lagsize is longer than 2 elements (i.e. don't trust the user
  ## to create his/her own lag vector properly).
  if (length(lagsize) > 2)
    error (["\nLag vector must be a scalar, or a 2-element vector\n", \
	    "that indicates the minimum and maximum lag indices.\n"]);
  endif

  if (length(lagsize) == 1)
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
  if ((rows(i_ts) - rows(o_ts)) != length(lags)-1 ||
      (length(lags) >= rows(i_ts)) )
    error (["\nInput and Output time series don\'t match,\n ", \
	    "or they are not long enough to determine a \n", \
	    "filter recursively.\n"]);
  endif


  ## A bunch of potentially useful static variables
  llags   = rows(lags);          ## number of lags
  ni_ts   = columns(i_ts);       ## number of different inputs
  no_ts   = columns(o_ts);       ## number of different outputs
  lo_ts   = rows(o_ts);          ## number of output data points
  li_ts   = rows(i_ts);          ## number of input data points


  ## Vectors and/or matrices should be initialized if possible
  yhat = zeros(no_ts,1);
  PHI = zeros(llags*ni_ts,1);


  ## These are matrices designed to store the error, state, and co-
  ## variance matrix at each time step.  Only initialize them if the
  ## user wants them returned because they suck up a LOT OF MEMORY!
  if (nargout > 3)
    errs = zeros(lo_ts,no_ts);
  endif
  if (nargout > 4)
    theta_mtx = zeros(lo_ts,llags*ni_ts*no_ts );
  endif
  if (nargout > 5)
    Pbar_mtx = zeros(lo_ts, (llags*ni_ts*no_ts)^2);
  endif


  ## If an initial value theta_0 was provided, 1) is it the correct size?
  ## 2) convert scalars to properly sized matrix of the same value; and
  ## 3) copy the value to the first element of theta.
  if (is_scalar (theta_0))
    theta = theta_0 * ones(columns(i_ts)*rows(lags), no_ts);
  elseif (is_vector(theta_0))
    if (length(theta_0) != (ni_ts*llags*no_ts))
      error (["\nFilter's initial value has wrong dimensions!\n"]);
    endif
    ## This is horrible, because we will simply reshape it again later,
    ## but for now parameterize theta according to form #1 in the Help section.
    theta = reshape(theta_0, [no_ts, ni_ts*llags])';
  elseif (is_matrix(theta_0))
    if (rows(theta_0) / ni_ts != llags || \
	columns(theta_0) != no_ts )
      error (["\nFilter's initial value has wrong dimensions!\n"]);
    endif
    theta = theta_0;
  else
    error (["\nFilter is wrong data-type, it must be a scalar or properly\n",\
	    "sized vector or matirx.\n"]);
  endif


  ## Initialize Pbar
  if (is_scalar(Pbar_0))
    Pbar = [eye(ni_ts*llags*no_ts) .* Pbar_0];
  elseif (is_vector(Pbar_0))
    if (length(Pbar_0) == (ni_ts*llags*no_ts))
      Pbar = diag(Pbar_0);
    elseif (length(Pbar_0) == (ni_ts*llags*no_ts)^2)
      Pbar = reshape(Pbar_0, (ni_ts*llags*no_ts), (ni_ts*llags*no_ts))';
    else
      error(["\nCovariance's initial value has wrong dimension!\n"]);
    endif
  elseif (is_matrix(Pbar_0))
    if (size(Pbar_0) != size(eye(ni_ts*llags*no_ts)) )
      error (["\nWrong size Covariance matrix!\n"]);
    endif
    Pbar = Pbar_0;
  else
    error (["\nCovariance matrix is wrong data-type, it must be a scalar or ",\
	    "properly sized vector or matirx.\n"]);
  endif


  ## Make sure the covariance matrix R1 is set up properly.
  tv_R1 = 0;
  if (is_scalar(R1))
    R1 = eye(llags*ni_ts*no_ts) * R1;
  elseif (is_vector(R1))
    if (length(R1) == (llags*ni_ts*no_ts))
      ## This pertinent for constant variances, but potentially multiple outputs
      R1 = diag(R1);
    elseif (size(R1) == [lo_ts, llags*ni_ts*no_ts])
      ## This is pertinent only if the desired state is a single scalar value
      ## with a time-dependent covariance...an extremely rare situation.
      ## Change nothing in R1, but set a flag so that the algorithm knows
      ## how to handle it later on.
      tv_R1 = 1;
    else
      error (["\n\"R1\" is not an appropriate dimension, read the help file.\n"]);
    endif
    
  elseif (is_matrix (R1))
    if (size(R1) == [lo_ts, (llags*ni_ts*no_ts)])
      ## this is pertinent if the diagonal process noise covariance is time-
      ## dependent..
      ## Change nothing in R2, but set a flag so that the algorithm knows how
      ## to handle it later on.
      tv_R1 = 1;
    elseif (size(R1) == [llags*ni_ts*no_ts,llags*ni_ts*no_ts])
      ## do nothing, R1 is OK as is
    else
      error (["\n\"R1\" is not an appropriate dimension, read the help file.\n"]);
    endif
  else
    error (["\n \"R1\" covariance matrix is wrong data-type, it must be a scalar or ",\
	    "properly sized vector or matirx.\n"]);    
  endif



  ## Make sure the covariance matrix R2 is set up properly.
  tv_R2 = 0;
  if (is_scalar(R2))
    R2 = eye(no_ts) * R2;
  elseif (is_vector(R2))
    if (length(R2) == no_ts)
      ## This pertinent for multiple outputs, with constant weighting
      R2 = diag(R2);
    elseif (size(R2) == size(o_ts))
      ## This is pertinent for single outputs, with time-varying weighting.
      ## Change nothing in R2, but set a flag so that the algorithm knows
      ## how to handle it later on.
      tv_R2 = 1;
    else
      error (["\"R2\" is not an appropriate dimension, read the help file.\n"]);
    endif

  elseif (is_matrix(R2)) 
    if (size(R2) == size(o_ts))
      ## this is pertinent for multiple outputs, with time-varying weighting.
      ## Change nothing in R2, but set a flag so that the algorithm knows how
      ## to handle it later on.
      tv_R2 = 1;
    elseif (size(R2) == size(eye(no_ts)))
      ## do nothing, R2 is OK as is
    else
      error (["\"R2\" is not an appropriate dimension, read the help file.\n"]);
    endif
  else
    error (["\n \"R2\" covariance matrix is wrong data-type, it must be a scalar or ",\
	    "properly sized vector or matirx.\n"]);    
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
      
      jj = ceil(j/ni_ts);
      PHI(j:j+ni_ts-1) = i_ts(abs(jj - (llags + 1) ) + (i-1), :)';

    endfor


    ## This converts PHI into a form compatible with the parameterization
    ## of theta described by form #2 in the Help section.
    PHI_tmp = kron(PHI,eye(no_ts));

    yhat = PHI_tmp' * theta_tmp; 
    err = o_ts(i,:)' - yhat;


    ## Check whether R1 is time-varying or not.
    if (tv_R1 == 1)
      R1_tmp = diag(R1(i,:));
    else
      R1_tmp = R1;
    endif


    ## Check whether R2 is time-varying or not.
    if (tv_R2 == 1)
      R2_tmp = diag(R2(i,:));
    else
      R2_tmp = R2;
    endif


    ## Return "rcond" to suppress warnings about singular matrices;
    ## maybe someday I'll actually check the "rcond".
    [denom,rcond] = inverse(PHI_tmp' * Pbar * PHI_tmp + R2_tmp );


    gamma = Pbar * PHI_tmp * denom;
    Pbar = Pbar + R1_tmp - gamma * PHI_tmp' * Pbar;
    theta_tmp = theta_tmp + gamma * err;

    #keyboard('Check Kalman Gain (gamma):> ')

    ## Store errs, state, and covariance matrices for each time step.
    ## Only store these if the user wants them returned, otherwise they
    ## suck up a WHOLE LOT OF MEMORY (especially the covariances)!
    if (nargout > 3)
      errs(i,:) = err';
    endif
    if (nargout > 4)
      theta_mtx(i,:) = theta_tmp';
    endif
    if (nargout > 5)
      Pbar_mtx(i,:) = reshape(Pbar', 1, rows(Pbar)*columns(Pbar));
    endif


    ## This probably slows things down slightly, but that's OK for now.
    printf ("\r %d of %d  -- %d",i,lo_ts,toc);
    fflush (stdout);

    
    if (nargout == 0)

      if(floor(i/200) == (i/200))
	keyboard ('Check P :>')
      endif
      
      #if (i == 1)
#	multiplot(ni_ts, no_ts)
      #endif

#      for k=1:no_ts
#	for l=1:ni_ts  
#	  #clearplot;
#	  grid ("on");
#	  mplot(lags, theta_tmp( k + (l-1)*no_ts  :ni_ts*no_ts:length(theta_tmp)),';;-' );
#	  replot
#	endfor
#      endfor

#      pause;
      
    endif

  endfor


  ## reshape the last state vector into the proper form
  theta = reshape(theta_tmp, size (theta'))';


  printf("\n");
  fflush(stdout);

endfunction
