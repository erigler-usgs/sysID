###
### fir_cor_ls.m
###
### This routine will produce a simple Finite Impulse Response (FIR) linear 
### filter based on a given input and output time series.  This filter is the 
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
### [filt, lag, crosscov, autoin, autoout] = fir_cor_ls (its, ots[, lagsize])
###
### filt     - vector of optimal filter coefficients
### lag      - vector of corresponding lag values (units are whatever the
###            sampling interval is)
### crosscov - full cross-correlation vector for input/output time series
### autoin   - full auto-correlation vector for the input time series
### autoout  - full auto-correlation vector for the output time seriess
###

function [filt, lag, crosscov, autoin, autoout] = fir_cor_ls (its, ots, lagsize)

  if ( nargin < 2 || nargin > 3)
    usage (["\n\n[filter, crosscov, autocovin, autocovout] = ", \
	    "lpf (input_ts, output_ts [, lagsize])\n"]);
  endif

  ## Set a default lagsize of 10
  if nargin == 2
    lagsize = 10;
  endif

  ## Make sure the input and output time series vectors are column vectors
  if (rows(its)<columns(its))
    its = its';
  endif
  if (rows(ots) < columns(ots))
    ots = ots';
  endif

  ## Calculate auto-correlations
  autoin = xcov (its',"none");
  autoout = xcov (ots',"none");

  ## Calculate cross-correlations
  [crosscov,lag] = xcov (its', ots',"none");

  ## Set up the regression matrix
  ## (While it's not necessary with SVD, we'll generate a square matrix...
  ##  this means that we need the correlation vectors to be at least
  ##  (2*lagsize+1) long...that's something that should be checked for)
  covl = length (autoin); # length of correlation matrices
  auto_in_mtrx (1,:) = fliplr (autoin(round(covl/2)-(2*lagsize):
				      round(covl/2) ) );
  loops = 2*lagsize;
  for i=1:loops
    auto_in_mtrx (i+1,:) = fliplr (autoin(round(covl/2)-(2*lagsize)+i:
					  round(covl/2)+i ) );
  endfor


  ## Now use Single Value Decomposition to solve the system of equations
  [U,S,V] = svd (auto_in_mtrx);
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
  filt = V*(S\(U'*crosscov(round(covl/2)-lagsize:round(covl/2)+lagsize)' ) );
  lag = lag (round(covl/2)-lagsize:round(covl/2)+lagsize)';

  ## ...and that's all folks!

endfunction
