## Calculates parameters for simple iir filter.
##
## Usage: [theta,alength,blength] = iir (input, output, alength, blength)
##
## where...
##
## input   : input time series
## output  : measured output, driven by input time series
## flength : for now A and B lengths are the same and equal to flength
## theta   : resulting filter with both causal and recursive
##           coefficients
## alength : recursive coefficients
## blength : causal coefficients
##

function [filta, filtb]=iir(alength, output, blength, input)

  if (nargin < 2 || nargin == 3 || nargin > 4)
    usage ("[filta, filtb]=iir(alength, output, blength, input)");
  endif

  if (nargin == 2)
    blength = 0;
    input = 0;
  endif


  [nr,nc] = size (input);
  if nr ~= 1 & nc ~=1
    error ("Only works with 1-D input vectors for now\n");
  endif
  if nr < nc, input=transpose(input); endif
  
  [nr,nc] = size(output);
  if nr ~= 1 & nc ~=1
    error ("Only works with 1-D output vectors for now\n");
  endif
  if nr < nc, output=transpose(output); endif
  

  ##maxlag = max([alength,blength]);
  maxlag = alength + blength + 1;

  ## Auto-covariance of input
  [auto_in, lag] = xcov (output, maxlag);

  ## Cross-covariance of input and output
  if (nargin == 4)
    cross_inout = xcov(output, input, maxlag);
  endif

  nrows = maxlag - 1;


  ## Calculate the LHS matrix
  ## (include zerolag for input)
  lhs_matrix=zeros(nrows,nrows);
  if blength ~= 0
    for m = 1:blength
      for l = 1:nrows
	lhs_matrix (l,m) = cross_inout ( (maxlag + 1) + (l-m) + 1);	
      endfor
    endfor
  endif


  if alength ~= 0
    for n=1:alength
      for l=1:nrows
	lhs_matrix (l,n+blength) = -1 * auto_in ( (maxlag + 1) + (l-n) );
      endfor
    endfor
  endif


  ## Calculate the RHS matrix (should be cross correlations one timestep
  ## in the future
  rhs_matrix = auto_in (maxlag + 2: maxlag + 2 + nrows-1);


  ## Use SVD and back substitution to solve
  [U,S,V] = svd (lhs_matrix);
  ## take out very small values of S
  S = diag (S);
  small = find (S/max(S) <= 1e-9);

  if small
    S(small) = 0;
  endif
  S = diag (S);

  theta = V*(S\(U'*rhs_matrix));

  ## Calculate the filter vector from overdetermined system via LS
  ## method
  ##
  ##theta=inv(lhs_matrix' * lhs_matrix) * lhs_matrix' * rhs_matrix ;


  filtb = theta (1:blength);
  filta = theta (blength+1:length (theta));

endfunction