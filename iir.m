## Calculates parameters for simple iir filter.
##
## Usage: [filta, filtb]=iir(alength, output, blength, input [, etype])
##
## where...
##
## alength : recursive coefficients
## input   : input time series
## blength : causal coefficients
## output  : measured output, driven by input time series
## etype   : error type (equation error, or output error)
##

function [filta, filtb] = iir (alength, output, blength, input, \
			       etype)

  if (nargin < 2 || nargin == 3 || nargin > 5)
    usage ("[filta, filtb]=iir(alength, output, blength, input [, etype])");
  endif

  if (nargin == 2)
    blength = [0];
    input = [0];
  endif

  if (nargin == 4)
    ## Don't do a repeated least squares to get the output error, unless
    ## explicitly requested...if we want this, the user must pass blength
    ## and input, even if they are only equal to zero!
    etype = 'ee';
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
  [auto_in, lag] = xcov (output, maxlag, "unbiased");


  ## Cross-covariance of input and output
  if (nargin >= 4)
    ## This is to kludge weirdness in xcov.m
    if (blength == 0)
      input = [];
    endif

    cross_inout = xcov(output, input, maxlag, "unbiased");

    ## This is to kludge weirdness in xcov.m
    if (blength == 0)
      input = 0;
    endif

  endif

  nrows = maxlag - 1;

  ## 
  ## There HAS to be a faster way to set up the lhs and rhs matrices!
  ## 

  ## Calculate the LHS matrix
  lhs_matrix=zeros(nrows,nrows);
  if blength ~= 0
    for m = 1:blength
      for l = 1:nrows
	lhs_matrix (l,m) = cross_inout ( (maxlag + 1) + (l-m) );
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

  filtb = theta (1:blength);
  filta = theta (blength+1:length (theta));


  if (nargin >= 5)

    ## If called with etype of "oe", or "output", or even "output error",
    ## filter input and output, and recursively call this function again
    ## to get the appropriate oe filter.
    if (strcmp (etype,"oe") || strcmp (etype,"output") || 
	strcmp (etype,"output error"))

      ## Maximum number of iterations
      for i=1:500
	
	## Prewhiten the signals with filta
	#outputf = filter ([filta], -[1], output);
	#inputf = filter ([filta], -[1], input);

	## This is the one that I like the best!
	## "filtfilt" can be used to remove the phase shift inherent to any
	## AR filter...this can cause problems with marginally stable,
	## and unstable filters over extended data sets (the problem
	## exists if the regular filter function is used, but it's
	## not as extreme).
	outputf = filter (1, [1;filta], (output));
	inputf = filter (1, [1;filta], (input));

	## If we wanted to truly "whiten" a periodic signal, this is
	## what we would do.  I guess what we do above is color the
	## signal in such a way as to generate white residuals.
	#outputf = filter ([1;filta],1,output);
	#inputf = filter ([1;filta],1,input);

	#keyboard

	[filtaf, filtbf] = iir (alength, outputf, blength, inputf, 'ee');

	#keyboard
	
	if (blength == 0)

	  disp ([i,max (abs ( [filtaf - filta]))]);
	  if (max (abs ( [filtaf - filta])) <= 1e-4 * max(abs(filtaf)))
	    filta = filtaf;
	    filtb = filtbf;
	    break;
	  endif

	elseif (alength == 0)

	  disp ([i,max (abs ( [filtbf - filtb]))]);
	  if (max (abs ( [filtbf - filtb])) <= 1e-4 * max(abs(filtbf)))
	    filta = filtaf;
	    filtb = filtbf;
	    break;
	  endif

	else

	  disp ([i,max (abs ( [filtaf - filta; filtbf - filtb]))]);
	  if (max (abs ( [filtaf - filta; filtbf - filtb])) <= 
	      1e-4 * max (abs ([filtaf;filtbf])) )
	    filta = filtaf;
	    filtb = filtbf;
	    break;
	  endif

	endif
	
	#keyboard

	filta = filtaf;
	filtb = filtbf;

	#hold off
	#plot (filta)
	#hold on
	#plot (outputf)

	## uncomment if we want this "recursive"
	#output = outputf;
	#input = inputf;
	
	
      endfor
      
    endif
    
  endif
  
endfunction