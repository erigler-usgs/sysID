## Calculates optimal OE IIR filter through repeated least squares
## techniques, as described in Nelles, 2001
##
## Usage: [theta]=repls(finput,foutput,theta)
##
## where...
##
## finput   - filtered input
## foutput  - filtered output
## theta    - 

function [theta2,length_A,lengthB]=repls(input_ts,output_ts,theta,length_A,length_B,tol)

  if nargin < 5
    printf ("\nNot enough arguments...quitting!\n\n");
    theta2 = -1;
    return
  elseif nargin < 6
    tol = 1e-6;
  endif

  ## Copy feedback filter to temporary filter, and put it in the
  ## denominator 
  filt = -1*theta(1:length_A);
  
  ## Make certain input and output are column vectors
  [nr,nc] = size (input_ts);
  if nr ~= 1 & nc ~=1
    error ("Only works with 1-D input vectors for now\n");
  endif
  if nr < nc, input_ts = input_ts'; endif

  [nr,nc] = size(output_ts);
  if nr ~= 1 & nc ~=1
    error ("Only works with 1-D output vectors for now\n");
  endif
  if nr < nc, output_ts = output_ts'; endif


  ## Make certain filter is a column vector
  [nr,nc] = size(filter);
  if nr ~= 1 & nc ~= 1 
    error ("Only works with theta as a vector for now\n");
    theta2=-1;
    return;
  endif
  if nr < nc, filt = filt'; endif



  ## Initialize the output vector
  finput=zeros(1,length(input_ts));
  foutput=zeros(1,length(input_ts));

  ## Begin iteration loop.  10 is an arbitrary number...if this doesn't
  ##converge in this time, give up.
  for loop=1:10
    
    printf("Starting iteration #%d \n",loop);

    ## Starting at t=max([length_A,length_B]), increment and update the
    ## output vector.  The "+1" is because we are assuming that the
    ## filters are "strictly proper", meaning that input and/or feedback
    ## at  time "t" cannot instantaneously affect the output.
    for t=max([length_A,length_B]):length(input_ts);
      
      if length_B == 1
	## a simple rot90 does the same as a transpose(rot90(rot90()))
	phi = (rot90(output_ts(t-length_A+1:t) ) );
      elseif length_A == 1
	## a simple rot90 does the same as a transpose(rot90(rot90()))
	phi = (rot90(input_ts(t-length_B+1:t) ) ) ;
      else
	phi = [(rot90((output_ts(t-length_A+1:t)) ) );(rot90( (input_ts(t-length_B+1:t)) ) )];
      endif
      
      
      ## filtering by the feedback filter only...I don't really understand
      ## this, but it's what's written in Nelles, 2001
      foutput(t) = transpose(filter(1:length_B)) * phi(1:length_B);
      finput(t) = transpose(filter(1:length_A)) * phi(length_B+1:length(phi));
      
    endfor # end creation of filtered input and output
    
    ## Recalculate filter
    [theta2,length_A,length_B]=iir3(finput,foutput,length_A,length_B);
    
    ## Check to see if delta-theta is less than convergence tolerance
    ## (use maximum delta between thetas as criterion)

    keyboard
    
    if max(abs(theta-theta2)) <= tol
      return;
    endif

    ## If tolerance wasn't met, update filter and copy theta2 to theta
    filter = -1 * theta2(1:length_A);
    theta=theta2;
    

  endfor # end iteration loop
  
  ## Shouldn't get here if we got a solution
  printf("\nCould not converge in %d iterations...quitting!\n\n",loop);
  
endfunction