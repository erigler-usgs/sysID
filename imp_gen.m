## imp_gen.m
##
## THIS ROUTINE NEEDS TO BE BROKEN UP INTO SEPARATE FUNCTIONS.
## IT SHOULD REALLY BE NO MORE THAN A WRAPPER SCRIPT.
##
## Pre-whitens imput and output time-series, then generates a matrix of
## discrete impulse response functions and their lag indices
##

function [lagnfilt, lagnfilt_struct] = imp_gen (input, \
						output, \
						lagsize, \
						subds)
  
  if (nargin < 3 || nargin > 4)
    usage (["\n\n [lagnfilt, in_prw, out_prw, out_prw_prw] = ",\
	    "imp_gen (input, output, lagsize, [subds])\n"]);
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
  if (nargin == 4 )

    if (subds != floor(subds))
      error ("subds must be an integer\n");
    endif

    loops = subds;
    ll = length (input) / subds;
  else
    ## Default
    subds = 1;
  endif

  ## Make sure input and output vectors are same length
  if (rows (output) != length (input))
    error ("Input and Output dimenions don't match\n");
  endif
 
  ## Initialize temporary matrices...AR filters are of arbitrary
  ## length, and aren't known until arfit is run.
  in_prw = zeros (ll,1);
  out_prw = zeros (ll,columns(output));
  out_prw_prw = out_prw;

  ## Initialize lagnfilt matrix
  lagnfilt = zeros (2*lagsize+1,columns (output) + 1);
  lagnfilt (:,1) = [-lagsize:1:lagsize]';
  tmp_filt = lagnfilt;

  ## Initialize list of matrices
  ## (lists aren't saved in octave or matlab binary files...
  ##  strangely, structures are)
  ##lagnfilt_list = list ();
  
  for i = 0:loops-1
    
    #if i != 4 
    #  continue
    #endif

    ## Prewhiten I/O with optimal AR filter (assumes Input causes
    ## Output)
    [w,prw_io] = arfit (input (i*ll+1 : (i+1)*ll) - 
			mean (input (i*ll+1 : (i+1)*ll)),
			1, 100, 'zero');

    in_prw = filter ([1;-prw_io'], 1, input (i*ll+1 : (i+1)*ll) -
		     mean (input (i*ll+1 : (i+1)*ll)) );
    
    ## Remove any remaining autocorrelation from the output time series
    ## since this is either diurnal variability or some sort of
    ## magnetospheric ringing...neither of which should impact the
    ## impulse response function.
    for j = 1:columns (output)
      
      #if j != 40
      #  continue
      #endif

      out_prw(:,j) = filter ([1;-prw_io'], 1, output (i*ll+1 : (i+1)*ll, j) -
			     mean (output (i*ll+1 : (i+1)*ll, j)) );
 
      #[w, prw_oo] = arfit (out_prw (:,j) - mean (out_prw (:,j)), 1, 100, \
	#		   'zero');

      #out_prw_prw (:,j) = filter ([1;-prw_oo'], 1, out_prw (:,j) - \
	#			  mean (out_prw (:,j)) );


      ## Calculate filter, and update output matrices and "list" 
      ## of intermediate filters.
      tmp_filt (:,j+1) = fir_cor_ls (in_prw, out_prw_prw (:,j), lagsize);
      #tmp_filt (:,j+1) = lpf (in_prw, out_prw (:,j), lagsize);
      #tmp_filt (:,j+1) = lpf (input, output (:,j), lagsize);
      lagnfilt (:,j+1) = lagnfilt (:,j+1) + tmp_filt (:,j+1);

      printf ("\r %d %d",i,j);
      fflush (stdout);

      #keyboard

    endfor

    ## append the current filter matrix to the "list"
    ## (lists won't save to either octave or matlab binary files...
    ##  strangely structures will save to matlab files)
    ##lagnfilt_list = append (lagnfilt_list, tmp_filt);

    ## save intermediate filter matrices to a structure
    ## (structures require strings to identify their subelements,
    ##  and I don't want to deal with multi-character strings...
    ##  so we won't allow more than 52 filter matrices, allowing
    ##  for subelements named [A..Z,a..z])
    if (nargout > 1)
      if (loops > 52)
	error ("Too many intermediate matrices requested\n");
      endif
      eval (sprintf ("lagnfilt_struct.%s = tmp_filt;",setstr(i+65)));
    endif
      
  endfor

  lagnfilt (:,2:columns(lagnfilt)) = lagnfilt (:,2:columns(lagnfilt)) \
      ./ subds;


endfunction
