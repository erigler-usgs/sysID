###
### rpe_iss.m
###
### This function applies an recursive prediction error identification
### algorithm to the innovations form state space model (ISS, see Ljung
### and Soderstram, 1983, Chapter 3, Eqs 3.145) in order to estimate both 
### the optimal state as well as the linear parameters that comprise
### properly formed polynomial matrices that can be used as follows:
### 
###
###   x(t+1) = A*x(t) + B*u(t) + K*err(t)
###     y(t) = C*x(t) + D*u(t) + err(t)
###
### ...or, as an optimal 1-step predictor...
###
###   x(t+1) = (A-KC)*x(t) + (B-KD)*u(t) + K*y(t)
###   _
###   y(t+1) = C*x(t+1) + D*u(t+1)
###
###
### Usage:
###
### {[Y, ISS, X, Qs, Retries, ISS_Mtrx, X_Mtrx, Qs_Mtrx]} = \
###                      rpe_iss (uobs, yobs {, ISS0, X0, Qs0, maxRetries} )
###
### ...where "{}" indicates optional parameters...
###
### INPUTS:
###
### uobs       - Input time-series matrix.  Columns are different inputs
###              and rows are sequential input observations.
###
### yobs       - Output time-series matrix.  Columns are different outputs
###              and rows are sequential output observations.
###
###              Note0:   If only the first two input parameters are provided, 
###                     it is assumed that the user wants to identify a system
###                     of order equal to the number of outputs that fits the
###                     canonical observable form, that ISS0 consists of zeros
###                     for both the polynomial matrices and shadow matrices
###                     (implying that a time stationary system is to be
###                      identified)
###
### ISS0       - Structure containing initial matrices for an innovations
###              form state-space model.  Only a basic check for dimensional
###              consistency will be performed.  The following example is 
###              for a 4th order, 2 input, 3 output system in the canonical 
###              observable form necessary to allow conversion to traditional
###              multi-parameter I/O filters (thus the zeros and ones):
###
###
###     ISS0.A = [ 0   1   0   0 ] ISS0.B = [b01 b02] ISS0.K = [k01 k02 k03]
###              [a01 a02 a03 a04]          [b03 b04]          [k04 k05 k06]
###              [a05 a06 a07 a08]          [b05 b06]          [k07 k08 k09]
###              [a09 a10 a11 a12]          [b07 b08]          [k10 k11 k12]
###
###     ISS0.C = [ 1   0   0   0 ] ISS0.D = [d01 d02]
###              [ 0   1   0   0 ]          [d03 d04]
###              [ 0   0   1   0 ]          [d05 d06]
###
###
###               Each polynomial matrix may have a corresponding "shadow
###              matrix", which is only relevant when this function is being
###              used for system identification.  Each element of the shadow
###              matrix is either a scalar describing the variance of the
###              parameter, and will be ultimately placed in a diagonal
###              parameter process noise covariance matrix that corresponds
###              to the vector of adjustable parameters that will augment
###              state vector, OR it can be a NaN, indicating that that 
###              parameter is NOT adjustable at all. For example if time-
###              stationary versions of the above matrices needed to
###              be determined, the shadow matrices would be:
###
###
###     ISS0.Av = [NaN NaN NaN NaN] ISS0.Bv = [ 0   0 ] ISS0.Kv = [ 0   0   0 ]
###               [ 0   0   0   0 ]           [ 0   0 ]           [ 0   0   0 ]
###               [ 0   0   0   0 ]           [ 0   0 ]           [ 0   0   0 ]
###               [ 0   0   0   0 ]           [ 0   0 ]           [ 0   0   0 ]
###
###     ISS0.Cv = [NaN NaN NaN NaN] ISS0.Dv = [ 0   0 ]
###               [NaN NaN NaN NaN]           [ 0   0 ]
###               [NaN NaN NaN NaN]           [ 0   0 ]
###
###
###              *   Eventually, if yobs(t) = NaN (i.e. no available output),
###                we might have the parameters sort of decay exponentially
###                to the average time-stationary values.  This constitutes
###                a kind of Gauss-Markov process, but where the background
###                state is still determined by a set of dynamical difference
###                equations, not some static or quasi-static climatology.
###                  For now, we simply do not update the parameters when there
###                are no observations from which to calculate residuals, and
###                thus determine the appropriate update.
###                -EJR (10/21/04)
###
###         
### X0         - Initial state vector of the system.  If this vector is not
###              provided, it is set to an appropriate length zero vector.
###
### Qs0        - Structure containing initial derivative matrices W and Psi, 
###              as well as the innovations covariance (G, for big gamma), and
###              the covariance matrix P, which is equal to Gss*inv(R).
###              If any or all of these are missing, default initializers
###              will be generated which should eventually converge to optimal  
###              values.  In addition, a shadow matrix for G can be specified (Gss)
###              specifying the step size (little gamma) of the recursive
###              algorithm used to determine the optimal G independently of
###              the RPE algorithm (P, W, and Psi are determined in the algorithm
###              proper).  For the sample 4th order, 2 input, 3 output system
###              described thus far, the matrices might be:
###
###             Qs0.G = [g01 g02 g03 ]       Qs0.Gss = [.1   0   0 ]
###                     [g04 g05 g06 ]                 [ 0  .1   0 ]
###                     [g07 g08 g09 ]                 [ 0   0  .1 ]
###
###            Qs0.W =  [ 0   0 .(# adjustable coefficients). 0 ]
###                     [ 0   0 ............................. 0 ]
###                     [ 0   0 ............................. 0 ]
###                     [ 0   0 ............................. 0 ]
###
###           Qs0.Psi = [ 0   0   0 ]
###                     [ :       : ]
###                     [ #       : ]
###                     [ c       : ]
###                     [ o       : ]
###                     [ e       : ]
###                     [ f       : ]
###                     [ f       : ]
###                     [ i       : ]
###                     [ c       : ]
###                     [ i       : ]
###                     [ e       : ]
###                     [ n       : ]
###                     [ t       : ]
###                     [ s       : ]
###                     [ :       : ]
###                     [ 0   0   0 ]
###
###             Qs0.P = [ 0   0 .(# adjustable coefficients). 0 ]
###                     [ 0                                   : ]
###                     [ :                                   : ]
###                     [ :                                   : ]
###                     [ :                                   : ]
###                     [ :                                   : ]
###                     [ :                                   : ]
###                     [ :              ( Square )           : ]
###                     [ :                                   : ]
###                     [ :              ( Matrix )           : ]
###                     [ :                                   : ]
###                     [ :                                   : ]
###                     [ :                                   : ]
###                     [ :                                   : ]
###                     [ :                                   : ]
###                     [ 0   0 ............................. 0 ]
###
###
###              Note1:   G essentially scales the correction of the model
###                     coefficients at each time step.  If G is allowed to 
###                     vary it is entirely possible that the identified
###                     coefficients will not converge to a time-stationary
###                     solution, even if the values in the shadow matrices
###                     of the ISS structure are equal to zero.
###
###              Note2:   Gss should be between 0 (implies G is constant) and 1 
###                     (implies that G is the latest instantaneous error^2,
###                     and previous G's are disregarded completely).  The
###                     recursive algorithm is quite simple:
###                                                       T
###                     G(t) = G(t-1) + Gss .* (err(t)*err (t) - G(t-1))
###
###                       Gss may also be a string describing a function of time-
###                     step (t).  For example, if Gss='1/t', then the value 1/t
###                     will be used at each time step.  If Gss is not defined
###                     a simple hard-coded algorithm is applied that forces
###                     GSS-->0 as t-->infinity.
###
###              Note3:   Shadow matrices Pss, Wss, and Psiss may actually be
###                     set but only have meaning if Qs_Mtrx is requested as an
###                     output.  If this is the case, each P, W, and Psi
###                     will be accumulated and returned to the user.  THIS
###                     MIGHT EASILY USE UP ALL AVAILABLE MEMORY, BE CAREFUL!
###
### maxStabCheck- The function is designed to check whether or not the matrix
###              [A(t)-K(t)C(t)] is exponentially stable (are its eigen values
###               all within the unit circle?).  If it is not, the correction
###              just applied to the adjustable parameters is halved and applied
###              again, and the eigen values are checked once more.  This process
###              is repeated until a stable matrix is produced, or the number of
###              stability checks reaches maxStabCheck.  Then the correction is
###              dropped and the coefficients from the previous time step are used.
###
###              Note4:   Obviously this requires that ISS0 start out with matrices
###                     that fall within the subset of the exponentially stable
###                     [A(t)-K(t)C(t)].
###
###              Note5:   If maxStabCheck=0, no stability check is performed, and 
###                     the update is unconditionally applied.  This might be
###                     dangerous, but it makes sense when one realizes that
###                     it's the stability check that sucks up the most CPU
###                     time; halving the correction factor and readjusting
###                     the coefficients is trivial by comparison.  If one is
###                     very concerned about the stability of [A(t)-K(t)C(t)],
###                     one should always set maxRetries>0.
###
### OUTPUTS:
###
### Y          - Time series of predicted output.
###
### ISS        - Structure containing last set of polynomial matrices.  This 
###              can be passed as-is back into this function as ISS0 along with
###              sequential input and output observations, picking up where the
###              last run left off.
###
### X            The last state vector.  This can be passed back into the function
###              to pick up state estimation exactly where it left off.
###
### Qs         - Structure containing las set of derivative matrices.  This can
###              be passed as-is back into this function as Qs0 to pick up where
###              the last run left off.
###
### StabCheck  - Each row of this vector corresponds to an update of the ISS
###              adjustable parameters.  The number of checks corresponds to
###              how many times it was necessary to halve the correction factor
###              before the matrix [A(t)-K(t)C(t)] was exponentially stable.
###              If the eigen values could not be pushed back inside the unit
###              circle after maxStabCheck attempts, a value of -maxStabCheck
###              is placed here.
###
###              Note6:   If a negative value is returned, this does not necessarily
###                      mean that the filter was unstable, since the stability
###                      check is only performed once for each loop.  It is possible
###                      that the very last correction managed to push roots back
###                      inside the unit circle.  It is up to the user to check the
###                      eigen values of the polynomial matrices corresponding to
###                      time steps with negative StabCheck values.
###
###
###
### Iss_Mtrx   - This "column" vector holds a time series of structures containing the 
###              state space polynomial matrices, and is only accumulated if asked
###              for since it can potentially use a tremendous amount of memory.
###              It's last element should correspond to the single structure ISS.
###
### X_Mtrx     - This "column" vector holds a time series of structures of state vectors.
###              It is probably not necessary to represent it this way (structures with
###              single state vectors as elements) but is done so for consistency with 
###              the other *_Mtrx outputs.
###
### Qs_Mtrx    - This "column" vector holds a time series of structures containing the
###              covariance matrices, and is only accumulated if asked for since it
###              can potentially use a tremendous amount of memory.  In fact, these
###              matrices can get SO big, the structures will only contain P, W, and Psi
###              IF otherwise useless corresponding shadow matrices Pss, Wss, and Psiss
###              were included in Qs0.
###              
###              Note7:    Really, I mean it, be careful with Pss, Wss, and Psiss!  For
###                      the sample 4th order, 2 input, 3 output system described above,
###                      each time step will result in an ISS_Mtrx element with
###                           (4x4)+(4x2)+(4x3)+(3x4)+(3x2) = 42 elements.  
###                      Each time step would result in a Qs_mtrx element with 
###                           (3x3)+(4x38)+(38x3)+(38x38) = 1719 elements!!!
###                      (and remember, we weren't estimating every element of the
###                       state space polynomial matrices, just those required for
###                       canonical observable form)
###

function [Y, ISS, X, Qs, StabCheck, ISS_Mtrx, X_Mtrx, Qs_Mtrx,Gss_Mtrx] = \
        rpe_iss (uobs, yobs, ISS0, X0, Qs0, maxStabCheck)


        if (nargin < 2)
            error("\nNot input enough parameters\n");
        elseif (nargin > 6)
            error("\nToo many input parameters\n");
        elseif (nargin == 2)
                
        ##
        ##  Assume default ISS0, X0, Qs0, and maxStabCheck values
        ##
            
            nu = columns(uobs);         # number of different inputs
            ny = columns(yobs);         # number of different outputs
            n = ny;                     # system order (default)
            
            ISS0.A = zeros(n,n);
            #for i=1:n for j=1:n if j==i+1 ISS0.A(i,j)=1; endif; endfor; endfor;
            ISS0.Av = zeros(n,n);       # no process noise for A
            #for i=1:n for j=1:n if i~=n ISS0.Av(i,j)=NaN; endif; endfor; endfor;
            ISS0.B = zeros(n,nu);
            ISS0.Bv = zeros(n,nu);      # no process noise for B
            ISS0.K = zeros(n,ny);
            ISS0.Kv = zeros(n,ny);      # no process noise for K
            ISS0.C = eye(ny,n);
            ISS0.Cv = nan*eye(ny,n);    # fix parameters for C
            ISS0.D = zeros(ny,nu);
            ISS0.Dv = nan*zeros(ny,nu); # fix parameters for D            
            
        endif
            
        if (~ (isfield(ISS0,'A') & isfield(ISS0,'B') & isfield(ISS0,'K') & isfield(ISS0,'C') & isfield(ISS0,'D') ) )
            error("\nA, B, K, C, and D matrices must reside in the ISS0 structure\n");
        else
            
            nu = columns(uobs);         # number of different inputs
            ny = columns(yobs);         # number of different outputs

            ##
            ## Check the dimensions of ISS0.*
            ##
            if (issquare(ISS0.A)) 
                n=rows(ISS0.A); 
            else
                error("\nThe A matrix must be square, with dimension equal to the order of the system being modeled\n");
            endif
            
            if (rows(ISS0.B) ~= rows(ISS0.A))
                error("\nThe number of rows in A and B matrices must match\n");
            endif
            
            if (columns(ISS0.B) ~= nu)
                error("\nThe number of columns in B matrix must match the number of inputs\n");
            endif
            
            if (n ~= rows(ISS0.K))
                error("\nThe number of rows in K matrix must match the order of the square A matrix\n");
            endif
            
            if (ny ~= columns(ISS0.K))
                error("\nThe number of columns in K matrix must match the number of outputs\n");
            endif
            
            if (n ~= columns(ISS0.C))
                error("\nThe number of columns in C matrix must match the order of the square A matrix\n");
            endif
            
            if (ny ~= rows(ISS0.C))
                error("\nThe number of rows in C matrix must match the number of outputs\n");
            endif
            
            if (size(ISS0.D) ~= [ny,nu])
                error(["\nThe D matrix must have number of rows equal to number of outputs",\
                       "\nand number of columns equal to number of inputs\n"]);
            endif
            
            
            ##
            ## Check for "shadow matrices" in ISS0.  If they exist make sure they have
            ## the right dimension; if they don't exist, assume that this function was
            ## not called to identify optimal parameters, but simply to filter some data,
            ## and set all shadow elements to NaNs
            ##
            if (isfield(ISS0,'Av') )
                if any(size(ISS0.A) ~= size(ISS0.Av) )
                    error("\nThe A shadow matrix must have same dimension as A matrix\n");
                endif
            else
                ISS0.Av = nans*zeros(size(ISS0.A));
            endif
            
            if (isfield(ISS0,'Bv') )
                if any(size(ISS0.B) ~= size(ISS0.Bv))
                    error("\nThe B shadow matrix must have same dimension as B matrix\n");
                endif
            else
                ISS0.Bv = nans*zeros(size(ISS0.B));
            endif
            
            if (isfield(ISS0,'Kv') )
                if any(size(ISS0.K) ~= size(ISS0.Kv))
                    error("\nThe K shadow matrix must have same dimension as K matrix\n");
                endif
            else
                ISS0.Kv = nans*zeros(size(ISS0.K));
            endif
            
            if (isfield(ISS0,'Cv') )
                if any(size(ISS0.C) ~= size(ISS0.Cv))
                    error("\nThe C shadow matrix must have same dimension as C matrix\n");
                endif
            else
                ISS0.Cv = nans*zeros(size(ISS0.C));
            endif
            
            if (isfield(ISS0,'Dv') )
                if any(size(ISS0.D) ~= size(ISS0.Dv))
                    error("\nThe D shadow matrix must have same dimension as D matrix\n");
                endif
            else
                ISS0.Dv = nans*zeros(size(ISS0.D));
            endif


        endif
        
        
        ##
        ## Now, it's necessary to determine the number of adjustable parameters
        ## (there must be a better way, but I don't know it)
        nv=0; 
        for [val,key] = ISS0 
            if (key(end)=='v')
                nv = nv+numel(find(~ isnan(val)));
            endif; 
        endfor
        

        ##
        ## Simply make X0 a vector of zeros if it wasn't passed in
        ##
        if (~ exist('X0'))
            X0=zeros(n,1);
        elseif (~ ismatrix(X0))
            X0 = zeros(n,1);
        endif
        if (size(X0) ~= [n,1])
            error("\nX0 must be a column vector with length that matches the order of square A matrix\n");
        endif
        
        
        ##
        ## Check existence of Qs0 structure
        ##
        if (~ (exist('Qs0')) )
            Qs0=struct;
        elseif (~ isstruct(Qs0))
            Qs0=struct;
        endif
        
        ##
        ## Check existence and dimensions of initial innovations covariance matrix
        ##
        if (~ (isfield(Qs0,'G')) )
            Qs0.G = zeros(ny,ny);
        elseif (size(Qs0.G)~=[ny,ny])
            error("\nThe innovations covariance matrix G must must be square, with dimension equal to the number of outputs\n");
        endif
        

        
        ##
        ## Check existence and dimensions of initial derivative matrices
        ##
        if (~ (isfield(Qs0,'W')) )
            Qs0.W = zeros(n,nv);
        elseif (size(Qs0.W)~=[n,nv])
            error(["\nThe derivative matrix W should have a number of rows equal to the system order",\
                   "\nand a number of columns equal to the number of adjustable parameters\n"]);
        endif
        
        if (~ (isfield(Qs0,'Psi')) )
            Qs0.Psi = zeros(nv,ny);
        elseif (size(Qs0.Psi)~=[nv,ny])
            error(["\nThe derivative matrix Psi should have a number of rows equal to the number of",\
                   "\nadjustable parameters, and a number of columns equal to the number of outputs\n"]);
        endif
        
        
        ##
        ## Check existence and dimensions of initial covariance matrix P
        ##
        if (~ (isfield(Qs0,'P')) )
            # This is a typical choice for initial guess, since a matrix of
            # zeros would result in non-adjustable parameters
            Qs0.P = 100 * max(nanstd(yobs).^2) * eye(nv,nv);
        elseif (size(Qs0.P)~=[nv,nv])
            error("\nThe covariance matrix P must be square, with dimension equal to the number of adjustable parameters\n");
        endif
        

        ##
        ## If maxStabCheck is not set, assume that no stability monitoring will be performed
        ##
        if (~ exist('maxStabCheck'))
            maxStabCheck=0;
        endif


        ##
        ## These gain sequence adjusters are hard coded because they
        ## are only really relevant in the transient phase of identificaton,
        ## and only if Qs0.Gss is not set by the user.  They force the
        ## algorithm to start with a forgetting factor of 0.9 and exponentially
        ## approach 1 at a rate defined by lam0.
        ##
        lam0=0.95;
        lam=0.9;
        Gss=.1;
        


        ##
        ## OK, now that all the inputs have been checked and/or set to defaul values
        ## extract theta vector and parameter variances from ISS0
        ##
        Av_idx = find(~ isnan(ISS0.Av));
        Bv_idx = find(~ isnan(ISS0.Bv));
        Kv_idx = find(~ isnan(ISS0.Kv));
        Cv_idx = find(~ isnan(ISS0.Cv));
        Dv_idx = find(~ isnan(ISS0.Dv));

        theta = [ISS0.A(Av_idx);
                 ISS0.B(Bv_idx);
                 ISS0.K(Kv_idx);
                 ISS0.C(Cv_idx);
                 ISS0.D(Dv_idx)];

        ##
        ## Set up a process noise covariance matrix using non-NaN values from ISS0 matrices
        ##
        Qv =    [ISS0.Av(Av_idx);
                 ISS0.Bv(Bv_idx);
                 ISS0.Kv(Kv_idx);
                 ISS0.Cv(Cv_idx);  # need to fix it so that size(*_idx) always returns [1 0]
                 ISS0.Dv(Dv_idx)]; # instead of [0 1] for zero-matrices
        Qv = diag(Qv);             # perhaps someday we can modify this to allow non-diagonal
                                   # Qv matrices, but this should suffice for now

                 
        
        ##
        ##  Convert fortran-style indices to column wise indexes, outside of
        ## loop for speed, for later use in determining derivative matrices
        ##
        Av_cidx = mod(Av_idx-1,n)+1;
        Bv_cidx = mod(Bv_idx-1,n)+1;
        Kv_cidx = mod(Kv_idx-1,n)+1;
        
        Cv_cidx = mod(Cv_idx-1,ny)+1;
        Dv_cidx = mod(Dv_idx-1,ny)+1;
        
        
        ##
        ## Preallocate matrices for derivatives
        ##
        Mt_A = zeros(n,length(Av_cidx));
        Mt_B = zeros(n,length(Bv_cidx));
        Mt_K = zeros(n,length(Kv_cidx));
        Mt = zeros(n,nv);
        
        Dt_C=zeros(ny,length(Cv_cidx));
        Dt_D=zeros(ny,length(Dv_cidx));
        Dt = zeros(ny,nv);
        
        
        ##
        ## Initialize matrix for predicted output
        ##
        Y = zeros(size(yobs));
        
        
        ##
        ## Initialize other output matrices if they were called for
        ##
        if (nargout > 9)
            error("\nToo many output arguments\n");
        endif
        
        if (nargout > 4)
            StabCheck = zeros(rows(yobs),1);
        endif

        if (nargout > 5)
            ISS_Mtrx = struct('A',ISS0.A,'B',ISS0.B,'K',ISS0.K,'C',ISS0.C,'D',ISS0.D);
            ISS_Mtrx(1:rows(yobs),1) = ISS_Mtrx; 
        endif
        
        if (nargout > 6)
            X_Mtrx = struct('X',X0);
            X_Mtrx(1:rows(yobs),1) = X_Mtrx;
        endif
        
        if (nargout > 7)
            Qs_Mtrx.G = Qs0.G;
            if (isfield(Qs0,'Wss') )
                Qs_Mtrx.W = Qs0.W;
            endif
            if (isfield(Qs0,'Psiss') )
                Qs_Mtrx.Psi = Qs0.Psi;
            endif
            if (isfield(Qs0,'Pss') )
                Qs_Mtrx.P = Qs0.P;
            endif
            Qs_Mtrx(1:rows(yobs),1) = Qs_Mtrx;
        endif
        
        if (nargout > 8)
            Gss_mtrx = zeros(rows(yobs),1);
        endif
        
        
        ##
        ## Start loop
        ##
        
        keyboard ('Before Loop > ');
        
        for t=1:rows(yobs)
                        
            obs_avail = ~ any(isnan(yobs(t,:)) );
            ins_avail = ~ any(isnan(uobs(t,:)) );
            
            
            ##
            ## Determine observables and calculate error vector
            ##
            if (ins_avail)
                yhat = ISS0.C * X0 + ISS0.D * uobs(t,:)';
            else
                yhat = ISS0.C * X0;
            endif
            
            if (obs_avail)
                err = yobs(t,:)' - yhat;
            endif
            
            
            ##
            ## Update innovations covariance (determine gain (GSS)
            ##  and forgetting factor (lam) first)
            ##
            if (obs_avail)
                
                if (isfield(Qs0,'Gss'))
                    if isstr(Qs0.Gss)
                        Gss=eval(Qs0.Gss);
                        lam=1-Gss;
                    else
                        Gss=Qs0.Gss;
                        lam=1-Gss;
                    endif
                else
                    lam=lam0*lam+(1-lam0);
                    Gss = 1/(1+(lam/Gss));
                endif
                Qs0.G = Qs0.G + Gss * (err*err' - Qs0.G);
                Gss_Mtrx(t) = Gss;
            endif


            ##
            ## Determine the derivative Mt (d/dtheta(Ax+Bu+Ke))
            ## (If this gets really slow, it might be a good idea to
            ##  extract those parts that can be calulated outside of this
            ##  loop; this is really only Mt_B and Dt_D)
            ##
            if (obs_avail)
                for i=1:length(Av_cidx)
                    Mt_A(Av_cidx(i),i) = X0(ceil(Av_idx/n)(i));
                endfor
                if ismatrix(Mt_A)  # If zero-length matrix, don't do anything
                    Mt(1:numel(Mt_A)) = Mt_A(:);
                endif
                
                if (ins_avail)
                    for i=1:length(Bv_cidx)
                        Mt_B(Bv_cidx(i),i) = uobs(t,ceil(Bv_idx/n)(i));
                    endfor
                else
                    Mt_B(:) = 0; # Not sure if this is correct approach or not
                endif
                
                if ismatrix(Mt_B)  # If zero-length matrix, don't do anything
                    Mt(numel(Mt_A)+1:numel(Mt_A)+numel(Mt_B)) = Mt_B(:);
                endif

                for i=1:length(Kv_cidx)
                    Mt_K(Kv_cidx(i),i) = err(ceil(Kv_idx/n)(i));
                endfor
                if ismatrix(Mt_K)  # If zero-length matrix, don't do anything
                    Mt(numel(Mt_A)+numel(Mt_B)+1:numel(Mt_A)+numel(Mt_B)+numel(Mt_K)) = Mt_K(:);
                endif
            endif

            
            ##
            ## Update L gain matrix
            ##
            if (obs_avail)
                S = Qs0.Psi' * Qs0.P * Qs0.Psi + lam * Qs0.G;
                L = Qs0.P * Qs0.Psi * inv(S);
            endif

            ##
            ## Calculate theta(t) and update ISS0.*
            ##
            if (obs_avail)
                theta_old = theta;
                theta = theta + L * err;

                if ismatrix(Av_idx)
                    ISS0.A(Av_idx) = theta(find(Av_idx));
                endif
                if ismatrix(Bv_idx)
                    ISS0.B(Bv_idx) = theta(numel(Av_idx)+find(Bv_idx));
                endif
                if ismatrix(Kv_idx)
                    ISS0.K(Kv_idx) = theta(numel(Av_idx)+numel(Bv_idx)+find(Kv_idx));
                endif
                if ismatrix(Cv_idx)
                    ISS0.C(Cv_idx) = theta(numel(Av_idx)+numel(Bv_idx)+numel(Kv_idx)+find(Cv_idx));
                endif
                if ismatrix(Dv_idx)
                    ISS0.D(Dv_idx) = theta(numel(Av_idx)+numel(Bv_idx)+numel(Kv_idx)+numel(Cv_idx)+find(Dv_idx));
                endif



                for i=1:maxStabCheck
                    if (~ is_stable([ISS0.A-ISS0.K*ISS0.C], 200*eps, 1))


                        if (i==1) printf ("\n"); endif
                        printf ("\r t = %d;  i = %d",t,i);
                        fflush (stdout);
                        #eig([ISS0.A-ISS0.K*ISS0.C])
                        #keyboard("\ninloop; unstable > ");

                        theta = theta_old + (.5 * (theta-theta_old));

                        if ismatrix(Av_idx)
                            ISS0.A(Av_idx) = theta(find(Av_idx));
                        endif
                        if ismatrix(Bv_idx)
                        ISS0.B(Bv_idx) = theta(numel(Av_idx)+find(Bv_idx));
                        endif
                        if ismatrix(Kv_idx)
                        ISS0.K(Kv_idx) = theta(numel(Av_idx)+numel(Bv_idx)+find(Kv_idx));
                        endif
                        if ismatrix(Cv_idx)
                        ISS0.C(Cv_idx) = theta(numel(Av_idx)+numel(Bv_idx)+numel(Kv_idx)+find(Cv_idx));
                        endif
                        if ismatrix(Dv_idx)
                        ISS0.D(Dv_idx) = theta(numel(Av_idx)+numel(Bv_idx)+numel(Kv_idx)+numel(Cv_idx)+find(Dv_idx));
                        endif

                        StabCheck(t) = -i;

                    else

                        StabCheck(t) = i;
                        break;

                    endif
                endfor

            endif # (another obs_avail block)
            
            ##
            ## Update P
            ##
            if (obs_avail)
                Qs0.P = 1/lam * (Qs0.P - L * S * L') + Qv;
            endif
            
            ##
            ## Calculate next state X(t+1)
            ##
            if (obs_avail)
                if (ins_avail)
                    X0 = ISS0.A * X0 + ISS0.B * uobs(t,:)' + ISS0.K * err;
                else
                    X0 = ISS0.A * X0 + ISS0.K * err;
                endif
            else
                if (ins_avail)
                    X0 = ISS0.A * X0 + ISS0.B * uobs(t,:)';
                else
                    X0 = ISS0.A * X0;
                endif
            endif        
            
            ##
            ## Update W 
            ## 
            if (obs_avail)
                Qs0.W = (ISS0.A-ISS0.K*ISS0.C) * Qs0.W + Mt - ISS0.K*Dt;
            endif

            ##
            ## Determine the derivative Dt (d/dtheta(Cx+Du))
            ## (I don't really understand this, but calculate Dt here instead of
            ##  *before* the Qs0.W update)
            ##
            if (obs_avail)
                for i=1:length(Cv_cidx)
                    Dt_C(Cv_cidx(i),i) = X0(ceil(Cv_idx/ny)(i));
                endfor
                if ismatrix(Dt_C)  # If zero-length matrix, don't do anything
                    Dt((ny*nv)-numel(Dt_D)-numel(Dt_C)+1:(ny*nv)-numel(Dt_D)) = Dt_C(:); 
                endif
                
                if (ins_avail)
                    for i=1:length(Dv_cidx)
                        Dt_D(Dv_cidx(i),i) = uobs(t,ceil(Dv_idx/ny)(i));
                    endfor
                else
                    Dt_D(:) = 0; # Not sure if this is correct approach or not
                endif
                
                if ismatrix(Dt_D)  # If zero-length matrix, don't do anything
                    Dt((ny*nv)-numel(Dt_D)+1:(ny*nv)) = Dt_D(:);
                endif
            endif
            
            ##
            ## Update Psi
            ##
            if (obs_avail)
                Qs0.Psi = Qs0.W' * ISS0.C' + Dt';
            endif
           
            #if (mod(t,100)==0)
            #    sys=ss(ISS0.A,ISS0.B,ISS0.C,ISS0.D,1);
            #    impulse(sys,1,25);
            #    eval(sprintf("keyboard(\"\\nInside Loop (%d)> \")",t)); 
            #endif
            
            
            ##
            ## Accumulate requested output matrices
            ##
            Y(t,:)=yhat';
            
            ISS=ISS0;
            X=X0;
            Qs=Qs0;
            
            
            if (nargout > 5)
                ISS_Mtrx(t) = struct('A',ISS0.A,'B',ISS0.B,'K',ISS0.K,'C',ISS0.C,'D',ISS0.D);
            endif
            
            if (nargout > 6)
                X_Mtrx(t) = struct('X',X0);
            endif
            
            if (nargout > 7)
                Qs_Mtrx(t).G = Qs0.G;
                if (isfield(Qs0,'Wss') )
                    Qs_Mtrx(t).W = Qs0.W;
                endif
                if (isfield(Qs0,'Psiss') )
                    Qs_Mtrx(t).Psi = Qs0.Psi;
                endif
                if (isfield(Qs0,'Pss') )
                    Qs_Mtrx(t).P = Qs0.P;
                endif
                
            endif


        endfor
        
        keyboard("\nEKF_ISS (after loop) > ");


endfunction
