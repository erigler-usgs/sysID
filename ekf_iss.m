###
### ekf_iss.m
###
### This function implements an Extended Kalman Filter in the
### innovations form of the state space model (ISS, see "Asymptotic
### Behavior of the Extended Kalman Filter as a Parameter Estimator for
### Linear Systems", IEEE Trans. Automat. Contr., Ljung, 1979) in order
### to determine both the optimal state estimate and linear parameters
### for properly formed matrices that can be used as follows:
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
### [Y_mtrx, X_mtrx, Theta_mtrx] = ekf_iss (its, ots, iss, X0, Theta0, iforget)
###
### ...where...
###
### INPUTS:
###
### its        - Input time-series matrix.  Columns are different inputs
###              and rows are sequential input observations.
###
### ots        - Output time-series matrix.  Columns are different outputs
###              and rows are sequential output observations.
###
### iss        - Structure of state space matrices in innovations
###              form.If the individual matrix elements are actual
###              numbers, then they are fixed values.  If the elements
###              are to be adjusted for identification, they should be
###              filled with NaNs.  The NaNs will then be filled with
###              values from the parameter vector, Theta, in the
###              following manner (example is for a 4th order, 2 input,
###              3 output system in canonical observable form):
###
###      iss.A = [ 0   1   0   0 ] iss.B = [b01 b02] ssi.K = [k01 k02 k03]
###              [a01 a02 a03 a04]         [b03 b04]         [k04 k05 k06]
###              [a05 a06 a07 a08]         [b05 b06]         [k07 k08 k09]
###              [a09 a10 a11 a12]         [b07 b08]         [k10 k11 k12]
###
###      iss.C = [ 1   0   0   0 ] iss.D = [d01 d02]
###              [ 0   1   0   0 ]         [d03 d04]
###              [ 0   0   1   0 ]         [d05 d06]
###
###              In addition, the two following structure elements
###              correspond to the innovations covariance and the
###              diagonal of the parameter process noise covariance
###              matrices:
###
### iss.Lambda = [l11 l12 l13]   iss.Qp = diag ([q01  0  ...  0 ]
###              [l21 l22 l23]                  [ 0  q02 ...  0 ]
###              [l31 l32 l33]                  [ 0   0   .   0 ]
###                                             [ 0   0  ... qnT] )
###                                            (nT = length(Theta))
###
###              These parameters are not "user-adjustable", however a
###              forgetting factor can be specified for Lambda (see
###              description of "iforget" input parameter).  Qp will
###              be implicitly set to zero for time-steps corresponding
###              to bad or missing I/O data (NaNs in its or ots)*.
###
###              * I'm not sure if this is best or not; it might be better
###                to have the parameters sort of decay exponentially to
###                the average time-stationary values.  This constitutes a
###                kind of dynamic Gauss-Markov process, but the background
###                "state" is still determined by a dynamical set of
###                equations, not some static or quasi-static climatology.  
###                Lets get this working first with Qp=0 for missing data,  
###                and return to this question later -EJR (10/13/04)
###
### Theta0     - This vector of parameters comprises the "extended" state
###              vector.  It is comprised of values that correspond to the
###              NaNs in the system matrices, in the order shown above:
###
###      Theta = [a01 a02 ... b01 b02 ... k01 k02 ... c01 c02 ... d01 d02 ...]'           
###         
### X0         - Initial state of the system.  If this vector is not
###              provided, it is set to a vector of zeros.
###
### iforget    - Innovations covariance "forgetting factor".  The algorithm
###              employed here uses the following recursive relationship to
###              determine Lambda from the data/model:
###                                                   T
###  Lambda(t) = Lambda(t-1) + 1/iforget * (err(t)*err (t) - Lambda(t-1))
###
###              "iforget" can be a scalar, or fixed non-normalized forgetting
###              factor, a vector of scalars for bootstrapping the algorithm
###              (iforget(i)=t, for example), or it can be Inf, which implies
###              that Lambda should be considered fixed in time [default].
###
###
### OUTPUTS:
###
###
### Y_mtrx     - Predicted output vector time series.  Columns correspond to
###              each Y(t+1), while rows correspond to multiple outputs.
###
### X_mtrx     - Propogated state vector time-series.  Columns correspond to
###              each X(t+1), while rows correspond to components of the Nth
###              order state vector.  
###
### Theta_mtrx - Theta vector time-series.  Columns correspond to each 
###              Theta(t), while rows correspond to the Theta components
###              described above for Theta0.  In other words, each column
###              of Y_mtrx and X_mtrx were predicted using the same column
###              of Theta_mtrx.  This means that Theta_mtrx will necessarily
###              have one less column than Y_mtrx and X_mtrx.
###
###              * I Need to develop a simple support routine to transform
###                Theta_mtrx into time-series of each polynomial matrix
###                (A,B,C,D,K) given the structure iss.  Or perhaps I should
###                just integrate that directly into this function, and take
###                A0, B0, C0, D0, and K0 as inputs, and provide A_mtrx, B_mtrx,
###                C_mtrx, D_mtrx, and K_mtrx as output (which of course will
###                 require considerably more storage space).


function [Th_mtrx, X_mtrx, I_mtrx] = ekf_iss (its, ots, iss, Th0, X0, iforget)



endfunction
