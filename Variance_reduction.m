% Variance Reduction

% useful site --> https://www.uni-ulm.de/fileadmin/website_uni_ulm/mawi.inst.110/lehre/ws15/MonteCarloMethods/Lecture_Notes.pdf

%The estimation of performance measures in Monte Carlo simulation can be made
%more efficient by utilizing known information about the simulation model. The
%more that is known about the behavior of the system, the greater the amount
%of variance reduction that can be achieved. The main variance reduction techniques discussed in this chapter are:
%1. Antithetic random variables.
%2. Control variables.
%3. Conditional Monte Carlo.
%4. Importance sampling.

%Estimating Probabilities

%ell=P(X in A)=EI(X in A)
%ell_hat=(1/n)sum_{i=1}^{N} I(X_i in A)
%RE=1/sqrt(N) (sqrt(ell(1-ell)))/ell
%P(|ell_hat-ell|/ell>eps)<=(1-ell)/(N*eps^2*ell), where eps=0.1
%ex1.(Estimating an exceedance probability (and estimating the relative error)).
%Suppose we wish to estimate ? = P(X > gamma), where gamma ? N and X ? Poi(lambda).

N = 10^4;
lambda = 10; gamma = 16;
 indics = zeros(N,1);
for i = 1:N
X = poissrnd(lambda);
indics(i) = (X > gamma);
 end
 ell = 1 - poisscdf(gamma, lambda) %ell=0.0270
 ell_hat = mean(indics)   %ell_hat=0.0217
 RE = sqrt(ell * (1 - ell)) / (ell * sqrt(N))  %RE=0.0600
 RE_hat = std(indics) / (ell_hat * sqrt(N))   %Re-_hat=0.0671

 %Variance Reduction Background: the simulation error for some
%parameter ? ? X¯ , depends on MSE = V ar(X¯ ) = Var(X)/n,
%so the simulation can be more efficient if Var(X) is reduced.


 % Antithetic Random Variables
%A pair of real-valued random variables (X1, X2) is called an antithetic
%pair if X1 and X2 have the same distribution and are negatively correlated.
%Key idea: if X1 and X2 are id RVs with mean ?,Var((X1 + X2)/2) = 1/4(Var(X1)+Var(X2)+ 2Cov(X1, X2)),
%so variance is reduced if X1 and X2 have Cov(X1, X2) ? 0.
%For many simulations, a ? estimator is X1 = h(U1, . . . , Un)for some h; 
%so consider the antithetic estimator X2 = h(1 ? U1, . . . , 1 ? Un).
%Combined estimator is (X1 + X2)/2.

%ex.2 
%Estimate the integral V =int_{0}^{inf} ln(1_x^2)e^(-x^2)dx.
%Use t = 1 ? e^(?x)// the function.

N = 1000; g = @(t)log(1+log(1-t).^2);
T = rand(1,N); X = g(T); % simple MC
disp([mean(X) std(X) 2*std(X)/sqrt(N)])
%0.66848 0.73278 0.046345
T = rand(1,N/2); X = ( g(T) + g(1-T) )/2;
disp([mean(X) std(X) 2*std(X)/sqrt(N/2)])
%0.67881 0.30933 0.027667

%ex.3
%Estimate V =int_{0}^{pi/4}int_{0}^{pi/4} x^2y^2sin(x + y) ln(x + y)dxdy;
%notice integrand is monotone increasing in x and y. Use x = ?u1/4, y =?u2/4.

N = 1000; U = rand(2,N);
f = @(x)prod(x).^2.*sin(sum(x)).*log(sum(x));
X = pi^2*f(pi*U/4)/16; % simple MC
disp([mean(X) std(X) 2*std(X)/sqrt(N)])
%0.0037452 0.012534 0.00079271
U = rand(2,N/2);
X = pi^2*( ( f(pi*U/4) + f(pi*(1-U)/4) )/2 )/16;
disp([mean(X) std(X) 2*std(X)/sqrt(N/2)])
%0.0037885 0.0085332 0.00076323

% Control variables
%Assume desired simulation quantity is ? = E[X].There is another simulation RV Y with known µY = E[Y].
%For any c, RV Z = X + c(Y ? µY ), is an unbiased estimator of ?, because E[Z] = E[X] + c(E[Y ] ? µY ) = ?. 
%Now, Var(Z) = Var(X+cY ) = Var(X)+c^2Var(Y)+2cCov(X,Y) is minimized when c = c*= ?Cov(X,Y)/Var(Y), and
%Var(X + c*Y)= Var(X)? Cov(X,Y)^2/Var(Y).
%Y is called a control variate for X; in order to reduce variance, choose a Y correlated with X.
%Implementation: Cov(X,Y), Var(Y) estimated from data.
%Cov(X,Y)?Covˆ(X,Y)=1/(n ? 1)sum_{i=1}^{n}(Xi ? X¯ )(Yi ? Y¯ ),
%Var(Y)?Varˆ(Y) = 1/(n-1)sum_{i=1}^{n}(Yi ? Y¯)^2, so
%c* ? cˆ* = ? {sum_{i=1}^{n}(Xi ? X¯)(Yi ? Y¯)}/{sum_{i=1}^{n}(Yi ? Y¯)^2},
%goal is to choose Y so that Y is ? X, with Y is easy to simulate and µY is easy to find.


%ex.4
%Estimate theta=int_{0}^{1}e^x dx. 
%Then muY=1/2 , Var(U)=1/12, Var(e^U)=.242 
%cov(e^U,U)=0.141086
%Var(X+c*Y)=0.0039

N = 1000; U = rand(1,N); X = exp(U);
disp([mean(X) std(X) 2*std(X)/sqrt(N)])
%1.7488 0.49766 0.031475
Y = U; muY = 1/2;
Xb = mean(X); Yb = mean(Y);
cs = -sum( (X-Xb).*(Y-Yb) )/sum( (Y-Yb).^2 );
Z = X + cs*( Y - muY );
disp([mean(Z) std(Z) 2*std(Z)/sqrt(N)])
%1.7189 0.063119 0.003992

 
%ex.5
%Estimate V=int_{0}^{2} e^(-x^2)dx=2int_{0}^{1}e^(-(2u)^2)du
% Try Y=2e^(-2u).

N = 1000; U = rand(1,N);
X = 2*exp(-(2*U).^2);
disp([mean(X) std(X) 2*std(X)/sqrt(N)])
%0.91666 0.68842 0.04354
Y = 2*exp(-2*U); muY = 1 - exp(-2);
Xb = mean(X); Yb = mean(Y);
cs = -sum( (X-Xb).*(Y-Yb) )/sum( (Y-Yb).^2 );
Z = X + cs*( Y - muY );
disp([mean(Z) std(Z) 2*std(Z)/sqrt(N)])
%0.88474 0.13402 0.0084762


%ex. 6
%Estimate V =int_{0}^{pi/4}int_{0}^{pi/4} x^2y^2sin(x + y) ln(x + y)dxdy.
%Use x^2y^2 to "approximate" the integrand.

N = 1000; U = rand(2,N);
f = @(x)prod(x).^2.*sin(sum(x)).*log(sum(x));
X = pi^2*f(pi*U/4)/16; Xb = mean(X);% simple MC
disp([mean(X) std(X) 2*std(X)/sqrt(N)])
%0.0035267 0.011023 0.00069718
g = prod(X).^2;
Y = pi^2*g(pi*U/4)/16; Yb = mean(Y);
cs = -sum( (X-Xb).*(Y-Yb) )/sum( (Y-Yb).^2 );
muY = (pi/4)^6/9; Z = X + cs*( Y - muY );
disp([mean(Z) std(Z) 2*std(Z)/sqrt(N)])
%0.0037289 0.0043002 0.00027197

 % Imporatnce sampling
 %ex.7
%suppose we want to estimate ? = P(X ? [5, ?)), where X ? Exp(1). We could make the event of interest more likely by choosing g to be the
%density of a Exp(1/5) random variable (so that EgX = 5). Suppose we have f(x) = e?x and g(x) = ?e??x, then
% f(x)/g(x)=1/? exp{(? ? 1)x}.
%This approach is implemented in Matlab for the general problem of estimating ? = P(X > ?),
%using an exponential distribution with ? = 1/?.

N = 10^5; gamma = 5;
lambda = 1/gamma;
results = zeros(N,1);
results_IS = zeros(N,1);
 for i = 1:N
 X = -log(rand);
 X_IS = -log(rand) / lambda;
results(i) = (X > gamma);
 results_IS(i) = exp((lambda-1)*X_IS) / lambda * (X_IS > gamma);
end
 ell = exp(-5)    %0.0067
 ell_hat = mean(results) %0.0064
 re_hat = std(results) / (sqrt(N) * ell_hat) %0.0394
 ell_hat_IS = mean(results_IS)  %0.0068
 re_hat_IS = std(results_IS) / (sqrt(N) * ell_hat_IS) %0.0081
 re_hat / re_hat_IS  %4.8720
 
 
%Conditional MC
% Let Y be a random variable and Z a random vector. Then
% Var(Y)=E Var(Y|Z)+Var(E[Y|Z]) and hence Var(E[Y|Z]<=Var(Y)).
% let l=EY,Y input from MC experiment, Z~g such the conditionaal
% expectation E[Y|Z=z] can be computed analytically. By l=EY=EE[Y|Z],
% it follows that E[Y|Z] is an unbiased estimator of l.

%https://www.mathworks.com/help/econ/monte-carlo-simulation-of-conditional-variance-models.html
