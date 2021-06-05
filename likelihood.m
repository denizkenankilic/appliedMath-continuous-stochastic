M = 100000; % M - number of Monte Carlo paths (number of W's)
r = 0.05; % r - interest rate
sigma = 0.2; % sigma - volatility
T = 1; % T - time interval
K = 100; % K - strike price
S0 = 100; % S0 - initial asset value
W = sqrt(T)*randn(M,1); % Brownian motion
S = S0*exp((r-0.5*sigma^2)*T+sigma*W); % integrated scalar geometric Brownian Motion 
f = exp(-r*T)*(0.5)*(1 + sign(S-K)); % payoff function
% delta
g = f .* W/(S0*sigma*T);
gbar = sum(g)/M;
sig2 = sum(g.^2)/M - gbar^2;
gexact = digital_call(r,sigma,T,S0,K,'delta');
disp(sprintf('Analytic Delta: %f',gexact))
disp(sprintf('LRM Delta: %f +/- %f',gbar,3*sqrt(sig2/M)))
% vega
g = f .* ( (W.^2-T)/(sigma*T) - W);
gbar = sum(g)/M;
sig2 = sum(g.^2)/M - gbar^2;
gexact = digital_call(r,sigma,T,S0,K,'vega');
disp(sprintf('Analytic Value: %f',gexact))
disp(sprintf('LRM Vega: %f +/- %f',gbar,3*sqrt(sig2/M)))
