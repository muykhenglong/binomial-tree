function Price = Binomial(S,K,T,r,vol,q,N,IsCall,IsAmer,Method)

% function `Binomial` returns the `Price` of an option
% Takes inputs:
    % - `S`: Current stock price
    % - `K`: Strike price of the option
    % - `T`: Time to expiration in years
    % - `r`: Risk-free rate
    % - `vol`: Volatility of the stock
    % - `q`: Dividend yield
    % - `N`: Number of time steps in the binomial model
    % - `IsCall`: Boolean indicating if the option is a call (1) or put (0)
    % - `IsAmer`: Boolean indicating if the option is American (1) or European (0)
    % - `Method`: String indicating which binomial model to use (`EQP`, `LR`, `CRR`, `TIAN`)

dT = T/N;

% Checks if the method specified is 'EQP' (Equal Probabilities)
if strcmp(Method,'EQP')
    %EQP specification
    u=exp((r-q-(vol^2)/2)*dT + vol*sqrt(dT));
    d=exp((r-q-(vol^2)/2)*dT - vol*sqrt(dT));
    p=0.5;

% Checks if the method specified is 'LR' (Leisen-Reimer)
elseif strcmp(Method,'LR')
    d1 = (log(S/K)+(r-q+vol^2/2)*T)/(vol*sqrt(T));
    d2 = d1 - vol * sqrt(T);
    p_prime = 0.5 + sign(d1) * (0.25 - 0.25 * exp(-((d1/(N+1/3))^2*(N+1/6))))^.5;
    p = 0.5 + sign(d2) * (0.25 - 0.25 * exp(-(d2/(N+1/3))^2*(N+1/6)))^.5;
    u = (p_prime * exp(r*dT)) / p;
    d = (exp(r*dT) - u * p) / (1 - p);

else
% Default method is 'CRR' (Cox-Ross-Rubinstein)
    %CRR specification
    u=exp(vol * sqrt(dT));
    d=1/u;
    
    if strcmp(Method,'TIAN')
        % Adjusts up and down factors if the method is 'TIAN'
        % Strike exactly on Node - TIAN
        j=ceil((log(K/S)- N*log(d))/(log(u/d)));
        tilt = (K/(S*(u^j)*(d^(N-j))))^(1/N);
        u=u*tilt;
        d=d*tilt;
    end
    p=(exp((r-q)*dT) - d)/(u-d);
end
Price = BinomialTree(S,K,T,r,u,d,p,N,IsCall,IsAmer);
end
