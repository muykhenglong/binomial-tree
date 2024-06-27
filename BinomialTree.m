function Price = BinomialTree(S,K,T,r,u,d,p,N,IsCall,IsAmer)

% function `BinomialTree` builds a pricing lattice and returns the `Price` of a European or American style call or put option.
%  - Takes inputs:
    % - `S`: Current stock price
    % - `K`: Strike price of the option
    % - `T`: Time to expiration in years
    % - `r`: Risk-free rate
    % - `u`: Up-factor in binomial model
    % - `d`: Down-factor in binomial model
    % - `p`: Probability of an up move in the binomial model
    % - `N`: Number of time steps in the binomial model
    % - `IsCall`: Boolean indicating if the option is a call (1) or put (0)
    % - `IsAmer`: Boolean indicating if the option is American (1) or European (0)

% Handle Call or Put
if IsCall
    Intrinsic = @(S) max(0,S-K);
else
    Intrinsic = @(S) max(0,K-S);
end

% CRR Parameters
dT = T/N;

%Allocate Memory
Lattice = zeros(N+1);

%Intrinsic Value at Expiration
for i = 0:N % Loops through each possible state at expiration
    Lattice(i+1,N+1) = Intrinsic(S*u^i*d^(N-i)); % Computes the intrinsic value at each node at expiration.
end

% Backwards Induction through Tree
for j = N:-1:1 % Loops backwards through the tree from time \(N-1\) to 0 to calculate the option values
    for i = 1:j % Loops through each node at each time step
        Lattice(i,j) = max(IsAmer * Intrinsic(S*u^(i-1)*d^(j-i)) , (p * Lattice(i+1,j+1) + (1-p) * Lattice(i,j+1)) * exp(-r*dT));
    end
end
Price = Lattice(1,1);

end
