function Price = BinomialTree(S,K,T,r,u,d,p,N,IsCall,IsAmer)
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
for i = 0:N
    Lattice(i+1,N+1) = Intrinsic(S*u^i*d^(N-i));
end

% Backwards through Tree
for j = N:-1:1
    for i = 1:j
        Lattice(i,j) = max(IsAmer * Intrinsic(S*u^(i-1)*d^(j-i)) , (p * Lattice(i+1,j+1) + (1-p) * Lattice(i,j+1)) * exp(-r*dT));
    end
end
Price = Lattice(1,1);

end