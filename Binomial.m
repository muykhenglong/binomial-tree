function Price = Binomial(S,K,T,r,vol,q,N,IsCall,IsAmer,Method)

dT = T/N;

if strcmp(Method,'EQP')
    %EQP specification
    u=exp((r-q-(vol^2)/2)*dT + vol*sqrt(dT));
    d=exp((r-q-(vol^2)/2)*dT - vol*sqrt(dT));
    p=0.5;
 
elseif strcmp(Method,'LR')
    d1 = (log(S/K)+(r-q+vol^2/2)*T)/(vol*sqrt(T));
    d2 = d1 - vol * sqrt(T);
    p_prime = 0.5 + sign(d1) * (0.25 - 0.25 * exp(-((d1/(N+1/3))^2*(N+1/6))))^.5;
    p = 0.5 + sign(d2) * (0.25 - 0.25 * exp(-(d2/(N+1/3))^2*(N+1/6)))^.5;
    u = (p_prime * exp(r*dT)) / p;
    d = (exp(r*dT) - u * p) / (1 - p);

else
    %CRR specification
    u=exp(vol * sqrt(dT));
    d=1/u;
    if strcmp(Method,'TIAN')
        %Strike exactly on Node - TIAN
        j=ceil((log(K/S)- N*log(d))/(log(u/d)));
        tilt = (K/(S*(u^j)*(d^(N-j))))^(1/N);
        u=u*tilt;
        d=d*tilt;
    end
    p=(exp((r-q)*dT) - d)/(u-d);
end
Price = BinomialTree(S,K,T,r,u,d,p,N,IsCall,IsAmer);
end