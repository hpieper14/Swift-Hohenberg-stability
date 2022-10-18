% returns a vector of coeffs (a_{-k},..., a_{-1}, a_0, a_1, .... a_k)
% order = n
function S = getFullFourierCoeffs(S)
    S = BKNormalForm4dim(S); 
    order = S.order;
    L = S.time;

    sol = S.nfData.sol;
    time = S.nfData.time;

    sol = sol(:,1); 


    coeffs=zeros(1, 2*order+1);
    for n = -order:order
        F = sol.*exp(-1i*n*pi.*time./L);
        coeffs(n+order+1) = real(trapz(time, F)/(2*L));
    end     
    S.fourierCoeffs = coeffs; 
end