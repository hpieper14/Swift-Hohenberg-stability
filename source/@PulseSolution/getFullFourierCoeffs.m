% returns a vector of coeffs (a_{-k},..., a_{-1}, a_0, a_1, .... a_k)
% order = n
function coeffs = getFullFourierCoeffs(S, u, time_vec)
    order = S.fourier.order;
    L = (time_vec(end) - time_vec(1))/2;
 


    coeffs=zeros(1, 2*order+1);
    for n = -order:order
        F = u.*exp(-1i*n*pi.*time_vec./L);
        coeffs(n+order+1) = real(trapz(time_vec, F)/(2*L));
    end     
end