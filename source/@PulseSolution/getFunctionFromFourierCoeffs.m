function sol = getFunctionFromFourierCoeffs(S,coeffs, interval_type)
order = S.fourier.order; 

if interval_type == "half"
    T=0:.05:S.time;
else
    T=-S.time:.05:S.time;
end

f = 0;
for n = -order:order
    f = f+real(coeffs(n+order+1)*exp(1i*n*pi.*T./S.time));
end

sol(:,1) = T';
sol(:,2) = f';

end
