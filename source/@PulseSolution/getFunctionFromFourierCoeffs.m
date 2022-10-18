function sol = getFunctionFromFourierCoeffs(S)
coeffs = S.fourierCoeffs;
order = S.order; 

T=-S.time:.05:S.time;
f = 0;
for n = -order:order
    f = f+real(coeffs(n+order+1)*exp(1i*n*pi.*T./S.time));
end

sol(:,1) = T';
sol(:,2) = f';

end
