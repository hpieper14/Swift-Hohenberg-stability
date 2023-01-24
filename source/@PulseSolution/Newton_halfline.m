function S = Newton_halfline(S) 
options=optimset('Display','iter','Jacobian','on','MaxIter',10000);     % Option to display output and use Jacobian

u = S.normalForm.sol(:,1); 
[uout,fval] = fsolve(@(u) S.fourierODE_halfline(u),u,options);  
coeffs = S.getFullFourierCoeffs(uout, S.normalForm.time);
S.fourier.half_coeffs = coeffs; 
T = 1:1:max(size(uout));
figure
plot(T, uout)
full_uout = [flip(uout); uout];
full_time = [-flip(S.normalForm.time); S.normalForm.time];
figure 
plot(full_time, full_uout)
title("full solution")
S.fourier.full_coeffs = S.getFullFourierCoeffs(full_uout, full_time);
end
