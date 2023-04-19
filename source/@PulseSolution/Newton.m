function S=Newton(S)
    s = S.fourier.full_coeffs;     
    % this should be placed as an attribute of S 
    tol = 1e-10; 

    side = 2;
    k=0;
    while k < 100
        fcn=S.fourierODE(s);
        DF=S.DFFourier(s);
        if vecnorm(fcn)< tol 
            break
        end
        s = s - (DF\(fcn'))';
        k=k+1;
    end
    S.fourier.full_coeff=s;

end