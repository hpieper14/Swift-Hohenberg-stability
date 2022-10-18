function test = testgetFullFourierCoeffs()
    order = 10;
    params.nu = 1.2;
    params.mu = .1;
    vfParams = params;
    nfBranch = 0; 
    time = 10; 

    S = PulseSolution(order, vfParams, nfBranch, time);
    
    coeffs = getFullFourierCoeffs(S);

    T=-L:.25:L;
    f = 0;
    for n = -order:order
        f = f+real(coeffs(n+order+1)*exp(1i*n*pi.*T./L));
    end
 
    
    figure
    plot(time, sol);
    hold on
    plot(T,f)
    legend('Sol via NF equation','Fourier approximation')

end