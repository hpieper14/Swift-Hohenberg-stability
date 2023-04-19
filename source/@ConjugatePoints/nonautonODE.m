function f1 = nonautonODE(C, T, Y)
    params = C.vfParams; 
    lam=params.lambda;
    mu=params.mu;
    nu=params.nu;

    f1=Y(2);
    f2=Y(3);
    f3=Y(4);
    f4= -2.*Y(3) - (mu+1).*Y(1) + nu.*Y(1).^2 - Y(1).^3;

    g1 = Y(8);
    g2 = Y(7) - 2.*Y(8);
    g3=(lam - 1 - mu + 2*nu.*Y(1) - 3.*Y(1).^2).*Y(5);
    g4=Y(6);

    f1=[f1;f2;f3;f4;g1;g2;g3;g4];


end
