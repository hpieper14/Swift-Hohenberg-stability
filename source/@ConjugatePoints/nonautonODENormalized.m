function f1 = nonautonODENormalized(C, T, Y)
    params = C.vfParams; 
    lam=params.lambda;
    mu=params.mu;
    nu=params.nu;

    f1=Y(2);
    f2=Y(3);
    f3=Y(4);
    f4= -2.*Y(3) - (mu+1).*Y(1) + nu.*Y(1).^2 - Y(1).^3;
    
    B_x = [0,0,0,1;0,0,1,-2;lam-1-mu+2*nu.*Y(1)-3.*Y(1).^2,0,0,0;0,1,0,0];

    
    barw = [Y(5); Y(6); Y(7); Y(8)];
    
    vec = B_x*barw - barw.*dot(B_x*barw, barw);


    f1=[f1;f2;f3;f4;vec];
end
