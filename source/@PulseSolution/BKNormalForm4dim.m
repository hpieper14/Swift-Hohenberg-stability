function S = BKNormalForm4dim(S)


    r = S.vfParams.mu; 
    nu = S.vfParams.nu;
    gam = 38*nu^2/9-3;


    phi = S.nfBranch;
    
    
    x=(-S.time:.01:S.time).';
   
    
        
    N=max(size(x));
    sol=zeros(N,4);
    
 
    sol(:,1) = 2*sqrt(2.*r./gam).*sech(x.*sqrt(r)./2).*cos(x+phi);
    sol(:,2) = - (2.*sin(phi + x).*(2.*r/gam)^(1/2))./cosh((r^(1/2).*x)./2)...
        - (r^(1/2).*sinh(r^(1/2).*x./2).*cos(phi + x).*(2.*r./gam)^(1/2))./cosh(r^(1/2).*x./2).^2;
    sol(:,3) = (r.*sinh(r^(1/2).*x./2).^2.*cos(phi + x).*(2.*r./gam)^(1/2))./cosh(r^(1/2).*x./2).^3 ...
        - (r.*cos(phi + x).*(2.*r./gam)^(1/2))./(2.*cosh(r^(1/2).*x./2)) ...
        - (2.*cos(phi + x).*(2.*r./gam)^(1/2))./cosh(r^(1/2)*x./2) ...
        + (2.*r^(1/2).*sinh(r^(1/2).*x./2).*sin(phi + x).*(2.*r./gam)^(1/2))./cosh(r^(1/2).*x./2).^2;
    sol(:,4) = (2.*sin(phi + x).*(2.*r./gam)^(1/2))./cosh(r^(1/2).*x./2) ...
        + (3.*r.*sin(phi + x).*(2.*r./gam)^(1/2))./(2.*cosh(r^(1/2).*x./2)) ...
        + (3.*r^(1/2).*sinh(r^(1/2).*x./2).*cos(phi + x).*(2.*r./gam)^(1/2))./cosh((r^(1/2).*x)./2).^2 ...
        + (5.*r^(3/2).*sinh(r^(1/2).*x./2).*cos(phi + x).*(2.*r./gam)^(1/2))./(4.*cosh(r^(1/2).*x./2).^2) ...
        - (3.*r.*sinh(r^(1/2).*x./2).^2.*sin(phi + x).*(2.*r./gam)^(1/2))./cosh(r^(1/2).*x./2).^3 ...
        - (3.*r^(3/2).*sinh(r^(1/2).*x./2).^3.*cos(phi + x).*(2.*r./gam)^(1/2))./(2.*cosh(r^(1/2).*x./2).^4);

    S.nfData.time = x; 
    S.nfData.sol = sol; 
end
