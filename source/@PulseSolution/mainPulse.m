% this function computes the fourier approximation of the pulse solution,
% the initial conditions for U_\phi and the unstable eigenvalues of the
% pulse solution
function S = mainPulse(S)
    S = S.BKNormalForm4d_halfline();
    S = S.trimNFSol_halfline();
    S = BKNormalForm4d_halfline(S); 
    
    
    S = S.Newton_halfline(); 
    S = S.Newton();
    
    DF=S.DFFourier(S.fourier.full_coeffs);
    
    D=eig(DF);
    [i,j] = find(D>1e-12);
    vals = [];
    if min(size(D(i,j))) == 1
        vals = [vals, D(i,j)];
    else 
        A = D(i,j);
        vals = [vals, A(:, 1)];
    end
    
    S.fourier.unstable_eigs = vals;
end