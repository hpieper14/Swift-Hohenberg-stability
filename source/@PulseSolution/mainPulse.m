% ADDME  Add two values together.
%   C = ADDME(A) adds A to itself.
%
%   C = ADDME(A,B) adds A and B together.
%
%   See also SUM, PLUS.
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