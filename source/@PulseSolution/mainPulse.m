% MAINPULSE  Computes Fourier approximation of pulse solution and
% approximates the unstable eigenvalues associated to the pulse solution.
%   S = S.mainPulse()
%   S = mainPulse(S) 
% 
function S = mainPulse(S)
    % get pulse approximation from normal form. Determines where to
    % truncate spatial domain so $\varphi'(L) \approx 0$. 
    S = S.BKNormalForm4d_halfline();
    S = S.trimNFSol_halfline(); 
    S = BKNormalForm4d_halfline(S); 
    
    % refine Fourier coefficients with Newton's method
    S = S.Newton_halfline(); 
    
    % Redefine the Fourier Coefficients
    if S.vfParams.mu == .2
        load('StableCoeff.mat','StableCoeff')
        % StableCoeff
        N_fourier = 2*S.fourier.order+1; 
        N_stored = length(StableCoeff);
        if N_fourier  == N_stored 
            % Everything is good
        elseif N_fourier  > N_stored 
            % Pad Stored with zeros
            n_padding = (N_fourier -N_stored )/2;
            StableCoeff = [zeros(1,n_padding)  , StableCoeff , zeros(1,n_padding)];
        elseif N_fourier  < N_stored 
            % Cut Stored down to size
            n_padding = (N_stored -N_fourier  )/2;
            StableCoeff =StableCoeff (1+n_padding:end-n_padding);
        end
        S.fourier.full_coeffs = StableCoeff;
    end

    S = S.Newton();

    % get Jacobian 
    DF = S.DFFourier(S.fourier.full_coeffs);
    
    % compute eigenvalues of Jacobian and save positive eigenvalues
    D = eig(DF);
    [i,j] = find(D>1e-12);
    vals = [];
    if length(i) >0
        if min(size(D(i,j))) == 1
            vals = [vals, D(i,j)];
        else 
            A = D(i,j);
            vals = [vals, A(:, 1)];
        end
    end
    
    S.fourier.unstable_eigs = vals;
end