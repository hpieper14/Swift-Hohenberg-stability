classdef PulseSolution 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Class properties.                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties( Access = public )

        fourierCoeffs = [];         % Fourier coefficients of approximation     
        order = [];   % Order of Fourier approximation
        vfParams = [];         % Parameters for the vector field       
        nfBranch = [];                % Branch for normal form solution
        time = [];              % L for [-L,L] time interval to compute approximation on 
        nfData = [];
                               
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             Constructor.                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
    
    methods( Access = public, Static = false )
        
        function S = PulseSolution(order, vfParams, nfBranch, time)
         
            % Store approximation order.
            S.order = order;
        
            % Store parameters for the vector field.
            S.vfParams = vfParams;

            % Store value of normal form branch 
            allowedBranchValues = [0,pi];
            if (ismember(nfBranch, allowedBranchValues) == 1) == 0
                error('The branch value is invalid. Please choose branch value 0 or pi.')
            end
            S.nfBranch = nfBranch;
            
            % Store temporal domain, should be an array with right/left
            % endpts
            S.time = time;
                                    
            
        end
        
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              Methods.                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods( Access = public, Static = false )
        
        % Normal form solution via Burke and Knobloch 
        S = BKNormalForm4dim(S);

        % Get normal form solution with Dirichlet/Neumann BC 
        S = trimNFSol(S)

        % Get fourier coefficients from BK normal form solution 
        S = getFullFourierCoeffs(S)

        % 
        sol = getFunctionFromFourierCoeffs(S)
                                    
    end   
    methods( Access = public, Static = true )
                
        
                                    
    end   
    
end

