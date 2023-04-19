classdef PulseSolution 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Class properties.                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties( Access = public )

        fourier = [];         % struct containing the fourier coefficients 
                                % and associated parameters for Newtons method    
        vfParams = [];         % Parameters for the vector field       
        time = [];              % L for [-L,L] time interval to compute approximation on 
        normalForm = [];        % struct containing normal form data and branch value
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             Constructor.                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
    
    methods( Access = public, Static = false )
        
        function S = PulseSolution(fourier, vfParams, normalForm, time)
         
            % Store approximation order.
            S.fourier = fourier;
        
            % Store parameters for the vector field.
            S.vfParams = vfParams;

            % Store value of normal form branch 
            allowedBranchValues = [0,pi];
            if (ismember(normalForm.branch, allowedBranchValues) == 1) == 0
                error('The branch value is invalid. Please choose branch value 0 or pi.')
            end
            S.normalForm = normalForm;
            
            % Store temporal domain, should be a positive number for domain
            % [-L, L] 
            S.time = time;
                                    
            
        end
        
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              Methods.                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods( Access = public, Static = false )

        S = mainPulse(S)
        
        % Normal form solution via Burke and Knobloch 
        S = BKNormalForm4d_halfline(S)

        % Get normal form solution with Dirichlet/Neumann BC 
        S = trimNFSol_halfline(S);

        % Get fourier coefficients from BK normal form solution 
        S = getFullFourierCoeffs(S, sol, interval_type)

        % Perform Newton's method 

        S = Newton_halfline(S); 
        S = Newton(S)

        % get pulse derivative initial conditions in symplectic coordinates
        % 
        S = getPulseDerivIC(S, t_0)

        % helper methods, could be made private 
        % helper method to recover function from Fourier coefficients 
        solution = getFunctionFromFourierCoeffs(S, coeffs, time_vec)
        jacobian = DFFourier(S, a)
        fourier_vf = fourierODE(S, a)


                                    
    end   
    methods( Access = public, Static = true )
                
        
                                    
    end   
    
end

