classdef ConjugatePoints 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Class properties.                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties( Access = public )

        conjPts = [];         % struct containing data associated to the conjugate points  
        vfParams = [];         % Parameters for the vector field
        Euminus = []; 
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             Constructor.                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
    
    methods( Access = public, Static = false )
        
        function C = ConjugatePoints(conjPts, vfParams)
         
            % Store approximation order.
            
            % Store parameters for the vector field.
            C.vfParams = vfParams;
            C.conjPts = conjPts; 
                                    
            
        end
        
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              Methods.                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods( Access = public, Static = false )
        
        mat = Binf(C); 
        [vectors, values]= getBinfEigs(C);
        isLagrangian = verifyLagrangian(C, frame);

        C = generateEuFrame(C);
            
        vals= calculateDeterminant(C);

        C = getConjPtLocations(C, vals); 

        [frame, time]= makeFrame(C, basis_1, basis_2)
        C = get4DimIC(C, S)


                                    
    end   
    methods( Access = public, Static = true )
                
        
                                    
    end   
    
end

