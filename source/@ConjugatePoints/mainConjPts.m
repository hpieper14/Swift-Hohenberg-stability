% this function computes the conjugate points for pulse solution S 
function [S, C] = mainConjPts(C, S) 
    L = C.conjPts.L; 
    S = S.getPulseDerivIC(L); 

    C = C.get4DimIC(S);
    ic = S.fourier.pulseDeriv.coords;
    C.Euminus.pulse_deriv_ic = ic/vecnorm(ic);
    
    
    [vectors, values]= C.getBinfEigs();
        
    C = C.generateEuFrame();
    C = C.calculateDeterminant();

end
