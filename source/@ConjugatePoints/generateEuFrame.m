% GENERATEEUFRAME  Generates the frame for $E^u_-(x;0)$ for various x
% values and saves it to the ConjugatePoints object C 
% 
%   C = generateEuFrame(C)
%   C = C.generateEuFrame()
function C = generateEuFrame(C)
    L = C.conjPts.L; 
    orig_ic = C.Euminus.pulse_ic';
    deriv_ic = C.Euminus.pulse_deriv_ic';
    params = C.vfParams; 
    
    % inialize initial condiction for the basis solution of $E^u_-$ that
    % is not constructed from the derivative of the pulse
    [vectors, values]= getBinfEigs(C);
    v1=real(vectors.u(:,1));    

    if C.Euminus.normalize == 1 
        deriv_ic = deriv_ic./vecnorm(deriv_ic);
        v1 = v1./vecnorm(v1);
    end

    new_ic1 = [orig_ic; v1];
    new_ic2 = [orig_ic; deriv_ic]; 

    [full_phi1, basis_1] = C.integrateDE(new_ic1,-L,L);
    [full_phi1, basis_2] = C.integrateDE(new_ic2,-L,L);

    time = full_phi1(:, 1);
    figure 
    tiledlayout(4,1)
    nexttile 
    plot(time, full_phi1(:, 2))
    nexttile 
    plot(time, full_phi1(:, 3))
    nexttile 
    plot(time, full_phi1(:, 4))
    nexttile 
    plot(time, full_phi1(:, 4))
    title("Pulse solution found by integrating the (non)autonomous system")



    
    
    [C.Euminus.frame, C.Euminus.timeVec] = C.makeFrame(basis_1, basis_2);
   
end    