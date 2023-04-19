function C = generateEuFrame(C)
    L = C.conjPts.L; 
    orig_ic = C.Euminus.pulse_ic';
    deriv_ic = C.Euminus.pulse_deriv_ic';
    params = C.vfParams; 

    [vectors, values]= getBinfEigs(C);
    v1=real(vectors.u(:,1));
    v2 = imag(vectors.u(:,1));
    
    v1 = v1./vecnorm(v1);
    v2 = v2./vecnorm(v2);

    if C.Euminus.normalize == 1 
        deriv_ic = deriv_ic./vecnorm(deriv_ic);
    end
 
    vec_ic = [v1,v2];

    new_ic1 = [orig_ic; vec_ic(:,1)];
    new_ic2 = [orig_ic; deriv_ic]; 

    [full_phi1, basis_1] = C.integrateDE(new_ic1,-L,L);
    [full_phi1, basis_2] = C.integrateDE(new_ic2,-L,L);
    
    [C.Euminus.frame, C.Euminus.timeVec] = C.makeFrame(basis_1, basis_2);
   
end    