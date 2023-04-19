% This function solves the ODE for the pulse and the basis solution for the
% unstable subspace

function [phisol, vsol] =integrateDE(C, ic,start,finish)  
    options=odeset('MaxStep',0.001, 'RelTol', 1e-10);

    if C.Euminus.normalize == 1
        [t,P]=ode45(@(t,P)C.nonautonODENormalized(t,P),[start finish],ic, options);
    else
        [t,P]=ode45(@(t,P)C.nonautonODE(t,P),[start finish],ic, options);
    end

    
    phisol=[t,P(:,1:4)];
    vsol=[t,P(:,5:8)];
end
