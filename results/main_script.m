% This script generates all results in Section 5 of the paper "The Maslov 
% Index, degenerate crossings, and the stability of pulse solutions to the 
% Swift-Hohenberg equation"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
% normal form branch -- set to 0, pi
normalForm.branch = 0; 

% vector field parameters
vfParams.nu = 1.6;
vfParams.mu = .2;  % Set to 0.2 for the stable pulse
vfParams.lambda = 0; 

% fourier approximation parameters
fourier.M = 1000; 
fourier.tol = 1e-14; 
fourier.order = 404; 

% fourier.M = 1000; 
% fourier.tol = 1e-14; 
% fourier.order  ||  Conjugate Points of the stable pulse
% 200   ||   2
% 250   ||  21
% 256   ||   1
% 348   ||   0
% 350   ||  14
% 352   ||   0
% 400   ||   0
% 402   ||  14
% 404   ||   2
% 500   ||   0
% 600   ||   1
% 700   ||   0
 
time = 100; 

% The stable pulse seems to be fairly sensitive to the integration time and
% the number of modes, in that it produces spurrious conjugate points ... 
% Things being finicky with adjusting time seems to make sense, as Newton
% seems to converge to a different profile. 
% However I am not sure why changing the number of Fourier modes is so sensitive. I find this a little concerning/confusing 



% Also, when we use the Newton's method to find the solution in Fourier,
% there appears to be a zero eigenvalue, which causes Newton's method to be
% ill conditioned. The zero-eigenvalue is likely due to a translation
% symmetry. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                PULSE SOLUTION APPROXIMATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% initialize PulseSolution
S = PulseSolution(fourier, vfParams, normalForm, time);

% perform Newton's method
S = S.mainPulse();


% Plot normal form approximation and pulse approximation obtained by
% performing Newton's method 

time = S.normalForm.time; 
sol = S.normalForm.sol(:,1);
full_sol = S.getFunctionFromFourierCoeffs(S.fourier.full_coeffs, "full");
figure 
hold on 
if vfParams.mu~=.2
    plot(full_sol(:, 1), full_sol(:, 2), color = 'b', LineWidth=1.25)
    plot([-flip(time), time], [flip(sol), sol], Color  = 'r', LineWidth=1.25)
    legend('Solution after Newtons method, $\bar \varphi^{(N)}$', 'Solution via Normal Form Eqn, $u_\phi$', Interpreter = 'latex')
    title('Pulse Approximations, $\phi =$ ' + string(normalForm.branch), Interpreter = 'latex')
else
    plot(full_sol(:, 1), full_sol(:, 2), color = 'b', LineWidth=1.25)
    legend('Solution after Newtons method, $\bar \varphi^{(N)}$', Interpreter = 'latex')
    title('Pulse Approximations, snaking', Interpreter = 'latex')
    hold off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 CONJUGATE POINT COMPUTATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set parameters
conjPts.L = 55; 
C = ConjugatePoints(conjPts, vfParams);
C.Euminus.normalize = 1;
C.Euminus.refPlane = [1,4]; 

% compute conjugate points
[S,C] = C.mainConjPts(S);




% Plot basis vectors for $E^u_-(x;0)$ 
time = C.Euminus.timeVec;
basis1 = reshape(C.Euminus.frame(:, 1,:), [4, max(size(time))]);
figure
tiledlayout(4,1)
nexttile
plot(time, basis1(1,:))
nexttile
plot(time, basis1(2,:))
nexttile
plot(time, basis1(3,:))
nexttile
plot(time, basis1(4,:))
title('First Basis Vector for $E^u_-(x; 0)$', Interpreter = 'latex')
xlabel('$x$', Interpreter = 'latex')

basis2 = reshape(C.Euminus.frame(:, 2,:), [4, max(size(time))]);
figure
tiledlayout(4,1)
nexttile
plot(time, basis2(1,:))
nexttile
plot(time, basis2(2,:))
nexttile
plot(time, basis2(3,:))
nexttile
plot(time, basis2(4,:))
title('Second Basis Vector (Derivative of the Pulse) for $E^u_-(x; 0)$', Interpreter='latex')
xlabel('$x$', Interpreter = 'latex')


% plot determinant function 

dets = C.conjPts.dets{:, 3}; 
time = C.conjPts.dets{:,2}; 

figure 
plot(time, dets, LineWidth=1.25)
if vfParams.mu~=.2
    title('Determinant $\det(A)$, $\phi =$ ' + string(normalForm.branch), Interpreter='latex')
else
    title('Determinant $\det(A)$, snaking ' , Interpreter='latex')
end
xlabel('$x$', Interpreter = 'latex')

% print eigenvalues and conjugate point locations 

disp('Eigenvalue approximation via Fourier modes: ')
disp(S.fourier.unstable_eigs)

conj_ind = find([0, diff(sign(dets))] ~= 0);
conj_pts = time(conj_ind);

disp('Conjugate points locations (x): ')
disp(conj_pts)


disp("We found " + string(max(size(S.fourier.unstable_eigs))) + " eigenvalue(s) and "...
    + string(length(conj_pts)) + " conjugate points for the pulse with phase " ...
    +'phi = ' + string(normalForm.branch))


