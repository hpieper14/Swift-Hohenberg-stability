vfParams.nu = 1.6;
vfParams.mu = .05;
vfParams.lambda = 0; 

normalForm.branch = pi; 
fourier.M = 1000; 
fourier.tol = 1e-14; 
fourier.order = 150; 
 
time = 100; 


S = PulseSolution(fourier, vfParams, normalForm, time);

S = S.mainPulse();


time = S.normalForm.time; 
sol = S.normalForm.sol(:,1);
full_sol = S.getFunctionFromFourierCoeffs(S.fourier.full_coeffs, "full");
figure 
hold on 
plot(full_sol(:, 1), full_sol(:, 2), color = 'b', LineWidth=1.25)
plot([-flip(time), time], [flip(sol), sol], Color  = 'r', LineWidth=1.25)
legend('Solution after Newtons method, $\bar \varphi^{(N)}$', 'Solution via Normal Form Eqn, $u_\phi$', Interpreter = 'latex')
title('Pulse Approximations, $\phi = \pi$', Interpreter = 'latex')



conjPts.L = 55; 
C = ConjugatePoints(conjPts, vfParams);
C.Euminus.normalize = 1;
C.Euminus.refPlane = [1,4]; 

[S,C] = C.mainConjPts(S);



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
title('First basis vector')

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
title('Second basis vector')



dets = C.conjPts.dets{:, 3}; 
time = C.conjPts.dets{:,2}; 

figure 
plot(time, dets, LineWidth=1.25)
title('Determinant $\det(A)$, $\phi =$ ' + string(normalForm.branch), Interpreter='latex')
xlabel('$x$', Interpreter = 'latex')


