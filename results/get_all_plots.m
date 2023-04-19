function plots = get_all_plots(branch)

vfParams.nu = 1.6;
vfParams.mu = .05;
vfParams.lambda = 0; 

normalForm.branch = branch; 
fourier.M = 1000; 
fourier.tol = 1e-12; 
fourier.order = 150; 
 
time = 100; 


S = PulseSolution(fourier, vfParams, normalForm, time);
S = S.BKNormalForm4d_halfline();
S = S.trimNFSol_halfline();
S = BKNormalForm4d_halfline(S); 


S = S.Newton_halfline(); 


fun = S.getFunctionFromFourierCoeffs(S.fourier.half_coeffs, "half");

time = S.normalForm.time; 
sol = S.normalForm.sol(:,1);
full_sol = S.getFunctionFromFourierCoeffs(S.fourier.full_coeffs, "full");
figure 
hold on 
plot(full_sol(:, 1), full_sol(:, 2), color = 'b', LineWidth=1.25)
plot([-flip(time), time], [flip(sol), sol], Color  = 'r', LineWidth=1.25)
legend('Solution after Newtons method, $\bar \varphi^{(N)}$', 'Solution via Normal Form Eqn, $u_\phi$', Interpreter = 'latex')
title('Pulse Approximations, $\phi = \pi$', Interpreter = 'latex')




S = S.Newton();

DF=S.DFFourier(S.fourier.full_coeffs);

D=eig(DF);
[i,j] = find(D>0);
vals = [];
if min(size(D(i,j))) == 1
    disp(D(i,j));
    vals = [vals, D(i,j)];
else 
    A = D(i,j);
    disp(A(:,1))
    vals = [vals, A(:, 1)];
end

S.fourier.unstable_eigs = vals;

conjPts.L = 60; 

C = ConjugatePoints(conjPts, vfParams); 

C = C.get4DimIC(S);

[vectors, values]= C.getBinfEigs();
isLagrangian = C.verifyLagrangian(vectors.u); 

C.Euminus.normalize = 1;

C = C.generateEuFrame();

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
title('First normalized basis vector')

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
title('Second normalized basis vector')




C = C.calculateDeterminant();

dets = C.conjPts.dets{:, 3}; 
time = C.conjPts.dets{:,2}; 

figure 
plot(time, dets, LineWidth=1.25)
title('Determinant normalized $\det(A)$, $\phi = \pi$', Interpreter='latex')
xlabel('$x$', Interpreter = 'latex')

%%%%%%%%%%%%%%%%%%%%%

C.Euminus.normalize = 1;
C.Esplus.normalize = 1; 

C = C.generateEuFrame();

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
title('First unnormalized basis vector')

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
title('Second unnormalized basis vector')




C = C.calculateDeterminant();

dets = C.conjPts.dets{:, 3}; 
time = C.conjPts.dets{:,2}; 

figure 
plot(time, dets, LineWidth=1.25)
title('Determinant unnormalized $\det(A)$, $\phi = \pi$', Interpreter='latex')
xlabel('$x$', Interpreter = 'latex')

plots = C;

end
