function [F,J]=fourierODE_halfline(S, u)

Lx= S.time;                % Domain trunction
% Differentiation matrices
% t = theta coordinate, ranging from 0 to 2*pi (M must be even)
  M = S.fourier.M;
  dt = 2*pi/M; t = dt*(1:M)'; M2 = M/2; % Fourier mesh
  
  column = [0 .5*(-1).^(1:M-1).*cot((1:M-1)*dt/2)];
  Dt = toeplitz(column,column([1 M:-1:2]));
  
  D2t = toeplitz([-pi^2/(3*dt^2)-1/6 ...
                 .5*(-1).^(2:M)./sin(dt*(1:M-1)/2).^2]);
  
  % rewrite matrix for 0..pi reflect
  semiDt = zeros(M2+1); 
  semiDt(:,1) = Dt(M2:M,M2);
  semiDt(:,2:M2) = Dt(M2:M,M2-1:-1:1)+Dt(M2:M,M2+1:M-1);
  semiDt(:,M2+1)=Dt(M2:M,M);
  
  semiD2t = zeros(M2+1);
  semiD2t(:,1) = D2t(M2:M,M2);
  semiD2t(:,2:M2) = D2t(M2:M,M2-1:-1:1)+D2t(M2:M,M2+1:M-1);
  semiD2t(:,M2+1)=D2t(M2:M,M);

  Dx  = (pi/S.time)*semiDt;         % Differentiation matrices on half line 
  D2x  = (pi/S.time)^2*semiD2t;
  x=S.time*(t-pi)/pi;               % Mesh on full line.

N = M2+1;
I = eye(N);
LN = -D2x^2 - 2*D2x - I; 

mu = S.vfParams.mu; 
nu = S.vfParams.nu; 


F = LN*u - mu*u + nu*u.^2 - u.^3;

if nargout > 1  
    J = LN + diag(-mu + 2*nu*u - 3*u.^2);
end
