clear all 
close all 

% nonsnaking region
mu = 0.05; 
nu = 1.6; 
generate_results(0, mu, nu)
generate_results(pi, mu, nu)

% snaking region 
% TODO Add comment to readme about scaling for stable mu
mu = .2; 
nu = 1.6; 
generate_results(0, mu, nu)