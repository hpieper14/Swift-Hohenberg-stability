function test = testBKNormalForm4dim()
    order = 10;
    params.nu = 1.2;
    params.mu = .1;
    vfParams = params;
    nfBranch = 0; 
    time = 10; 

    S = PulseSolution(order, vfParams, nfBranch, time);
    
    % get data for normal form solution
    S = BKNormalForm4dim(S); 

    time = S.nfData.time;
    sol = S.nfData.sol;

    figure 
    tiledlayout(4,1)
    nexttile
    plot(time, sol(:, 1))
    nexttile 
    plot(time, sol(:, 2))
    nexttile 
    plot(time, sol(:, 3))
    nexttile 
    plot(time, sol(:, 4))
    title("Normal Form Solution")
    
    % trim time domain so normal form solution satisfies Neumann BC 
    S = trimNFSol(S); 

    time = S.nfData.time;
    sol = S.nfData.sol;

    figure 
    tiledlayout(4,1)
    nexttile
    plot(time, sol(:, 1))
    nexttile 
    plot(time, sol(:, 2))
    nexttile 
    plot(time, sol(:, 3))
    nexttile 
    plot(time, sol(:, 4))
    title("Normal Form Solution with Neumann BC")

    test = 0; 
end