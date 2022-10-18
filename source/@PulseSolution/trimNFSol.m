function S = trimNFSol(S)
    old_sol = S.nfData.sol; 
    old_time = S.nfData.time; 
    index = find(diff(sign(old_sol(:,2)))); 
    sol = old_sol(index:end-index,:);
    time = old_time(index:end-index);

    S.nfData.sol = sol; 
    S.nfData.time  = time; 
end 