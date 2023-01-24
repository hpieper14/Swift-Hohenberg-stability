function S = trimNFSol_halfline(S)
    old_sol = S.normalForm.sol; 
    old_time = S.normalForm.time; 


    index = find(diff(sign(old_sol(:,2)))); 
    sol = old_sol(1:end-index,:);
    time = old_time(1:end-index);

    S.normalForm.sol = sol; 
    S.normalForm.time  = time; 
    S.time = time(end);
end 