function C = getConjPtLocations(C) 
    x_vals = C.conjPts.dets{1,2};
    dets = C.conjPts.dets{1,3};
    conj_point_index = diff(sign(dets));
    conj_point_locs = x_vals(conj_point_index);

    C.conjPts.locs = conj_point_locs; 
end
