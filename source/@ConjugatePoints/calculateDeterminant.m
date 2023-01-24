function C = calculateDeterminant(C)
        all_bases = C.Euminus.frame; 
        lambda = 0; 

        timevec = C.Euminus.timeVec;

        vals=cell(1,3);
        vals{1,2}=timevec;
        vals{1,1} = lambda;

        M=max(size(timevec));        
        det_vec=zeros(1,M);
        
        for j = 1:M
            frame_snapshot=all_bases(:,:,j);
            A1_snapshot=[frame_snapshot(1,:);frame_snapshot(4,:)];
            if isa(all_bases(1,1,1), 'intval')
                det_vec(j) = A1_snapshot(1,1)*A1_snapshot(2,2)...
                    - A1_snapshot(1,2)*A1_snapshot(2,1);
            else
            det_vec(j)=det(A1_snapshot);
            end
        end 
        vals{1,3}=det_vec;
    C.conjPts.dets = vals; 
end