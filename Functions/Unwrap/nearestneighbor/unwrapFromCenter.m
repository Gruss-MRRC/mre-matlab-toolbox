function ph_out = unwrapFromCenter(ph1,unwrapDir,centPosition)
% unwrap from center, correcting phase jumps

ph_out = ph1;

if unwrapDir == 1,
    for j = 1:size(ph1,2),
        for k = centPosition:size(ph1,unwrapDir),
            jump1 = (ph_out(k,j) - ph_out(k-1,j));
            if abs(jump1) > 2000 && ph1(k,j) > 0 && ph1(k-1,j) > 0,
                ph_out(k,j) = ph_out(k,j) - sign(jump1)*4096;
            end
        end
        for k = centPosition-1:-1:2,
            jump1 = (ph_out(k,j) - ph_out(k+1,j));
            if abs(jump1) > 2000 && ph1(k,j) > 0 && ph1(k+1,j) > 0,
                ph_out(k,j) = ph_out(k,j) - sign(jump1)*4096;
            end
        end
    end
else
    for j = 1:size(ph1,1),
        for k = centPosition:size(ph1,unwrapDir),
            jump1 = (ph_out(j,k) - ph_out(j,k-1));
            if abs(jump1) > 2000 && ph1(j,k) > 0 && ph1(j,k-1) > 0,
                ph_out(j,k) = ph_out(j,k) - sign(jump1)*4096;
            end
        end
        for k = centPosition-1:-1:2,
            jump1 = (ph_out(j,k) - ph_out(j,k+1));
            if abs(jump1) > 2000 && ph1(j,k) > 0 && ph1(j,k+1) > 0,
                ph_out(j,k) = ph_out(j,k) - sign(jump1)*4096;
            end
        end
    end
end