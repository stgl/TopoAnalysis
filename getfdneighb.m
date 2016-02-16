function [i_out, j_out, flag] = getfdneighb(i, j, dims, fd)
    
    %% find index of downstream neighbor based on flow direction (ESRI standard)

    flag = false; 
    i_out = nan; 
    j_out = nan;
    
    if(fd == 1)
        i_out = i;  
        j_out = j + 1;
        flag = ingrid(i_out, j_out, dims);
    elseif(fd == 2)
        i_out = i + 1;
        j_out = j + 1;
        flag = ingrid(i_out, j_out, dims);
    elseif(fd == 4)
        i_out = i + 1;
        j_out = j;
        flag = ingrid(i_out, j_out, dims);
    elseif(fd == 8)
        i_out = i + 1;
        j_out = j - 1;
        flag = ingrid(i_out, j_out, dims);
    elseif(fd == 16)
        i_out = i;
        j_out = j - 1;
        flag = ingrid(i_out, j_out, dims);
    elseif(fd == 32)
        i_out = i - 1;
        j_out = j - 1;
        flag = ingrid(i_out, j_out, dims);
    elseif(fd == 64)
        i_out = i - 1;
        j_out = j;
        flag = ingrid(i_out, j_out, dims);
    elseif(fd == 128)
        i_out = i - 1;
        j_out = j + 1;
        flag = ingrid(i_out, j_out, dims);
    end

    % internal function
    % -------------------------------------------------------------------------

    function flag = ingrid(i, j, dims)
        flag = true;
        
        if(((i < 0) | (i >= dims(1))) | ((j < 0) | (j >= dims(2))))
            flag = false;
        end

    end

end
