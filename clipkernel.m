function kern = clipkernel(i, j, dims, radius, de)
    
    %% Return grid values at indices within a square search kernel of specified
    %% radius

    idx = sub2ind(dims, i, j);
    
    ncells = floor(radius/de);
    row = i + [-ncells:ncells];
    col = j + [-ncells:ncells];
    
    [I, J] = meshgrid(row, col);
    mask = find(I > 0 & J > 0 & I <= dims(1) & J <= dims(2)); % avoid corner cases

    neighb = sub2ind(dims, I(mask), J(mask));
    neighb = neighb(neighb ~= idx);
    kern = [idx neighb'];

end
