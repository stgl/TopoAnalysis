function outgrid = window(ingrid, radius, functionname)
   
    %% Compute moving average of parameter on grid using window radius
    %% specified in grid resolution units
    %% 
    %% Input function takes central value as first argument and vector of values
    %% in search kernel as subsequent arguments

    [m,n] = size(ingrid.grid);
    outgrid = ingrid;
    outgrid.grid = nan*ones(m,n);

    for(i=1:m)
        for(j=1:n)
            neighb = clipkernel(i, j, [m,n], radius, ingrid.de);
            outgrid.grid(i,j) = feval(functionname, ingrid.grid(i,j), [ingrid.grid(neighb)]);
        end
    end

end

