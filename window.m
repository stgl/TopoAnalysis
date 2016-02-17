function outgrid = window(ingrid, radius, functionname)
   
    %% Compute arbitrary function on grid using a square kernel of radius
    %% specified in grid resolution units
    %%
    %% INPUT:   ingrid - input grid struct
    %%          radius - search radius in grid resolution units
    %%          functionname - function handle for desired calculation
    %%          takes one arguments: v, vector of values in window with central value first
    %%
    %% EXAMPLE: local relief over 100 m radius search kernel 
    %% lrelief = @(v) max(v) - min(v);
    %% outgrid = window(dem, 100, lrelief);
    %% 

    [m,n] = size(ingrid.grid);
    outgrid = ingrid;
    outgrid.grid = nan*ones(m,n);

    for(i=1:m)
        for(j=1:n)
            win = clipkernel(i, j, [m,n], radius, ingrid.de);
            outgrid.grid(i,j) = feval(functionname, win);
        end
    end

end

