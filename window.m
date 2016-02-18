function outgrid = window(ingrid, radius, functionname)
   
    %% Compute arbitrary function on grid using a square kernel of radius
    %% specified in grid resolution units
    %%
    %% INPUT:   ingrid - input grid struct
    %%          radius - search radius in grid resolution units
    %%          functionname - function handle for desired calculation
    %%          takes one arguments: v, vector of values in window with central value first
    %%
    %% EXAMPLE: local relief over 5 km radius search kernel 
    %% lrelief = @(v) max(v) - min(v);
    %% outgrid = window(dem, 5000, lrelief);
    %% 

    [m,n] = size(ingrid.grid);
    outgrid = ingrid;
    outgrid.grid = nan*ones(m,n);
    ncells = floor(radius/ingrid.de);
    nr = -ncells:1:ncells;
    nc = nr;

    for(i=1:m)
        for(j=1:n)
            ii = i + nr;
            jj = j + nc;
            ii = ii(ii > 0 & ii <= m & ii ~= i);
            jj = jj(jj > 0 & jj <= n & jj ~= j);
            ii = [i ii];
            jj = [j jj];

            win = ingrid.grid(ii, jj);
            outgrid.grid(i,j) = feval(functionname, win(:));
        end
    end

end

