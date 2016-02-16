function area = calcareagrid(ingrid)

    %% Compute area of grid cells in km^2

    re = 6371;

    dx = ingrid.de; % change this if grid has different x and y resolutions
    dy = dx;
    x0 = ingrid.xllcenter; % upper left coordinates
    y0 = ingrid.yllcenter + ingrid.ny*dy;

    [m, n] = size(ingrid.grid);
    area = nan*ones(m,n);

    for(i=1:m)
        lat1 = (y0 - dy*(i-1))*(pi/180);
        lat0 = (y0 - dy*i)*(pi/180);
        dlon = dx*(pi/180);
        % rows have equal latitudinal spacing, so compute areas by row
        area(i,:) = (sin(lat1) - sin(lat0))*dlon*re^2;
    end

end
