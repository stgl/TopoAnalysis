function out = calcaccumarea(fd, outfile)

    %% Compute accumulation area (km^2) based on grid of flow directions

    area = calcareagrid(fd);
    A = calcaccgrid(fd, area);

    out = fd;
    out.grid = A;

    mat2dem(out, outfile);

end
