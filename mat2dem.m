function [] = mat2dem(dem, filename)

%% Save MATLAB DEM as ESRI ASCII grid
%% Robert Sare 2014
%%
%% INPUT:       dem - dem to save
%%              filename - ASCII file 

% Load DEM info
fid = fopen(filename, 'w');

% Write grid parameters to file
fprintf(fid, 'ncols %i\n', dem.nx);
fprintf(fid, 'nrows %i\n', dem.ny);
fprintf(fid, 'xllcenter %10.3f\n', dem.xllcenter);
fprintf(fid, 'yllcenter %10.3f\n', dem.yllcenter);
fprintf(fid, 'cellsize %8.3f\n', dem.de);
fprintf(fid, 'nodata_value %10.3f\n', dem.nodata);

idx = find(isnan(dem.grid));
dem.grid(idx) = dem.nodata;

for(i = 1:1:dem.ny)
    fprintf(fid, '%10.3f', dem.grid(i,:));
    fprintf(fid, '\n');
end

fclose(fid);

end
