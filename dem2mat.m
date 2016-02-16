function dem = dem2mat(filename)

%% Load ESRI ASCII DEM as MATLAB structure
%% Robert Sare 2014
%%
%% INPUT:       filename - ASCII file containing DEM
%%
%% OUTPUT:      dem - dem in matlab grid struct

% Define DEM structure
dem = struct('nx', 0, 'ny', 0, 'xllcenter', 0, 'yllcenter', 0, 'de', 0, 'grid', [], 'nodata', NaN);

% Load DEM info
fid = fopen(filename);
C = textscan(fid,'%s %f',6);

% Get DEM dimension and resolution
dem.nx = C{2}(1);
dem.ny = C{2}(2);
dem.de = C{2}(5);

% Get DEM corner in UTM coordinates
dem.xllcenter = C{2}(3);
%dem.xllcenter = x1+(dem.nx/2)*dem.de;
dem.yllcenter = C{2}(4);
%dem.yllcenter = y1+(dem.ny/2)*dem.de;

% Load DEM grid
ndv = C{2}(6)
dem.grid = fscanf(fid,'%f',[dem.nx,dem.ny]);
dem.grid = (dem.grid')
%dem.nodata = ndv; % retain NDVs; useful for ESRI compatibility
dem.grid(dem.grid==ndv) = NaN; % set NaNs; useful for matlab processing

end
