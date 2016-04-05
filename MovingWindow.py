import Error
import numpy as np

class MovingWindow(object):
    
    function = None
    
    def __init__(self, *args, **kwargs):
        if kwargs.get('window_dimension') is None:
            raise Error.InputError('Window dimension','is a required parameter')
        self.window_dimension = kwargs.get('window_dimension')
        
        if self.__class__ is MovingWindow:
            raise Error.Error('MovingWindow is an abstract base class')
    
    def _build_search_kernel(self, dx):        
        return None, None
    
    def __adjust_kernel(self, row,col,grid,searchKernelRows,searchKernelCols):
        # Function to turn the rows and columns of the search kernel specified by searchKernelRows & ...Cols as just indices of the
        # data that is actually within the grid
    
        #Get the current absolute positions within the grid
        theseRows = searchKernelRows+row
        theseCols = searchKernelCols+col
    
        #Do any points in searchKernel extend to before the start of the grid? or off the end?
        gdIndices = (theseCols < grid.shape[1]) & (theseRows < grid.shape[0]) & (theseRows >= 0) & (theseCols >= 0)
    
        #Which rows are within the grid
        goodRows = theseRows[gdIndices]
        goodCols = theseCols[gdIndices]
    
        #Now that we are all within the grid, of those points within the grid, do any include
        #Do any points in searchKernel include nans ?
        gdIndices = np.logical_not(np.isnan(grid[goodRows, goodCols]))
    
        return goodRows[gdIndices], goodCols[gdIndices]

    def apply_moving_window(self, grid):
        # Function to scan the moving window specified by Kernel across the dem, a numpy grid, and apply the specified function
        # to each location. function is a method that returns a single value given any number of inputs, e.g.
        # movingWindow(i) = function(dem[Kernel@i])
    
        outgrid = np.zeros_like(grid)
    
        search_kernel_rows, search_kernel_cols = self._build_search_kernel(grid.dx())
        
        for i in range(grid.shape[0]):
            for j in range(grid.shape[1]):
                these_rows,these_cols = self.__adjustKernel(i,j,grid,search_kernel_rows,search_kernel_cols) #trim the searchKernel and return it as a list of indices
                outgrid[i,j] = self.function(grid[these_rows, these_cols]) #apply the specified function to the dem values for this search location
    
        return outgrid

class RectangularMovingWindow(MovingWindow):
    
    def __init__(self, *args, **kwargs):
        super(RectangularMovingWindow, self).__init__(*args, **kwargs)
        if self.__class__ is RectangularMovingWindow:
            raise Error.Error('RectangularMovingWindow has no bound function')
        
    def _build_search_kernel(self, dx):
        # Function to build a search kernel that is either square or circular
    
        pxlRadius = int(round(self.window_radius/dx)) #Is this right...
        relCoords = np.arange(1 + 2*pxlRadius)-pxlRadius #Relative coordinates of neighbors within this row, col distance
        searchKernelRow, serchKernelCol = np.meshgrid(relCoords,relCoords)
        return searchKernelRow.flatten(), serchKernelCol.flatten()
    
class CircularMovingWindow(MovingWindow):
    
    def __init__(self, *args, **kwargs):
        super(CircularMovingWindow, self).__init__(*args, **kwargs)
        if self.__class__ is CircularMovingWindow:
            raise Error.Error('CircularMovingWindow has no bound function')
        
    def _build_search_kernel(self, dx):
        # Function to build a search kernel that is either square or circular
    
        pxlRadius = int(round(self.window_radius/dx)) #Is this right...
        relCoords = np.arange(1 + 2*pxlRadius)-pxlRadius #Relative coordinates of neighbors within this row, col distance
        searchKernelRow, serchKernelCol = np.meshgrid(relCoords,relCoords)
        dists = np.sqrt(searchKernelRow**2 + searchKernelRow**2)
        searchKernelRow = searchKernelRow[dists<pxlRadius]
        serchKernelCol = serchKernelCol[dists<pxlRadius]
    
        return searchKernelRow.flatten(), serchKernelCol.flatten()