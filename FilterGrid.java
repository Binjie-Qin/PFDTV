class FilterGrid {
   public static void nPrintln(int rows, int cols) {
           
}



function [radius, u1, u2] = filtergrid(rows, cols)

    // Handle case where rows, cols has been supplied as a 2-vector
    if nargin == 1 & length(rows) == 2  
        tmp = rows;
        rows = tmp(1);
        cols = tmp(2);
    end
    
    //% Set up X and Y spatial frequency matrices, u1 and u2, with ranges
    //% normalised to +/- 0.5 The following code adjusts things appropriately for
    //% odd and even values of rows and columns so that the 0 frequency point is
    //% placed appropriately.
    if mod(cols,2)
        u1range = [-(cols-1)/2:(cols-1)/2]/(cols-1);
    else
        u1range = [-cols/2:(cols/2-1)]/cols; 
    end
    
    if mod(rows,2)
        u2range = [-(rows-1)/2:(rows-1)/2]/(rows-1);
    else
        u2range = [-rows/2:(rows/2-1)]/rows; 
    end
    
    [u1,u2] = meshgrid(u1range, u2range);
    
    //% Quadrant shift so that filters are constructed with 0 frequency at
    //% the corners
    u1 = ifftshift(u1);
    u2 = ifftshift(u2);
    
    //% Construct spatial frequency values in terms of normalised radius from
    //% center. 
    radius = sqrt(u1.^2 + u2.^2);  
