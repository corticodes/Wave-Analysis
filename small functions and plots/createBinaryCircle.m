function ringPixels = createBinaryCircle(imageSizeX,imageSizeY,centerX,centerY,radius,thickness,plotCircle)
%CREATEBINARYCIRCLE creates logical imageSizeY by imageSizeX matrix with
%ones on the circumfrance of a circle located at centerX,centerY 
%   plotCircle is an option to also plot the circle
 
[columnsInImage, rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);
% Next create the circle in the image.

outerRadius = radius+thickness;
array2D = (rowsInImage - centerY).^2 ...
    + (columnsInImage - centerX).^2;
ringPixels = array2D >= radius.^2 & array2D <= outerRadius.^2;
% ringPixels is a 2D "logical" array.
% Now, display it.
if plotCircle
    image(ringPixels) ;
    colormap([0 0 0; 1 1 1]);
    title('Binary Image of a Ring', 'FontSize', 25);
end
end

