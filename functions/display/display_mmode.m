function img = display_mmode(ax, mmode_matrix, mmode_imagethresh, mmode_x, mmode_y)
    % Convert the matrix to a normalized grayscale image using the given threshold
    img = mat2gray(mmode_matrix, [0 mmode_imagethresh]);

    % some image processing
    img = adapthisteq(img,'clipLimit',0.02,'Distribution','rayleigh');
    
    % Display the image in the specified axes
    imagesc(ax, mmode_x, mmode_y, img);
    
    % Set the colormap to grayscale so that higher values are brighter
    % colormap(ax, gray);
    
    % Adjust axes properties for proper display of the image
    % set(ax, 'YDir', 'normal'); % Optionally set y-axis direction to normal
    
    % make the image fit the axis
    axis(ax, 'tight');

    %remove this if you dont want depth limit
    ylim(ax, [0 30]);
end
