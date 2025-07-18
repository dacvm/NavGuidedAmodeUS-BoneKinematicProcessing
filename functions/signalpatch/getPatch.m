function patch = getPatch(I, init, window_depth)

    row_start = round(init(1)-(window_depth/2));
    row_end   = round(init(1)+(window_depth/2));
    col_start = init(2);
    col_end   = init(2) + 10;

    patch = I(row_start:row_end, col_start:col_end);
end

