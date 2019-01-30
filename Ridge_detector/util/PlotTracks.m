function PlotTracks(ax, tracks)
    % Plot reported 'tracks' on the given axis.
    % 'tracks' variable is a cell array, and each cell is an Nx3 matrix
    % in which rows 1, 2 & 3 contain a TF contour's time, frequency and
    % spectral intensity (in dB) values, respectively.
    % 'ax' is the second parameter that was passed to SetCallback() above.
    
    % Functionality to roll color selection
    persistent tid;
    if isempty(tid)
        tid = 0;
    end
    line_colors = hsv(8);
    
    num_tracks = length(tracks);

    hold(ax, 'on');
    for track_idx = 1:num_tracks
        plot(ax, tracks{track_idx}(1, :), tracks{track_idx}(2, :), '.-', ...
            'color', line_colors(mod(tid, 8)+1, :));
        tid = tid + 1;  % Update next color idx
    end
    hold(ax, 'off');
    drawnow;
    
    if num_tracks > 0
        fprintf('Added %i tracks to axis\n', num_tracks);
    end
end