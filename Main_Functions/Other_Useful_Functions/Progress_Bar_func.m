function Progress_Bar_func(state, varargin)
    
    persistent progress total_window_count
        
    progress_bar_length = 100;

    if strcmp(state, 'begin')
        if nargin == 2
            total_window_count = varargin{1};
        else
            error('ERROR! You did not enter total window count.');
        end
        progress = 0;
        fprintf('Progress Bar:\n');
        fprintf(['[', repmat('.', 1, progress_bar_length), ']', '\n\n']);
    elseif strcmp(state, 'ongoing')
        progress = progress + 1;
        progress_counter = floor( (progress / total_window_count) * progress_bar_length );
        fprintf([repmat('\b', 1, progress_bar_length + 4), '[', repmat('#', 1, progress_counter), repmat('.', 1, progress_bar_length - progress_counter), ']', '\n\n']);
    elseif strcmp(state, 'ended')
        progress = 0;
        total_window_count = 0;
        fprintf([repmat('\b', 1, progress_bar_length + 4), '[', repmat('#', 1, progress_bar_length), ']', '\n\n']);
    end

end
