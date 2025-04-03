% Read the image file
img = imread('Untitled.jpg');

% Convert to grayscale if the image is RGB
if size(img, 3) == 3
    img = rgb2gray(img);
end

% Normalize the image to a range of 0 to 1
img = double(img) / 255;

% Use the fractal_dimension function
D = fractal_dimension(img, true)
disp(['Fractal Dimension: ', num2str(D)])




function D = fractal_dimension(A, show)
    if nargin < 2
        show = false; % Default value for show
    end
    
    % Ensure A is binary (crystal presence or absence)
    threshold = 0.1; % Adjust as needed for structure detection
    A_binary = A > threshold;
    
    % Get size of the grid
    N = size(A, 1);
    
    % Define box sizes (powers of 2 for efficiency)
    box_sizes = 2.^(0:floor(log2(N))); 
    box_sizes = box_sizes(box_sizes <= N);
    
    % Initialize box count storage
    num_boxes = zeros(size(box_sizes));
    
    % Loop over different box sizes
    for i = 1:length(box_sizes)
        s = box_sizes(i);
        num_boxes(i) = sum(sum(blockproc(A_binary, [s, s], @(block) any(block.data(:)))));
    end
    
    % Remove zero counts for log-log fit
    valid = num_boxes > 0;
    log_sizes = log(1 ./ box_sizes(valid));
    log_counts = log(num_boxes(valid));
    
    % Perform linear regression on log-log data
    p = polyfit(log_sizes, log_counts, 1);
    D = p(1); % The slope gives the fractal dimension
    
    % Display results if show is true
    if show
        figure;
        plot(log_sizes, log_counts, 'o-', 'LineWidth', 2);
        xlabel('log(1/Box Size)');
        ylabel('log(Number of Boxes)');
        title(['Estimated Fractal Dimension: ', num2str(D, '%.3f')]);
        grid on;
    end
end
