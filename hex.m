clear; clc; close all;

iterations = 200;
beta = 0.6;             % Value for background level
gamma = 0.0035;        % Constant added to receptive sites
alpha = 1;
crystal_size = 151;     % Number of hex on top of center (note only use odd numbers)
mode = 2
if (mode ==1)
    %Simple test 
    tic  % Start timing

    A = initialize_crystal(crystal_size, beta);
    [X, Y, Dis] = generate_grid(crystal_size);

    A = compute_A(A, iterations, beta, gamma, alpha);
    
    toc

    tic
    clf

    drawing(A,X,Y,Dis,crystal_size)
    drawnow

    toc

    tic
    disp(fractal_dimension(A))
    toc
end
if (mode == 2)
    % Évaluation de la convergence de la dimension de corrélation, en fonction du nombre d'itérations
    N_iter = [50, 100, 150, 175,200,250,300,350,400,450,500]; % Number of iterations to evaluate
    D_iter = zeros(1, length(N_iter));
    
    for i = 1:length(N_iter)
        iterations = N_iter(i);
        A = initialize_crystal(crystal_size, beta);
        [X, Y, Dis] = generate_grid(crystal_size);
        
        A = compute_A(A, iterations, beta, gamma, alpha);
        
        D_iter(i) = fractal_dimension(A);
    end
    
    % Plot the results
    figure;
    plot(N_iter, D_iter, 'o-', 'LineWidth', 2);
    xlabel('Number of Iterations');
    ylabel('Fractal Dimension');
    title('Convergence of Fractal Dimension with Iterations');
    grid on;
end
if (mode == 4)
    %Évaluation de la convergence de la dimension de corrélation, en fonction de la taille du crystal

end
if (mode == 5)
    % Évaluation de la dimension fractale en fonction de beta
    tic
    N = 20; % Number of beta values to evaluate
    D = zeros(1, N); % Array to store fractal dimensions
    betas = linspace(0, 1, N); % Range of beta values
    
    for i = 1:N
        beta = betas(i); % Current beta value
        A = initialize_crystal(crystal_size, beta); % Initialize crystal
        [X, Y, Dis] = generate_grid(crystal_size); % Generate grid
        
        A = compute_A(A, iterations, beta, gamma, alpha); % Compute crystal growth
        
        D(i) = fractal_dimension(A); % Calculate fractal dimension
    end
    
    % Plot fractal dimension as a function of beta
    figure;
    plot(betas, D, 'o-', 'LineWidth', 2);
    xlabel('\beta (Background Level)');
    ylabel('Fractal Dimension');
    title('Fractal Dimension vs \beta');
    grid on;
    toc
end
% ---------------------------------------------------------------
% Functions
% ---------------------------------------------------------------

function A = initialize_crystal(crystal_size, beta)
    A = beta * ones(crystal_size*2 +1, crystal_size*2 + 1 ); % Array to store states
    A(crystal_size+1, crystal_size+1) = 1;  % Setting the center to ice
end

function [X, Y, Dis] = generate_grid(crystal_size)
    x = -crystal_size:1:crystal_size;  % Array to store x positions of hexagons
    y = crystal_size:-1:-crystal_size; % Array to store y positions of hexagons
    [X, Y] = meshgrid(x,y);

    % Add offset to alternate rows and columns
    for i = 1:2:crystal_size*2 +1
        X(i,:) = X(i,:) + 0.5;
    end

    Dis = sqrt(X.^2 + Y.^2);
end

function A = compute_A(A, iterations, beta, gamma, alpha)
    for i = 1 : iterations
        [odd_row, even_row] = neighbor_array(A, beta);

        receptive_odd = sum((odd_row  >= 1) , 3)>0;
        receptive_even = sum((even_row  >= 1) , 3)>0;
        receptive_logic = zeros(size(A));
        for k = 1:size(A,1)
            if mod(k,2) == 0
                receptive_logic(k,:) = receptive_even(k,:);
            else
                receptive_logic(k,:) = receptive_odd(k,:);
            end
        end
        unreceptive_logic = receptive_logic==0;
        receptive = receptive_logic .*  A;
        unreceptive = unreceptive_logic .*  A;

        receptive_update = (receptive_logic * gamma) + receptive;
        [odd_row, even_row] = neighbor_array(unreceptive, beta);

        unreceptive_ne = zeros(size(odd_row));
        for k = 1:size(A,1)
            if mod(k,2) == 0
                unreceptive_ne(k,:,:) = even_row(k,:,:);
            else
                unreceptive_ne(k,:,:) = odd_row(k,:,:);
            end
        end

        averaged_neighbors = ((unreceptive_ne(:,:,1)*0.5) + (unreceptive_ne(:,:,2)*(1/12)) + (unreceptive_ne(:,:,3)*(1/12)) ...
            + (unreceptive_ne(:,:,4)*(1/12)) + (unreceptive_ne(:,:,5)*(1/12)) + (unreceptive_ne(:,:,6)*(1/12)) + (unreceptive_ne(:,:,7)*(1/12)));

        final_update = averaged_neighbors + alpha * receptive_update;
        
        A = final_update;
    end
end



% ---------------------------------------------------------------
% Functions
% ---------------------------------------------------------------




function D = fractal_dimension(A, show)
    if nargin < 2
        show = false; % Default value for show
    end
    
    % Ensure A is binary (crystal presence or absence)
    threshold = 1; % Adjust as needed for structure detection
    A_binary = A < threshold;
    
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



function drawing(A, X, Y, Dis, crystal_size)

    hold on; % Keep hold for efficient plotting
    
    % Filter out hexagons outside the desired range
    validIdx = Dis < (crystal_size + 0.1);
    X_valid = X(validIdx);
    Y_valid = Y(validIdx);
    A_valid = A(validIdx);
    
    % Precompute hexagon shape (one-time calculation)
    hexCorners = hexagon([0, 0], 0.5);
    
    % Number of hexagons
    numHex = numel(X_valid);
    
    % Generate coordinates for all hexagons efficiently
    X_patches = repmat(hexCorners(1, :), numHex, 1) + X_valid(:);
    Y_patches = repmat(hexCorners(2, :), numHex, 1) + Y_valid(:);
    
    % Compute opacities with better contrast range
    opacities = 1 - min(max(A_valid / 2, 0), 0.5); % Ensures values are within range [0.5, 1]
    
    % Plot hexagons efficiently with better visibility
    patch(X_patches', Y_patches', 'k', ...
        'FaceVertexAlphaData', opacities, 'FaceAlpha', 'flat', ...
        'EdgeColor', 'w', 'LineWidth', 0.5); % White edges for visibility
    
    axis equal; % Keep aspect ratio correct
    hold off;
    
    end
    
    







function corners = hexagon(center, radius)

radius = radius / cosd(30);
angles = [30 90 150 210 270 330];

x = cosd(angles)*radius + center(1);
y = sind(angles)*radius + center(2);
corners = [x;y];

end







function [odd_row, even_row] = neighbor_array(A, beta)

top_left = [beta*ones(1,size(A,2)+1) ; [beta*ones(size(A,1),1) A]];
top_left = top_left(1:end-1 , 1:end-1);

top = [beta*ones(1, size(A,2)); A];
top = top(1:end-1 , :);

top_right = [beta*ones(1,size(A,2)+1) ; [A beta*ones(size(A,1),1)]];
top_right = top_right(1:end-1 , 2:end);

left = [beta*ones(size(A,1), 1) A];
left = left(: , 1:end-1);

right = [A beta*ones(size(A,1), 1)];
right = right(: , 2:end);

bottom_left = [[beta*ones(size(A,1),1) A] ; beta*ones(1,size(A,2)+1)];
bottom_left = bottom_left(2:end , 1:end-1);

bottom = [A ; beta*ones(1, size(A,2))];
bottom = bottom(2:end , :);

bottom_right = [[A beta*ones(size(A,1),1)] ; beta*ones(1,size(A,2)+1)];
bottom_right = bottom_right(2:end , 2:end);

odd_row = cat(3, A, top, top_right, left, right, bottom, bottom_right);
even_row = cat(3, A, top_left, top, left, right, bottom_left, bottom);

end