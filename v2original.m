iterations = 5000;
%beta = 0.95;
%gamma = 0.0035;
crystal_size = 101;

%A = generate_flake(iterations, beta, gamma, crystal_size);
%drawing(A, beta);
%fractal_dimension(A)
mode = 0;

if mode == 0
    tic
    seuil =200;
    beta = 0.35;
    gamma1 = 0.0001;
    gamma2 = 0.0005;
    A = generate_flake_radial(iterations, beta, gamma1,gamma2, crystal_size);
    drawing(A, beta);
    %fractal_dimension(A, true);
    saveas(gcf, sprintf('gamma%dgamma%d.png', gamma1, gamma2));

    toc
end

if mode ==3
    tic
    seuil =500;
    beta = 0.8;
    gamma1 = 0.0005;
    gamma2 = 0.001;
    A = generate_flake_seuil(iterations, beta, gamma1,gamma2,seuil, crystal_size);
    drawing(A, beta);
    %fractal_dimension(A, true);
    saveas(gcf, sprintf('seuil%d.png', iterations));

    toc
end

if mode == 1
    % Mode 1: Generate and display the flake
    tic
    beta = 0.3;
    gamma = 0.0001;
    A = generate_flake(iterations, beta, gamma, crystal_size);
    drawing(A, beta);
    %fractal_dimension(A, true);
    saveas(gcf, sprintf('%d.png', iterations));

    toc
end
if mode ==2

betas = linspace(0.3, 0.90, 10);
gammas = linspace(0, 0.05, 50); % Increase resolution for gamma
colors = jet(length(betas)); % Gradient of colors for each beta
figure;
hold on;

D_values = zeros(length(betas), length(gammas)); % Store fractal dimensions for each beta & gamma

for i = 1:length(betas)
    beta = betas(i);
    
    for j = 1:length(gammas)
        gamma = gammas(j);
        A = generate_flake(iterations, beta, gamma, crystal_size);
        
        % Compute fractal dimension
        D_values(i, j) = fractal_dimension(A, false);
        
        % Display progress
        progress = ((i-1) * length(gammas) + j) / (length(betas) * length(gammas)) * 100;
        fprintf('Beta: %.2f, Gamma: %.4f, Progress: %.2f%%\n', beta, gamma, progress);
    end
    
    % Plot results for each beta
    plot(gammas, D_values(i, :), 'Color', colors(i, :), 'LineWidth', 2, 'DisplayName', ['Beta = ', num2str(beta)]);
end


figure(1);
for i = 1:length(betas)-2
    beta = betas(i);
    hold on
     Plot results for each beta
    plot(gammas, D_values(i, :), 'Color', [0.5, 0.5, 0.5], 'LineWidth', 2, 'DisplayName', ['Beta = ', num2str(beta)]);
end

%hold off;
xlabel('\gamma');
ylabel('Dimension de recouvrement (D)');
title('Dimension de recouvrement D(\gamma) pour différents \beta');
legend('show');
grid on;
%colormap(jet);
%hold off;
%saveas(gcf,'fractal_dimension_gamma_beta_8.png');
%colorbar('Ticks', linspace(0, 1, length(betas)), 'TickLabels', arrayfun(@num2str, betas, 'UniformOutput', false));
end
function A = generate_flake(iterations, beta, gamma, crystal_size)
    % Initialisation de la matrice d'état
    A = beta * ones(crystal_size*2 + 1, crystal_size*2 + 1);  % Tableau pour stocker les états

    A(crystal_size +1, crystal_size + 1) = 1; % Centre de cristallisation
    %A(crystal_size +8, crystal_size + 1) = 1; % Centre de cristallisation

    
    % Pour simuler la grille hexagonale
    x = -crystal_size:1:crystal_size;
    y = crystal_size:-1:-crystal_size;
    [X, Y] = meshgrid(x, y);

    % Décalage sur les lignes impaires pour simuler l'empilement hexagonal
    for i = 1:2:crystal_size*2 + 1
        X(i,:) = X(i,:) + 0.5;
    end

    Dis = sqrt(X.^2 + Y.^2);

    for i = 1:iterations
        [odd_row, even_row] = neighbor_array(A, beta);

        receptive_odd = sum((odd_row >= 1), 3) > 0;
        receptive_even = sum((even_row >= 1), 3) > 0;
        receptive_logic = zeros(size(A));

        for k = 1:size(A, 1)
            if mod(k, 2) == 0
                receptive_logic(k, :) = receptive_even(k, :);
            else
                receptive_logic(k, :) = receptive_odd(k, :);
            end
        end

        unreceptive_logic = receptive_logic == 0;
        receptive = receptive_logic .* A;
        unreceptive = unreceptive_logic .* A;

        receptive_update = (receptive_logic * gamma) + receptive;
        [odd_row, even_row] = neighbor_array(unreceptive, beta);

        unreceptive_ne = zeros(size(odd_row));
        for k = 1:size(A, 1)
            if mod(k, 2) == 0
                unreceptive_ne(k, :, :) = even_row(k, :, :);
            else
                unreceptive_ne(k, :, :) = odd_row(k, :, :);
            end
        end

        averaged_neighbors = (unreceptive_ne(:,:,1) * 0.5) + (unreceptive_ne(:,:,2) * (1/12)) + ...
            (unreceptive_ne(:,:,3) * (1/12)) + (unreceptive_ne(:,:,4) * (1/12)) + ...
            (unreceptive_ne(:,:,5) * (1/12)) + (unreceptive_ne(:,:,6) * (1/12)) + ...
            (unreceptive_ne(:,:,7) * (1/12));

        final_update = averaged_neighbors + receptive_update;
        
        % Mise à jour de l'état du cristal
        A = final_update;
        
        % Arrêter les itérations si le cristal atteint le bord de la boîte
        if any(A(1,:) >= 1) || any(A(end,:) >= 1) || any(A(:,1) >= 1) || any(A(:,end) >= 1)
            fprintf('Crystal reached the boundary at iteration %d.\n', i);
            break;
        end
    end
end


function A = generate_flake_seuil(iterations, beta, gamma1, gamma2, threshold, crystal_size)
    % Initialisation de la matrice d'état
    A = beta * ones(crystal_size*2 + 1, crystal_size*2 + 1);
    A(crystal_size + 1, crystal_size + 1) = 1; % Centre de cristallisation

    x = -crystal_size:1:crystal_size;
    y = crystal_size:-1:-crystal_size;
    [X, Y] = meshgrid(x, y);

    % Décalage sur les lignes impaires pour simuler l'empilement hexagonal
    for i = 1:2:crystal_size*2 + 1
        X(i, :) = X(i, :) + 0.5;
    end

    figure; % Ouvrir une nouvelle figure pour l'affichage dynamique

    for i = 1:iterations
        % Choix du gamma selon l'itération
        if i <= threshold
            gamma = gamma1;
        else
            gamma = gamma2;
        end

        [odd_row, even_row] = neighbor_array(A, beta);

        receptive_odd = sum((odd_row >= 1), 3) > 0;
        receptive_even = sum((even_row >= 1), 3) > 0;
        receptive_logic = zeros(size(A));

        for k = 1:size(A, 1)
            if mod(k, 2) == 0
                receptive_logic(k, :) = receptive_even(k, :);
            else
                receptive_logic(k, :) = receptive_odd(k, :);
            end
        end

        unreceptive_logic = receptive_logic == 0;
        receptive = receptive_logic .* A;
        unreceptive = unreceptive_logic .* A;

        receptive_update = (receptive_logic * gamma) + receptive;
        [odd_row, even_row] = neighbor_array(unreceptive, beta);

        unreceptive_ne = zeros(size(odd_row));
        for k = 1:size(A, 1)
            if mod(k, 2) == 0
                unreceptive_ne(k, :, :) = even_row(k, :, :);
            else
                unreceptive_ne(k, :, :) = odd_row(k, :, :);
            end
        end

        averaged_neighbors = (unreceptive_ne(:,:,1) * 0.5) + (unreceptive_ne(:,:,2) * (1/12)) + ...
            (unreceptive_ne(:,:,3) * (1/12)) + (unreceptive_ne(:,:,4) * (1/12)) + ...
            (unreceptive_ne(:,:,5) * (1/12)) + (unreceptive_ne(:,:,6) * (1/12)) + ...
            (unreceptive_ne(:,:,7) * (1/12));

        final_update = averaged_neighbors + receptive_update;

        % Mise à jour de l'état du cristal
        A = final_update;

        % 🔹 **Affichage dynamique après chaque itération**
        drawing(A, beta);
        drawnow;

        % Vérifier si le cristal atteint le bord de la boîte
        if any(A(1,:) >= 1) || any(A(end,:) >= 1) || any(A(:,1) >= 1) || any(A(:,end) >= 1)
            fprintf('Crystal reached the boundary at iteration %d.\n', i);
            break;
        end
    end
end


function A = generate_flake_radial(iterations, beta, gamma1, gamma2, crystal_size)
    % Initialisation de la matrice hexagonale
    A = beta * ones(crystal_size*2 + 1, crystal_size*2 + 1);
    A(crystal_size + 1, crystal_size + 1) = 1; % Centre de cristallisation

    % Génération des coordonnées hexagonales
    x = -crystal_size:crystal_size;
    y = -crystal_size:crystal_size;
    [X, Y] = meshgrid(x, y);

    % Décalage sur les lignes impaires pour simuler une grille hexagonale
    for i = 1:2:crystal_size*2 + 1
        X(i, :) = X(i, :) + 0.5;
    end

    % Centre de la grille
    center_q = 0;
    center_r = 0;

    % Conversion des coordonnées cartésiennes en coordonnées axiales
    Q = X; % q-coordonnée
    R = Y - floor(X / 2); % r-coordonnée (axiale)

    % Calcul de la distance hexagonale Dist(q, r)
    D = (abs(Q - center_q) + abs(Q + R - (center_q + center_r)) + abs(R - center_r)) / 2;
    D_max = max(D(:)); % Distance max (coin de la grille)

    figure; % Ouvrir une nouvelle figure pour affichage dynamique

    for i = 1:iterations
        % 🔹 **Gamma dépend maintenant de la distance hexagonale**
        Gamma = gamma1 + (gamma2 - gamma1) .* (D / D_max);

        [odd_row, even_row] = neighbor_array(A, beta);

        receptive_odd = sum((odd_row >= 1), 3) > 0;
        receptive_even = sum((even_row >= 1), 3) > 0;
        receptive_logic = zeros(size(A));

        for k = 1:size(A, 1)
            if mod(k, 2) == 0
                receptive_logic(k, :) = receptive_even(k, :);
            else
                receptive_logic(k, :) = receptive_odd(k, :);
            end
        end

        unreceptive_logic = receptive_logic == 0;
        receptive = receptive_logic .* A;
        unreceptive = unreceptive_logic .* A;

        receptive_update = (receptive_logic .* Gamma) + receptive;
        [odd_row, even_row] = neighbor_array(unreceptive, beta);

        unreceptive_ne = zeros(size(odd_row));
        for k = 1:size(A, 1)
            if mod(k, 2) == 0
                unreceptive_ne(k, :, :) = even_row(k, :, :);
            else
                unreceptive_ne(k, :, :) = odd_row(k, :, :);
            end
        end

        averaged_neighbors = (unreceptive_ne(:,:,1) * 0.5) + (unreceptive_ne(:,:,2) * (1/12)) + ...
            (unreceptive_ne(:,:,3) * (1/12)) + (unreceptive_ne(:,:,4) * (1/12)) + ...
            (unreceptive_ne(:,:,5) * (1/12)) + (unreceptive_ne(:,:,6) * (1/12)) + ...
            (unreceptive_ne(:,:,7) * (1/12));

        final_update = averaged_neighbors + receptive_update;

        % Mise à jour de l'état du cristal
        A = final_update;

        % 🔹 **Affichage dynamique**
        drawing(A, beta);
        drawnow;

        % Vérifier si le cristal atteint le bord
        if any(A(1,:) >= 1) || any(A(end,:) >= 1) || any(A(:,1) >= 1) || any(A(:,end) >= 1)
            fprintf('Crystal reached the boundary at iteration %d.\n', i);
            break;
        end
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


function drawing(A, beta)

low = A < 1;
low_array = (low .* A) .* 0.5;
high = (low == 0);

temp_high = 1 - abs((high .* A) - 1);
logic_temp_high = (temp_high <= 0.5) .* high;
high_array = ((temp_high > 0.5).*temp_high) +  (0.5*logic_temp_high);

A = low_array + high_array;
A = kron(kron(A,ones(2,1)), ones(1,2));

back = beta*ones(size(A,1) , size(A,2)+1);
for ii = 1:size(A,1)
    choice = ceil(ii/2);
    if (mod(choice,2)==1) 
        back(ii, 2:end) = A(ii,:);
    else
        back(ii, 1:end-1) = A(ii,:);
    end
end
A = back;

imagesc(A)
colormap('sky')

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





