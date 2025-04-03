% Calcul de la moyenne sur i (pour chaque j)
D_mean = mean(D_values,1);

% Définition des indices j (qui représentent gamma γ)
j_values = gammas;

% Vérifier et filtrer les valeurs positives uniquement
valid_indices = (j_values > 0) & (D_mean > 0);
j_values_valid = j_values(valid_indices);
D_mean_valid = D_mean(valid_indices);

% Vérification si des données valides restent
if isempty(j_values_valid)
    error('Toutes les valeurs de gammas ou D_mean sont négatives ou nulles, impossible de faire un fit en loi de puissance.');
end

% Transformation logarithmique pour ajustement de loi de puissance
log_j = log(j_values_valid);
log_D_mean = log(D_mean_valid);

% Régression linéaire sur les données transformées : log(D) = B log(j) + log(A)
coeffs = polyfit(log_j, log_D_mean, 1); % Ajustement linéaire sur les logs

% Extraction des paramètres de la loi de puissance
B = coeffs(1);          % Exposant de la loi de puissance
A = exp(coeffs(2));     % Facteur de proportionnalité

% Calcul de l'ajustement sur les données valides
D_fit = A * (j_values_valid.^B);

% Affichage des résultats
% Affichage des résultats
plot(j_values, D_mean, 'bo', 'MarkerFaceColor', 'b'); % Points de données d'origine
hold on;
plot(j_values_valid, D_fit, 'r:', 'LineWidth', 4); % Fit en loi de puissance (dotted line)

% Activation du mode LaTeX pour les axes
xlabel('$\gamma$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Moyenne de $D$ sur $\beta$', 'Interpreter', 'latex', 'FontSize', 14);
title('Ajustement par une loi de puissance', 'Interpreter', 'latex', 'FontSize', 14);
legend('Moyenne des données', 'Fit: $A \gamma^B$', 'Location', 'best', 'Interpreter', 'latex');

grid on;
hold off;

% Affichage des coefficients de la loi de puissance
fprintf('Facteur de proportionnalité (A) : %f\n', A);
fprintf('Exposant (B) : %f\n', B);
