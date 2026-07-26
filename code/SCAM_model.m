clear all
close all

% STOCHASTIC CELLULAR AUTOMATA MODEL (SCAM)
% Simulation of avascular tumour growth under the effects of chemotherapy


% ----- Parameters ------------------------------------------------------------------------------------------------------


% Generic parameters for the model
L = 200;  % size of the domain (cells per side)
center = floor(L/2);
N_total_cells = L^2;  % number of total cells
N_PCcells = 1;  % number of proliferative tumor cells
N_ICcells = 0.001 * N_total_cells;  % number of immune cells
N_DRCcells = 0.1 * N_total_cells;  % number of drug-resistant cells
num_steps = 150;

% Growth variables
base_p_zero = 0.7;  % base probability of division of mutant PT cell
a = 0.42;  % base necrotic thickness
b = 0.53;  % base proliferative (living tumor) thickness
R_max = 37.5;  % maximum tumor extent
p_dt = 0.5;  % tumor death constant
p_di = 0.2;  % immune death constant

% Treatment variables
K_c = 1;  % chemotherapy effect on the division
Y_pc = 0.90;  % PC's resistance to treatment 
Y_qc = 0.7;  % QC's resistance to treatment 
Y_ic = 0.60;  % IC's resistance to treatment
Y_noc = 0.5;  % NoC's resistance to treatment
k_pc = 0.8;  % PC's death rate due to therapy
k_qc = 0.4;  % QC's death rate due to therapy
k_ic = 0.6;  % IC's death rate due to therapy
k_noc = 0.9;  % NoC's death rate due to therapy
c_i = 0.75;  % attenuation coefficient of a drug for any cell type
PK = 1;  % pharmacokinetics
t_ap = 60;  % start time of therapy infusion
t_per = 10;  % time interval between injections
n_d = 3;  % number of the cycles (infusions) of the medicine
tau = 20;  % time constant of each dose
g = 1;  % drug concentration
theta = rand * (1 - eps) + eps;
gamma_pc = theta * Y_pc;  % chemotherapy-induced resistance for PC cells
gamma_qc = theta * Y_qc;  % chemotherapy-induced resistance for QC cells
gamma_ic = theta * Y_ic;  % chemotherapy-induced resistance for IC cells
gamma_noc = theta * Y_noc;  % chemotherapy-induced resistance for IC cells



% ----- Simulation ------------------------------------------------------------------------------------------------------


% Initialize the domain and the relative grids 
[cells,age_cells,drug_resistant_cells,survival_steps, n_dead_steps]  = initialize_grid(L, center, N_ICcells, N_DRCcells, tau, n_d);

IC_successes = 0;
IC_failures = 0;


% % Plot for the initial configuration of the tissue
% figure;
% imagesc(cells);
% colormap([0 0 0;  % black for NoC/ES cells (state 0)
%           0.5 0.5 0.5;  % medium grey for PC cells (state 1)
%           0.75 0.75 0.75;  % light grey for QC cells (state 2)
%           0.25 0.25 0.25;  % dark grey for NeC cells (state 3)
%           0.1 0.1 0.1;  % very dark grey for US cells (state 4)
%           1 1 0.8;  % cream for IC cells (state 5)
%           1 1 1]);  % white for DC cells (state 6)
% title('Initial Tumor State');


% Create a list of all cell indices
indices = [];
for i = 1:L
    for j = 1:L
        indices = [indices; i, j];
    end
end

Number_PC_values = [];
Number_QC_values = [];
Number_NeC_values = [];


% Simulation loop
for step = 1:num_steps

    Number_T = sum(cells(:) == 1) + sum(cells(:) == 2) + sum(cells(:) == 3);  % total number of tumor cells (PC+QC+NeC) in the tissue
    Number_PC = sum(cells(:) == 1);  % total number of PC cells in the tissue
    Number_QC = sum(cells(:) == 2);  % total number of QC cells in the tissue
    Number_NeC = sum(cells(:) == 3);  % total number of NeC cells in the tissue
    Number_IC = sum(cells(:) == 5);  % total number of IC cells in the tissue

    Number_PC_values = [Number_PC_values; Number_PC];
    Number_QC_values = [Number_QC_values; Number_QC];
    Number_NeC_values = [Number_NeC_values; Number_NeC];

    if step < t_ap
        R_t = compute_average_radius(cells, L, center);
    else
        R_t = R_t;
    end

    delta_n = necrotic_layer(R_t, a);
    delta_p = proliferating_layer(R_t, b);
    R_n = necrotic_radius(R_t, delta_n, delta_p);

    number_newborn_IC_cells = get_growth_rate(IC_successes, IC_failures, Number_PC, Number_T);

    % initialize new immune cells at the corners of the lattice
    if number_newborn_IC_cells > 0
        for cell = 1:length(number_newborn_IC_cells)
            corner = randi([1, 4]);
            if corner == 1
                while true
                    rand_x = randi([1, 10]);
                    rand_y = randi([L-10, L]);
                    if cells(rand_x, rand_y) == 0
                        cells(rand_x, rand_y) = 5;
                        age_cells(rand_x, rand_y) = 0;
                        break;
                    end
                end
            elseif corner == 2
                while true
                    rand_x = randi([1, 10]);
                    rand_y = randi([1, 10]);
                    if cells(rand_x, rand_y) == 0
                        cells(rand_x, rand_y) = 5;
                        age_cells(rand_x, rand_y) = 0;
                        break;
                    end
                end
            elseif corner == 3
                while true
                    rand_x = randi([L-10, L]);
                    rand_y = randi([1, 10]);
                    if cells(rand_x, rand_y) == 0
                        cells(rand_x, rand_y) = 5;
                        age_cells(rand_x, rand_y) = 0;
                        break;
                    end
                end
            elseif corner == 4
                while true
                    rand_x = randi([L-10, L]);
                    rand_y = randi([L-10, L]);
                    if cells(rand_x, rand_y) == 0
                        cells(rand_x, rand_y) = 5;
                        age_cells(rand_x, rand_y) = 0;
                        break;
                    end
                end
            end
        end
    end

    % apply chemotherapy only if the step is within the interval for therapy cycle
    applied = false;
    for cycle = 0:(n_d-1) 
        start_time = t_ap + cycle * (tau + t_per);
        end_time = start_time + tau;
        if step >= start_time && step <= end_time
            [cells, age_cells, drug_resistant_cells, survival_steps, IC_successes, IC_failures] = update_grid(true, indices, step, cells, age_cells, ...
                drug_resistant_cells, survival_steps, n_dead_steps, L, center, R_t, delta_p, R_n, p_dt, p_di, base_p_zero, R_max, ...
                K_c, k_pc, k_qc, k_ic, k_noc, c_i, PK, n_d, tau, g, gamma_pc, gamma_qc, gamma_ic, gamma_noc, cycle);
            applied = true;
            break;
        end
    end
    if ~applied
        [cells, age_cells, drug_resistant_cells, survival_steps, IC_successes, IC_failures] = update_grid(false, indices, step, cells, age_cells, ...
            drug_resistant_cells, survival_steps, n_dead_steps, L, center, R_t, delta_p, R_n, p_dt, p_di, base_p_zero, R_max, ...
            K_c, k_pc, k_qc, k_ic, k_noc, c_i, PK, n_d, tau, g, gamma_pc, gamma_qc, gamma_ic, gamma_noc, cycle);
    end

    imagesc(cells);
    colormap([0 0 0;  % black for NoC/ES cells (state 0)
          0.5 0.5 0.5;  % medium grey for PC cells (state 1)
          0.75 0.75 0.75;  % light grey for QC cells (state 2)
          0.25 0.25 0.25;  % dark grey for NeC cells (state 3)
          0.1 0.1 0.1;  % very dark grey for US cells (state 4)
          1 1 0.8;  % cream for IC cells (state 5)
          1 1 1]);  % white for DC cells (state 6)
    title(['Tumor State at step ', num2str(step)]);
    pause(0.1);

end

% % Plot for the evolution of the total numbers of tumor cells
% figure;
% plot(1:num_steps, Number_PC_values, '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 2);  % Plot PC values in custom color
% hold on;
% plot(1:num_steps, Number_QC_values, '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 2);  % Plot QC values in custom color
% plot(1:num_steps, Number_NeC_values, '-', 'Color', [0 0 0], 'LineWidth', 2);  % Plot NeC values in custom color
% hold off;
% xlabel('Iteration');
% ylabel('Number of Cells');
% title('Tumor cells counts over iterations');
% legend('PC', 'QC', 'NeC');


% ----- Definition of the main function -------------------------------------------------------------------------------------------


% Function to update the system
function [new_cells, new_age_cells, new_drug_resistant_cells, new_survival_steps, IC_successes, IC_failures] = update_grid(apply_therapy, indices, step, cells, ...
    age_cells, drug_resistant_cells, survival_steps, n_dead_steps, L, center, R_t, delta_p, R_n, p_dt, p_di, base_p_zero, R_max, ...
    K_c, k_pc, k_qc, k_ic, k_noc, c_i, PK, n_d, tau, g, gamma_pc, gamma_qc, gamma_ic, gamma_noc, num_cycle)

    % Grids initialization
    new_cells = cells;
    new_age_cells = age_cells;
    new_drug_resistant_cells = drug_resistant_cells;
    new_survival_steps = survival_steps;

    % Number of successes and failures for IC cells in killing proliferative cells
    IC_successes = 0;
    IC_failures = 0;

    % Consider a random sequence for cell update
    shuffled_indices = indices(randperm(size(indices, 1)), :);
    for k = 1:size(shuffled_indices, 1)
        i = shuffled_indices(k, 1);
        j = shuffled_indices(k, 2);

        % Case for proliferative cell (PC cell)
        if cells(i,j) == 1
            r = sqrt((i - center).^2 + (j - center).^2);
            br = division_probability(apply_therapy, base_p_zero,r, R_max, K_c, gamma_pc, n_d, n_dead_steps(i,j));  % probability of division
            if rand() < br
                % proliferation
                neighbor_cells = get_neighbor_cells(i, j, L);
                normal_cells = get_normal_cells(new_cells, neighbor_cells, 0);
                if ~isempty(normal_cells)
                    daughter_cell = normal_cells{randi(length(normal_cells))};
                    new_cells(i,j) = 1;  % first daughter cell remains in the parent's position
                    new_age_cells(i,j) = 0; 
                    new_cells(daughter_cell(1), daughter_cell(2)) = 1;  % second daughter cell goes to a random '0' cell in the neighborhood
                    new_age_cells(daughter_cell(1), daughter_cell(2)) = 0; 
                    if drug_resistant_cells(i,j) == 1
                        % symmetric or asymmetric replication for drug resistant cells
                        if rand() < 0.5
                            new_drug_resistant_cells(i,j) = 1;
                            new_drug_resistant_cells(daughter_cell(1), daughter_cell(2)) = 1;
                        else
                            new_drug_resistant_cells(i,j) = 1;
                            new_drug_resistant_cells(daughter_cell(1), daughter_cell(2)) = 0;
                        end
                    else
                        % symmetric replication for for drug sensitive cells
                        new_drug_resistant_cells(i,j) = 0;
                        new_drug_resistant_cells(daughter_cell(1), daughter_cell(2)) = 0;
                    end
                else
                    new_cells(i,j) = 1;
                    new_age_cells(i,j) = age_cells(i,j) + 1;
                end
            else
                new_cells(i,j) = 1;
                new_age_cells(i,j) = age_cells(i,j) + 1;
            end
                
            if r < (R_t - delta_p)
                new_cells(i,j) = 2;
            end

            if apply_therapy == true
                if drug_resistant_cells(i,j) == 0
                    drug_killing_rate_PC = get_therapy_response_curve(step, k_pc, g, gamma_pc, num_cycle, PK, c_i, tau);
                    if rand() < drug_killing_rate_PC
                        % cell dies due to the chemotherapy
                        new_cells(i,j) = 6; 
                        new_age_cells(i,j) = 0;
                        new_drug_resistant_cells(i,j) = 0;
                        new_survival_steps(i,j) = 0;
                    else
                        new_survival_steps(i,j) = survival_steps(i,j) + 1;
                    end
                else
                    new_survival_steps(i,j) = survival_steps(i,j) + 1;
                end
                    
                if new_survival_steps(i,j) > n_dead_steps(i,j)
                    new_cells(i,j) = 4; 
                    new_age_cells(i,j) = age_cells(i,j);
                    new_drug_resistant_cells(i,j) = 0;
                    new_survival_steps(i,j) = 0;
                end
            end

        % Case for quiescent cell (QC cell)
        elseif cells(i,j) == 2
            r = sqrt((i - center).^2 + (j - center).^2);
            if r < R_n
                % turn into a NC due to the lack of nutrition
                new_cells(i,j) = 3;
                new_age_cells(i,j) = age_cells(i,j) + 1;
            elseif r > (R_t - delta_p) && r < R_t
                % turn into a PC due to accessing sufficient nutrient
                new_cells(i,j) = 1;
                new_age_cells(i,j) = age_cells(i,j) + 1;
            else
                new_cells(i,j) = 2;
                new_age_cells(i,j) = age_cells(i,j) + 1;
            end

            if apply_therapy == true
                if drug_resistant_cells(i,j) == 0
                    drug_killing_rate_QC = get_therapy_response_curve(step, k_qc, g, gamma_qc, num_cycle, PK, c_i, tau);
                    if rand() < drug_killing_rate_QC
                        % cell dies due to the chemotherapy
                        new_cells(i,j) = 6; 
                        new_age_cells(i,j) = 0;
                        new_drug_resistant_cells(i,j) = 0;
                        new_survival_steps(i,j) = 0;
                    else
                        new_survival_steps(i,j) = survival_steps(i,j) + 1;
                    end
                else
                    new_survival_steps(i,j) = survival_steps(i,j) + 1;
                end
                    
                if new_survival_steps(i,j) > n_dead_steps(i,j)
                    new_cells(i,j) = 4; 
                    new_age_cells(i,j) = age_cells(i,j);
                    new_drug_resistant_cells(i,j) = 0;
                    new_survival_steps(i,j) = 0;
                end
            end
    
        % Case for necrotic cell (NeC cell)
        elseif cells(i,j) == 3
            new_cells(i,j) = 3;
            new_age_cells(i,j) = age_cells(i,j) + 1;

        % Case for unstable state (US cell)
        elseif cells(i,j) == 4
            new_cells(i,j) = 4;
            new_age_cells(i,j) = age_cells(i,j) + 1;
            if new_age_cells(i,j) > n_dead_steps(i,j)
                new_cells(i,j) = 6;
            end

        % Case for immune cell (IC cell)
        elseif cells(i,j) == 5

            new_position = [i,j];

            if rand() < 0.2
                if i < center
                    new_position(1) = i + 1;
                elseif i > center
                    new_position(1) = i - 1;
                end
                if j < center
                    new_position(2) = j + 1;
                elseif j > center
                    new_position(2) = j - 1;
                end
    
                if new_cells(new_position(1), new_position(2)) == 0
                    new_cells(new_position(1), new_position(2)) = 5;
                    new_age_cells(new_position(1), new_position(2)) = age_cells(i,j) + 1;
                    new_cells(i, j) = 0;
                    new_age_cells(i, j) = 0;
                else
                    neighbor_cells = get_neighbor_cells(i, j, L);
                    normal_cells = get_normal_cells(new_cells, neighbor_cells, 0);
                    if ~isempty(normal_cells)
                        new_position = normal_cells{randi(length(normal_cells))};
                        new_cells(i,j) = 0;
                        new_age_cells(i,j) = 0;
                        new_cells(new_position(1), new_position(2)) = 5;
                        new_age_cells(new_position(1), new_position(2)) = age_cells(i,j) + 1;
                    end
                end
            end
           
            neighbor_cells = get_neighbor_cells(new_position(1), new_position(2), L);
            neighbor_ESs = get_normal_cells(new_cells, neighbor_cells, 0);
            neighbor_PCs = get_normal_cells(new_cells, neighbor_cells, 1);
            neighbor_immune_cells = get_normal_cells(new_cells, neighbor_cells, 5);

            r_i = antitumor_probability(p_dt, length(neighbor_immune_cells), length(neighbor_PCs));
            r_t = protumor_probability(p_di, length(neighbor_immune_cells), length(neighbor_PCs));

            if ~isempty(neighbor_PCs)
                success_in_killing = false;
                for time = 1:length(neighbor_PCs)
                    if rand() < r_i
                        proliferative_cell = neighbor_PCs{time};
                        new_cells(i,j) = 0;  
                        new_age_cells(i,j) = 0; 
                        new_cells(proliferative_cell(1), proliferative_cell(2)) = 5;  
                        new_age_cells(proliferative_cell(1), proliferative_cell(2)) = age_cells(i,j) + 1;
                        success_in_killing = true;
                        IC_successes = IC_successes + 1;
                        break;
                    end
                end
                if ~success_in_killing
                    IC_failures = IC_failures + 1;
                    if rand() < r_t
                        new_cells(i,j) = 0;  
                        new_age_cells(i,j) = 0; 
                    else
                        new_cells(i,j) = 5;  
                        new_age_cells(i,j) = age_cells(i,j) + 1; 
                    end
                end
            end

            if apply_therapy == true
                if drug_resistant_cells(i,j) == 0
                    drug_killing_rate_IC = get_therapy_response_curve(step, k_ic, g, gamma_ic, num_cycle, PK, c_i, tau);
                    if rand() < drug_killing_rate_IC
                        % cell dies due to the chemotherapy
                        new_cells(i,j) = 6; 
                        new_age_cells(i,j) = 0;
                        new_survival_steps(i,j) = 0;
                    else
                        new_survival_steps(i,j) = survival_steps(i,j) + 1;
                    end
                else
                    new_survival_steps(i,j) = survival_steps(i,j) + 1;
                end
                if new_survival_steps(i,j) > n_dead_steps(i,j)
                    new_cells(i,j) = 4; 
                    new_age_cells(i,j) = age_cells(i,j);
                    new_survival_steps(i,j) = 0;
                end
            end
        
        % Case for dead cell (DC cell)
        elseif cells(i,j) == 6
            neighbor_cells = get_neighbor_cells(i, j, L);
            normal_cells = get_normal_cells(new_cells, neighbor_cells, 0);
            if length(normal_cells) > 5
                new_cells(i,j) = 0;
                new_age_cells(i, j) = 0;
            else
                new_cells(i,j) = 6;
                new_age_cells(i, j) = age_cells(i, j) + 1;
            end

        % Case for normal cell (NoC cell)
        elseif cells(i,j) == 0
            if apply_therapy == true
                if drug_resistant_cells(i,j) == 0
                    drug_killing_rate = get_therapy_response_curve(step, k_noc, g, gamma_noc, num_cycle, PK, c_i, tau);
                    if rand() < drug_killing_rate
                        % cell dies due to the chemotherapy
                        new_cells(i,j) = 6; 
                        new_age_cells(i,j) = 0;
                        new_survival_steps(i,j) = 0;
                    else
                        new_survival_steps(i,j) = survival_steps(i,j) + 1;
                    end
                else
                    new_survival_steps(i,j) = survival_steps(i,j) + 1;
                end
                if new_survival_steps(i,j) > n_dead_steps(i,j)
                    new_cells(i,j) = 4; 
                    new_age_cells(i,j) = age_cells(i,j);
                    new_survival_steps(i,j) = 0;
                end
            end
        end
    end
end


% ----- Definition of other functions ---------------------------------------------------------------------------------------------


% Function for the initialization of the grid 
function [cells,age_cells,drug_resistant_cells,survival_steps,n_dead_steps]  = initialize_grid(L, center, N_ICcells, N_DRCcells, tau, n_d)
    cells = zeros(L);
    age_cells = zeros(size(cells));
    drug_resistant_cells = zeros(size(cells));
    survival_steps = zeros(size(cells));
    n_dead_steps = zeros(size(cells));

    % initialize proliferative tumor cells
    cells(center, center) = 1;

    % initialize random immune cells
    for i = 1:N_ICcells
        rand_x = randi(L);
        rand_y = randi(L);
        if cells(rand_x, rand_y) == 0 
            cells(rand_x, rand_y) = 5;
        end
    end

    % initialize drug-resistant cells
    for i = 1:N_DRCcells
        rand_x = randi(L);
        rand_y = randi(L);
        if cells(rand_x, rand_y) ~= 3
            drug_resistant_cells(rand_x, rand_y) = 1;
        end
    end

    % initialize n_dead values (number of time steps before death due to the treatment)
    for i = 1:L
        for j = 1:L
            n_dead_steps(i, j) = randi([1, (tau*n_d)]);
        end
    end
end

% Function to define the Moore neighborhood of a cell
function neighbor_cells = get_neighbor_cells(i, j, L)
    neighbor_cells = {};
    for dx = [-1, 0, 1]
        for dy = [-1, 0, 1]
            ni = mod(i + dx - 1, L) + 1;
            nj = mod(j + dy - 1, L) + 1;
            if ~isequal([ni, nj], [i, j])
                neighbor_cells{end+1} = [ni, nj];
            end
        end   
    end
end

% Function to define the Von Neumann neighborhood of a cell
function neighbors = get_von_neumann_neighbors(i, j, L)
    neighbors = {};
    if i > 1
        neighbors{end+1} = [i-1, j];
    else
        neighbors{end+1} = [L, j];
    end
    if i < L
        neighbors{end+1} = [i+1, j];
    else
        neighbors{end+1} = [1, j];
    end
    if j > 1
        neighbors{end+1} = [i, j-1];
    else
        neighbors{end+1} = [i, L];
    end
    if j < L
        neighbors{end+1} = [i, j+1];
    else
        neighbors{end+1} = [i, 1];
    end
end

% Function to store neighbors equal to a certain value in an array
function normal_cells = get_normal_cells(grid, neighborhood, value)
    normal_cells = {};
    for k = 1:numel(neighborhood)
        neighbor = neighborhood{k};
        ni = neighbor(1);
        nj = neighbor(2);
        if grid(ni, nj) == value
            normal_cells{end+1} = [ni, nj];
        end
    end
end

% Probability of division for proliferation
function br = division_probability(apply_therapy, base_p_zero,r, R_max, K_c, gamma_pc, n_d, n_dead)
    if apply_therapy == false
        br = base_p_zero * (1 - (r / (R_max - K_c)));
    else
        p_zero = (base_p_zero * gamma_pc) / (n_d ^ (1/n_dead));
        br = p_zero * (1 - (r / (R_max - K_c)));
    end
end

% Average radius of the tumor
function R_t = compute_average_radius(cells, L, center)
    edge_cells = [];
    for i = 1:L
        for j = 1:L
            if cells(i, j) == 1 
                neighbors = get_von_neumann_neighbors(i, j, L);
                normal_cells = get_normal_cells(cells, neighbors, 0);
                if ~isempty(normal_cells)
                    edge_cells = [edge_cells; i, j];
                    break;
                end
            end
        end
    end
    distances = sqrt((edge_cells(:, 1) - center).^2 + (edge_cells(:, 2) - center).^2);
    R_t = mean(distances);
end

% Thickness of necrotic cells (determines necrotic fraction)
function delta_n = necrotic_layer(R_t, a)
    delta_n = a * (R_t ^ (2/3));
end

% Thickness of proliferating cancerous cells (determines proliferative fraction)
function delta_p = proliferating_layer(R_t, b)
    delta_p = b * (R_t ^ (2/3));
end

% Average necrotic layer radius
function R_n = necrotic_radius(R_t, delta_n, delta_p)
    R_n = R_t - (delta_n + delta_p);
end

% Probability of turning from IC to unstable state
function r_i = antitumor_probability(p_dt, neighbor_immune_cells, neighbor_PC_cells)
    r_i = p_dt * (neighbor_immune_cells / neighbor_PC_cells);
end

% Probability of dying for IC cells
function r_t = protumor_probability(p_di, neighbor_immune_cells, neighbor_PC_cells)
    r_t = p_di * (neighbor_PC_cells / neighbor_immune_cells);
end

% Positive nonlinear growth rate of recruiting ICs
function IC_rate_new_cells = get_growth_rate(n_successes, n_failures, n_PC_cells, n_tumor_cells)
    IC_rate_new_cells = (n_successes - n_failures) * (n_PC_cells / n_tumor_cells);
end

% Exponential treatment response curve
function Fg_i = get_therapy_response_curve(t, kill_rate, g, resistance_factor, num_cycle, PK, attenuation_coefficient, tau)
    l_i = (kill_rate * g) / (resistance_factor * (num_cycle+1) + 1);
    Fg_i = l_i * PK * exp(- attenuation_coefficient * (t - ((num_cycle+1) * tau)));
end









    