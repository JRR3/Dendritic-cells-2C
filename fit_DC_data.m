%%%Source        : Houston Methodist Research Institute
%%%Location      : Houston, TX.
%%%Origin        : September 28, 2018
%%%PI            : Mauro Ferrari
%%%Supervisor    : Vittorio Cristini
%%%Collaborator  : Prashant Dogra
%%%Developer     : Javier Ruiz Ramirez

%%%================================================================

function fit_DC_data
%%%Driver

clc;
close all;
% labels = {'LN','Tumor'};
s = read_xls_data;

%Extract data
[time_vector, data, standard_deviation] = extract_experimental_data(s);

%Data is in microliters
volume_of_individual_popliteal_lymph_node = 0.78;
volume_of_individual_tumor                = ...
    3.5 * volume_of_individual_popliteal_lymph_node;

v_lymph_node = volume_of_individual_popliteal_lymph_node * 2;
v_tumor      = volume_of_individual_tumor                * 2;
v_total      = v_lymph_node + v_tumor;

compartment_volume    = [v_lymph_node, v_tumor];
initial_conditions    = [0, 0];

%Initial guess for the optimization routine
t_off        = 4;
g_rate       = 1;
M            = 1;

guess(1) = t_off;
guess(2) = g_rate;
guess(3) = M;

list_of_parameters{1} = compartment_volume;
list_of_parameters{2} = initial_conditions;
list_of_parameters{3} = time_vector;
list_of_parameters{4} = data;
list_of_parameters{5} = standard_deviation;

solve_with_statistics(list_of_parameters, guess);

%%%================================================================

function solve_with_statistics(list_of_parameters, guess)

objective_function = @(p) weighted_function(p, list_of_parameters);

labels_for_parameters = {'t_off', 'g_rate', 'M'};
labels = {'LN','Tumor'};

data = list_of_parameters{4};

n_observations = numel(data);

n_parameters = length(guess);
[optimal_p, resnorm, residual, exitflag, output, lambda, jacobian] =...
    lsqnonlin(...
    objective_function, ...
    guess);

t_off        = optimal_p(1);
g_rate       = optimal_p(2);
M            = optimal_p(3);

jacobian  = full(jacobian);
i_hessian = jacobian.' * jacobian;
hessian   = i_hessian \ eye(n_parameters);
mse       = resnorm / (n_observations - n_parameters);
covariance_matrix = hessian * mse;

format long e
se_manual = sqrt(diag(covariance_matrix));

%Iterate over each of the optimization parameters and compute the standard 
%error
for k = 1:length(se_manual)
    disp([labels_for_parameters{k}, '=']);
    optimal_p(k)
    disp(['Standard error: ']);
    se_manual(k)  
    disp(['-------------------------'])
end

format short g;

disp(['Residual norm: ', num2str(resnorm)]);

v = list_of_parameters{1};
v_lymph_node = v(1);
v_tumor      = v(2);
v_total      = v_lymph_node + v_tumor;

w = g_rate * (1/v_tumor + 1/v_lymph_node)
M_over_v = M / v_total
% a_constant = M_over_v * (1 + v_tumor / v_lymph_node)
% b_constant = M * v_tumor / (v_total * v_lymph_node)
plot_numerical_data(optimal_p, list_of_parameters, labels);

%%%TODO: Compare numerical results with physical lengths and flow speeds
mouse_pixel_length = 871;
mouse_length       = 10;
pixels_to_cm       = mouse_length / mouse_pixel_length;
injection_site_to_tumor = 266 * pixels_to_cm;

%%%================================================================

function r = weighted_function(decision_variables, list_of_parameters)

compartment_volume = list_of_parameters{1};
initial_conditions = list_of_parameters{2};
time_vector        = list_of_parameters{3};
expected_data      = list_of_parameters{4};

time_vector        = [0, time_vector];
n_tissues          = size(expected_data, 2);
% Note that an additional column has to be created due to the fact
% that a zero was added to the time vector
expected_data      = [zeros(1,n_tissues); expected_data];

f = @(t,y) ode_sys(t, y,...
    compartment_volume,...
    decision_variables);

[~,y] = ode15s(f, time_vector, initial_conditions);

r = (y - expected_data);

%%%================================================================

function dy = ode_sys(t, y, compartment_volume, variable)

t_off        = variable(1);
g_rate       = variable(2);
M            = variable(3);

v_lymph_node = compartment_volume(1);
v_tumor      = compartment_volume(2);

%Smoothing the step function using a sigmoidal analog
% F     = @(t) M * ( 1 - ( (t - t_off) > 0 ) );
F     = @(t) M * sigmoidal(t,t_off);
dy    = y*0;

dy(1) = F(t) - g_rate*(y(1) - y(2));
dy(2) = g_rate*(y(1) - y(2));

dy = dy ./ [v_lymph_node; v_tumor];

%%%================================================================

function test_sigmoidal()
t_off    = 10;
t        = linspace(0,24,100);
plot(t,sigmoidal(t,t_off),'linewidth',2);

%%%================================================================

function y = sigmoidal(t, t_off)
k           = 6;
argument    = -k * (t - t_off);
denominator = 1 + exp(argument);
y           = 1 - 1./denominator;

%%%================================================================

function [time_vector, data] = ...
    old_extract_experimental_data(s)

% Choose which type of manipulation is going to be applied to the data
% data_manipulation_mode = 'sum';
data_manipulation_mode = 'average';

group    = 'control';
clusters = {{'lnol','lnor'},{'tl','tr'}};
n_clusters = length(clusters);
time_vector = s.times;
% tissues = fieldnames(s.(group));
n_times   = length(time_vector);
data      = zeros(n_times, n_clusters);
counter   = 0;
% raw_data  = cell(n_times,n_clusters);
for c = 1:n_clusters
%     raw_data = [];
    cluster_size = length(clusters{c});
    for tissue_index = 1:cluster_size
        tissue = clusters{c}{tissue_index};
        data(:,c) = data(:,c) + s.(group).(tissue).mean.';
    end
    if strcmp(data_manipulation_mode, 'average')
        data(:,c) = data(:,c) / cluster_size / 100;
    end
    counter = counter + 1;
end

%%%================================================================

function [time_vector, data, standard_deviation] = ...
    extract_experimental_data(s)

%Currently only one group is avai
group          = 'control';
mouse_selector = [1,2,3];

clusters = {{'lnol','lnor'},{'tl','tr'}};
% clusters = {{'lnol'},{'tl'}};
% clusters = {{'lnor'},{'tr'}};

n_clusters = length(clusters);
time_vector = s.times;
% tissues = fieldnames(s.(group));
% n_tissues = length(tissues);
n_times   = length(time_vector);
data               = zeros(n_times, n_clusters);
standard_deviation = zeros(n_times, n_clusters);
raw  = cell(n_times,n_clusters);
factor    = 1/100;

for c = 1:n_clusters
    cluster_size = length(clusters{c});
    for tissue_index = 1:cluster_size
        tissue = clusters{c}{tissue_index};
        for time_index = 1:n_times
            raw{time_index, c} = [raw{time_index,c},...
                s.(group).(tissue).raw_data{time_index}(mouse_selector)];
        end
    end
    for time_index = 1:n_times
        data(time_index,c) = ...
            mean(raw{time_index,c});
        standard_deviation(time_index,c) = ...
            std(raw{time_index,c});
    end
    
end

data = data * factor;
standard_deviation = standard_deviation * factor;


%%%================================================================

function plot_experimental_data(time_vector, data, labels)
hold on;
color_vector = {'k','r','r','g','b','b','k'};
color_vector_size = length(color_vector);
marker_vector = {'p','<','>','*','<','>','^','x'};
data_size = size(data);
N = data_size(2);

for k = 1:N
    plot(time_vector, data(:,k),...
        'Marker', marker_vector{k},...
        'Color', color_vector{k},...
        'LineStyle', 'None',...
        'MarkerSize',6,...
        'LineWidth', 2);
end
legend(labels,'location','northwest');

set(gca, 'FontSize', 16)
xlabel('Time (hours)','FontSize',18);
ylabel('%ID/g','FontSize',18);


%%%================================================================

function plot_numerical_data(optimal_p, list_of_parameters, labels)

compartment_volume = list_of_parameters{1};
initial_conditions = list_of_parameters{2};
time_vector        = list_of_parameters{3};
expected_data      = list_of_parameters{4};
standard_deviation = list_of_parameters{5};

%Units are hours
t_interval         = [0, max(time_vector)];

f = @(t,y) ode_sys(t, y,...
    compartment_volume,...
    optimal_p);

[t,y] = ode15s(f, t_interval, initial_conditions);

hold on;
color_vector = {'k','r','r','g','b','b','k'};
marker_vector = {'p','<','>','*','<','>','^','x'};

N = size(y,2);

hold on;

legend_list = cell(1,2*N);
counter     = 0;
factor      = 100;
% Numerical approximation
for k = 1:N
    plot(t, y(:,k) * factor,...
        'LineWidth', 2,...
        'Color', color_vector{k});
    
    label = [labels{k}, '-Num'];
    counter = counter + 1;
    legend_list{counter} = label;
end


% Experimental data
for k = 1:N
    
    plot(time_vector, expected_data(:,k) * factor,...
        'Marker', marker_vector{k},...
        'Color', color_vector{k},...
        'LineStyle', 'None',...
        'MarkerSize',6,...
        'LineWidth', 2);
     
    label = [labels{k}, '-Exp'];
    counter = counter + 1;
    legend_list{counter} = label;
end


% Error bars
for k = 1:N
    
    errorbar(time_vector, ...
        expected_data(:,k) * factor, ...
        standard_deviation(:,k) * factor, ...
        'vertical',...
        'Marker', marker_vector{k},...
        'Color', color_vector{k},...
        'LineStyle', 'None',...
        'MarkerSize',6,...
        'LineWidth', 2);
    
end


legend(legend_list, 'location', 'northwest')
set(gca, 'FontSize', 16);


%%%Choose label based on the units
xlabel('Time (hours)','FontSize',18);
% ylabel('%ID/g','FontSize',18);
ylabel('Cells / 100 uL','FontSize',18);
% ylabel('Cells / uL','FontSize',18);

return;

t_0  = optimal_p(1);
g    = optimal_p(2);
M    = optimal_p(3);
V    = sum(compartment_volume)
V_L  = compartment_volume(1);
V_T  = compartment_volume(2);
w    = g * (1/V_T + 1/V_L);

% % % Derivative plot
S = @(t) exp(-w * t);
x_a = @(t) M/V * (t + (1-S(t)) * (V_T/g - 1/w));
x_a_p = @(t) M/V * (1 + S(t) * (V_T/V_L) );
x_b = @(t) M/V * ( (S(t-t_0) - S(t)) * (V_T/g - 1/w) + t_0 );
x_b_p = @(t) - M / V * (V_T / V_L) * (S(t-t_0) - S(t));
x_c = @(t) M / V * (t - (1-S(t))/w );
x_c_p = @(t) M / V * (1 - S(t));
x_d = @(t) M/V * ( t_0 - (S(t-t_0) - S(t)) / w);
x_d_p = @(t) M/V * (S(t-t_0) - S(t));

f_L       = @(t) x_a(t) .* (t < t_0) + x_b(t) .* (t_0 < t);
f_T       = @(t) x_c(t) .* (t < t_0) + x_d(t) .* (t_0 < t);
f_L_prime = @(t) x_a_p(t) .* (t < t_0) + x_b_p(t) .* (t_0 < t);
f_T_prime = @(t) x_c_p(t) .* (t < t_0) + x_d_p(t) .* (t_0 < t);

r_fun = @(t) (f_T_prime(t) .* f_L(t) - f_L_prime(t) .* f_T(t)) ./ (f_L(t).^2); 

f_cell = {f_L_prime, f_T_prime};

figure();
hold on;
counter = 0;
legend_list = cell(1,N);

for k = 1:N
%     dy_dx = gradient(y(:,k),t);
%     list_of_derivatives{k} = dy_dx;
    plot(t, f_cell{k}(t),...
        'LineWidth', 2,...
        'Color', color_vector{k});
    
    label = [labels{k}, '-rate'];
    counter = counter + 1;
    legend_list{counter} = label;
end

legend(legend_list, 'location', 'northeast')

set(gca, 'FontSize', 16)
xlabel('Time (hours)','FontSize',18);
% ylabel('%ID/g','FontSize',18);
% ylabel('Cells / 100 uL','FontSize',18);
ylabel('Cells / (h uL)','FontSize',18);

% % % Plot ratio of derivatives 
figure();
ratio = (x_c_p(t) ./ x_a_p(t)) .* (t < t_0) + ...
    (x_d_p(t) ./ x_b_p(t)) .* (t_0 < t);
% ratio((abs(ratio) > 1)) = nan;
plot(t, ratio,...
        'LineWidth', 2,...
        'Color', 'b');

set(gca, 'FontSize', 16)
xlabel('Time (hours)','FontSize',18);
% ylabel('%ID/g','FontSize',18);
% ylabel('Cells / 100 uL','FontSize',18);
% ylabel('$$x''_T(t) / x''_L(t)$$','FontSize',18, 'Interpreter', 'latex');
ylabel('x''_T(t) / x''_L(t)','FontSize',18);
 

% % % Plot for the ratio of concentration 
figure();
ratio = y(:,2) ./ y(:,1);
plot(t, ratio,...
        'LineWidth', 2,...
        'Color', 'b');

set(gca, 'FontSize', 16)
xlabel('Time (hours)','FontSize',18);
% ylabel('%ID/g','FontSize',18);
% ylabel('Cells / 100 uL','FontSize',18);
% ylabel('$$x''_T(t) / x''_L(t)$$','FontSize',18, 'Interpreter', 'latex');
ylabel('x_T(t) / x_L(t)','FontSize',18);



% % % Plot for the derivative of the ratio of concentration 
figure();
% ratio((abs(ratio) > 1)) = nan;
plot(t, r_fun(t),...
        'LineWidth', 2,...
        'Color', 'b');

set(gca, 'FontSize', 16)
xlabel('Time (hours)','FontSize',18);
% ylabel('%ID/g','FontSize',18);
% ylabel('Cells / 100 uL','FontSize',18);
% ylabel('$$x''_T(t) / x''_L(t)$$','FontSize',18, 'Interpreter', 'latex');
ylabel('[x_T(t) / x_L(t)]''','FontSize',18);
 
%%%================================================================

function compute_residual

% Compute manually the residual norm
[t,y] = ode15s(f, [0, time_vector], initial_conditions);
z = y - [zeros(1, size(expected_data,2)); expected_data];
m_error = norm(z,2);
z = z(:);
error = norm(z,2);

disp(['L2 error: ', num2str(error)]);
disp(['L2 m_error: ', num2str(m_error)]);

%%%================================================================

function plot_mouse_data(s)

% Choose which type of manipulation is going to be
% applied to the data

group          = 'control';
mouse_selector = [1,2,3];
n_mice         = length(mouse_selector);
% clusters = {{'lnol','lnor'},{'tl','tr'}};
% clusters = {{'lnol'},{'tl'}};
% clusters = {{'lnor'},{'tr'}};

clusters = {'lnol','lnor','tl','tr'};

n_clusters = length(clusters);
time_vector = s.times;
% tissues = fieldnames(s.(group));
n_times   = length(time_vector);

marker_vector = {'s','o','x'};
color_vector = {'r','k','b'};



for tissue_index = 1:n_clusters
    tissue = clusters{tissue_index};
    data = zeros(n_times, n_mice);
    
    figure();
    hold on;
    
    for time_index = 1:n_times
        for mouse_index = 1:n_mice
            y = s.(group).(tissue).raw_data{time_index}(mouse_index);
            data(time_index, mouse_index) = y;
        end
    end
    
    for mouse_index = 1:n_mice 
        plot(time_vector, data(:,mouse_index),...
            'Marker', marker_vector{mouse_index},...
            'LineStyle', '-', ...
            'Color', color_vector{mouse_index},...
            'MarkerSize',6,...
            'LineWidth', 2);
    end
    set(gca, 'FontSize', 16)
    xlabel('Time (hours)','FontSize',18);
    ylabel('Cells / \muL','FontSize',18);
    legend('M1','M2', 'M3','location','northwest')
    title(tissue)
    
end

























