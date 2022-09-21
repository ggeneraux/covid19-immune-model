%% Run U1 - Nominal Response

% load parameters

data_dictionary = get_data_dictionary();
p = data_dictionary.parameters; % obtain parameters
%%% base case - no change to parameters 


%%%% [parameter struct, viral load, damaged cells]
% turn off Ab production
p.tau_Ab = 1e9; % antibody molecules/mL (set to initial value but neutralizing effect currently off in model)
p.k_Ab_V = 0; 
data_dictionary.parameters = p;
[T,X] = function_run_model(data_dictionary,1e6,0,1);

S = table;
S.time = T;

for i=1:data_dictionary.n_species
    S.(data_dictionary.species_names(i,2)) = X(:,i);
end


%% Run U2 - No innate Response

data_dictionary = get_data_dictionary();
p = data_dictionary.parameters; % obtain parameters

% No innate immunity
 
p.a_ifnb_i = 0;
p.a_ifnb_d = 0;
p.a_ifnb_dc = 0;
% virus and infected cell activation of DC and M1 turned off
p.k_v = 0.0; ...
p.k_I = 0.0; ...
p.k_dAT = 0.0; ...
p.a_Ab = 0;
p.k_N_IFNg = 0;
p.k_N_TNFa = 0;
p.k_N_GMCSF = 0;
p.k_N_IL17c = 0;
p.a_il6_i = 0;
p.a_tnf_i = 0;
% turn off Ab production
p.tau_Ab = 1e9;
p.k_Ab_V = 0; 

data_dictionary.parameters=p;

[T,X] = function_run_model(data_dictionary,1e6,0,2);

S = table;
S.time = T;

for i=1:data_dictionary.n_species
    S.(data_dictionary.species_names(i,2)) = X(:,i);
end
 
 
%% Run U3 - Sterile Immune Response
 
data_dictionary = get_data_dictionary();
p = data_dictionary.parameters; % obtain parameters

% turn off Ab production
p.tau_Ab = 1e9;
p.k_Ab_V = 0; 
[T,X] = function_run_model(data_dictionary,0,1.5e8,3);

S = table;
S.time = T;

for i=1:data_dictionary.n_species
    S.(data_dictionary.species_names(i,2)) = X(:,i);
end

%% Run U4 - Sustained Inflammation

data_dictionary = get_data_dictionary();
p = data_dictionary.parameters; % obtain parameters
p_original = p;

% turn off Ab production
p.tau_Ab = 1e9;
p.k_Ab_V = 0; 
damage_sensitivity = 12; %Increase sensitivity of inflammatory cell death
p.k_damage_TNFa = p_original.k_damage_TNFa * damage_sensitivity;
p.k_damage_IL6 = p_original.k_damage_IL6 * damage_sensitivity;
p.k_damage_IL1b = p_original.k_damage_IL1b * damage_sensitivity;
p.k_damage_IFNg = p_original.k_damage_IFNg * damage_sensitivity;
data_dictionary.parameters=p;
 
[T,X] =  function_run_model(data_dictionary,1e6,0,4);

S = table;
S.time = T;

for i=1:data_dictionary.n_species
    S.(data_dictionary.species_names(i,2)) = X(:,i);
end

%% Calculate algebraic rates for inspection post-simulation


time_index = 1:height(S);

rates_width = data_dictionary.n_rates;
rates_length = length(time_index);

rates = zeros(rates_length, rates_width);

%% Calculate algebraic rates

for i = 1:rates_length
    rates(i,:) = calculate_rates(S(time_index(i),:),data_dictionary);
end


%% Create table for rates

R = table;
R.time = T(time_index,1);

for i=1:data_dictionary.n_rates
    R.(data_dictionary.rates_names(i,2)) = rates(:,i);
end

%% figure formating
 figuremods();