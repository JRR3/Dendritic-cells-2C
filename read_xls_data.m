%%%Source        : Houston Methodist Research Institute
%%%Location      : Houston, TX.
%%%Origin        : September 28, 2018
%%%PI            : Mauro Ferrari
%%%Supervisor    : Vittorio Cristini
%%%Collaborator  : Prashant Dogra
%%%Developer     : Javier Ruiz Ramirez

%%%================================================================
function S = read_xls_data()
clc;
format compact;
path = '.';
fname = 'mouse_data.xlsx';
fname = [path,'\',fname];
[~,~,raw] = xlsread(fname);
%---------------------------Default variables
list_of_group_symbols = {'control','gene'};
list_of_tissues       = {};
%---------------------------Counters
s           = size(raw);
nrows       = s(1);
ncols       = s(2);
current_row = 1;
%---------------------------Variables
flag        = false;
%---------------------------keywords
time_keyword    = 'time';
control_keyword = 'control';

list_of_group_names = {control_keyword, 'gene'};
map_group_symbol_to_name = ...
    containers.Map('KeyType','char','ValueType','char');

%In case the notation used for the group symbol contains
for k = 1:length(list_of_group_symbols)
    symbol = list_of_group_symbols{k};
    map_group_symbol_to_name(symbol) = list_of_group_names{k};
end

%---------------------------Containers
S             = [];
time_vector   = [];
group         = [];
tissue        = [];
temporal_data = [];

%---------------------------Read time
%Assumption: Time keyword is followed by time value.
while isempty(time_vector)
    for k = 1:ncols
        cell_field = lower(raw{current_row, k});
        if flag == true
            time_vector = [time_vector, cell_field];
            flag = false;
        end
        if strcmp(cell_field,time_keyword)
            flag = true;
        end
    end
    current_row = current_row + 1;
end
%---------------------------Save time vector
S.times = time_vector;

%---------------------------Read groups
while current_row <= nrows
    %disp(['Currently working on row: ', num2str(current_row)]);
    for k = 1:ncols
        %disp(['Currently working on col: ', num2str(k)]);
        cell_field = lower(raw{current_row, k});
        
        %Is this a new group?
        if string_is_inside(cell_field,list_of_group_symbols)
            group = map_group_symbol_to_name(cell_field);
            S.(group) = [];
            break;
        end
        
        %Is this a new tissue?
        %if string_is_inside(cell_field, list_of_tissues)
        [logic_value, label] = is_this_a_new_tissue(cell_field);
        if logic_value    
            tissue = label;
            temporal_data = [];
            S.(group).(tissue).raw_data = [];
            S.(group).(tissue).mean     = [];
            S.(group).(tissue).std      = [];
            continue;
        end
        
        %Is this new data?
        if ~isnan(cell_field)
            temporal_data(end+1) = cell_field;
        else
            if isempty(temporal_data)
                continue;
            else
                S = store_data_in_structure(S,group,tissue,temporal_data);
                temporal_data = [];
            end
        end
        
    end
    
    if ~isempty(temporal_data)
        S = store_data_in_structure(S,group,tissue,temporal_data);
        temporal_data = [];
    end
    current_row = current_row + 1;
end

%Compute significance
significance_value = 0.05;
control = map_group_symbol_to_name(control_keyword);
for i = 2:length(list_of_group_names)
    group = list_of_group_names{i};
    
    if ~isfield(S, group)
        continue;
    end
    
    for k = 1:length(list_of_tissues)
        tissue = list_of_tissues{k};
        
        if ~isfield(S.(group), tissue)
            continue;
        end
        
        for r = 1:length(S.(group).(tissue).raw_data)
            x_sample = S.(control).(tissue).raw_data{r};
            x_mean   = S.(control).(tissue).mean(r);
            y_sample = S.(group).(tissue).raw_data{r};
            y_mean   = S.(group).(tissue).mean(r);
            %Assuming homoscedasticity
            [h,p] = ttest2(x_sample, y_sample,...
                'VarType', 'equal',...
                'Alpha', significance_value, ...
                'Tail', 'both');
            S.(group).(tissue).h_value(r) = h;
            S.(group).(tissue).p_value(r) = p;
            S.(group).(tissue).ratio_to_control(r)   = y_mean / x_mean;
            if h == 0
                disp('Data was not significant in:');
                disp(['Group : ', group]);
                disp(['Tissue: ', tissue]);
                disp(['Index : ', num2str(r)]);
                disp(['Time  : ', num2str(time_vector(r))]);
                disp(['p val.: ', num2str(p)]);
                disp(['--------------------']);
            end
        end
        %Ratio of means
    end
end
%End of main program

%%%================================================================

function r = string_is_inside(y,L)

r = false;

for k = 1:length(L)
    if strcmp(y,L{k})
        r = true;
        return;
    end
end

%%%================================================================

function s = store_data_in_structure(s,group,tissue,temporal_data)
s.(group).(tissue).raw_data{end+1} = temporal_data;
s.(group).(tissue).mean(end+1) = mean(temporal_data);
s.(group).(tissue).std(end+1) = std(temporal_data);

%%%================================================================

function [logic_value, label] = is_this_a_new_tissue(cell_field)

logic_value = true;

if sum(isnan(cell_field)) || isa(cell_field, 'double')
    logic_value = false;
    label = [];
    return;
end

if ~isempty(regexpi(cell_field, '([.0-9])|(control)', 'match'))
    logic_value = false;
    label = [];
    return;
end

regexp = '[+- -]';
rep    = 'o';
label  = regexprep(cell_field, regexp, rep);











