%% Load whisking and neural time series struct
clear
load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\Raw\excitatory_clean.mat') %L5b excitatory cells recorded by Jon and Phil

data_directory = 'C:\Users\jacheung\Dropbox\LocationCode\DataStructs\';

load([data_directory 'Raw\excitatory_clean.mat']);
feature_list = {'angle','phase'};
whisk_struct = CompileWStruct(data_directory, feature_list);
touch_struct = CompileTStruct(data_directory, feature_list);
%% Plots
% proportion of neurons that intersect pie (A)
hilbert_var = 'phase'; 
pie_comparison(touch_struct,whisk_struct,hilbert_var)

% whisk and touch tuning curves (B)
hilbert_var = 'angle'; 
plot_tuning_curves(touch_struct,whisk_struct,hilbert_var) 

% Whisk x touch modulation depth (C)
hilbert_var = 'phase'; 
population_modulation_comparison(touch_struct,whisk_struct,hilbert_var)

% Shape correlation of tuning curves (D)
hilbert_var = 'angle'; 
intersect_correlation_v2(touch_struct,whisk_struct,hilbert_var); 

% scatter of tuning preference of whisk and touch (E)
hilbert_var = 'angle';
tuning_preference_scatter(touch_struct,whisk_struct,hilbert_var)

%% Heatmap (SF)
population_heatmap_builder(touch_struct.(hilbert_var),whisk_struct.(hilbert_var),hilbert_var)

saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig6\';
figure(50)
fn = [hilbert_var '_intersect_heat_raw.eps'];
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])

figure(51)
fn = [hilbert_var '_intersect_heat_squish.eps'];
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])

figure(52)
fn = [hilbert_var '_histograms.eps'];
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])
