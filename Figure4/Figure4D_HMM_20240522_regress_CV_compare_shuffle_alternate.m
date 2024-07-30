
function Figure4D_HMM_20240522_regress_CV_compare_shuffle_alternate

[OFC_p, OFC_neuron] = temp_HMM_20240522_regress_CV_compare_shuffle('altern_OFC_20240427');
[Hippo_p, Hippo_neuron] = temp_HMM_20240522_regress_CV_compare_shuffle('altern_Hippo_20240427');
[AC_p, AC_neuron] = temp_HMM_20240522_regress_CV_compare_shuffle('altern_AC_20240427');
[PPC_p, PPC_neuron] = temp_HMM_20240522_regress_CV_compare_shuffle('altern_PPC_20240427');
[M1_p, M1_neuron] = temp_HMM_20240522_regress_CV_compare_shuffle('altern_M1_20240427');
[STR_p, STR_neuron] = temp_HMM_20240522_regress_CV_compare_shuffle('altern_STR_20240427');
close all

% length_neuron = [sig_before_sound, N_before_sound;
%                 sig_all_sound,     N_all_sound;
%                 sig_choice2,       N_choice2];
disp('plot before_sound')
plot_all_zigzag(1, OFC_p,OFC_neuron,Hippo_p,Hippo_neuron,AC_p,AC_neuron,PPC_p,PPC_neuron,M1_p,M1_neuron,STR_p,STR_neuron);
title('Before sound')
disp('plot sound')
plot_all_zigzag(2, OFC_p,OFC_neuron,Hippo_p,Hippo_neuron,AC_p,AC_neuron,PPC_p,PPC_neuron,M1_p,M1_neuron,STR_p,STR_neuron);
title('During sound')
disp('plot choice')
plot_all_zigzag(3, OFC_p,OFC_neuron,Hippo_p,Hippo_neuron,AC_p,AC_neuron,PPC_p,PPC_neuron,M1_p,M1_neuron,STR_p,STR_neuron);
title('During outcome')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_all_zigzag(number, OFC_p,OFC_neuron,Hippo_p,Hippo_neuron,AC_p,AC_neuron,PPC_p,PPC_neuron,M1_p,M1_neuron,STR_p,STR_neuron)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y_ofc = (1/OFC_neuron(number,2):1/OFC_neuron(number,2):1)';
y_Hippo = (1/Hippo_neuron(number,2):1/Hippo_neuron(number,2):1)';
y_AC = (1/AC_neuron(number,2):1/AC_neuron(number,2):1)';
y_PPC = (1/PPC_neuron(number,2):1/PPC_neuron(number,2):1)';
y_M1 = (1/M1_neuron(number,2):1/M1_neuron(number,2):1)';
y_STR = (1/STR_neuron(number,2):1/STR_neuron(number,2):1)';


OFC = sort(OFC_p(number).matrix); 
Hippo = sort(Hippo_p(number).matrix); 
AC = sort(AC_p(number).matrix); 
PPC = sort(PPC_p(number).matrix); 
M1 = sort(M1_p(number).matrix); 
STR = sort(STR_p(number).matrix); 

OFC = remove_inf_OFC(OFC);
Hippo = remove_inf_OFC(Hippo);
AC = remove_inf_OFC(AC);
PPC = remove_inf_OFC(PPC);
M1 = remove_inf_OFC(M1);
STR = remove_inf_OFC(STR);

thre = -log10(0.025);

figure
plot(OFC,y_ofc,'r')
hold on
plot(Hippo,y_Hippo,'g')
hold on
plot(AC,y_AC,'b')
hold on
plot(PPC,y_PPC,'c')
hold on
plot(M1,y_M1,'m')
hold on
plot(STR,y_STR,'y')
hold on
plot([thre,thre],[0,1],'k')

disp([OFC_neuron(number,1:2); Hippo_neuron(number,1:2); AC_neuron(number,1:2); ...
    PPC_neuron(number,1:2); M1_neuron(number,1:2); STR_neuron(number,1:2)])

disp([OFC_neuron(number,3); Hippo_neuron(number,3); AC_neuron(number,3); ...
    PPC_neuron(number,3); M1_neuron(number,3); STR_neuron(number,3)])

return




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function OFC = remove_inf_OFC(OFC)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
temp = find(OFC == inf);
temp_OFC = OFC;
temp_OFC(temp) = [];
max_OFC = max(temp_OFC);
OFC(temp) = max_OFC;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p_values, length_neuron] = temp_HMM_20240522_regress_CV_compare_shuffle(folders)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch nargin
    case 0
        hoge
    case 1
        disp('OK to analyze')
    otherwise
        hoge
end

[analysis_dir,depth_def] = eval(folders);

prob_before_sound = [];
prob_all_sound = [];
prob_choice2 = [];
for i = 1:length(analysis_dir)
    disp([i,length(analysis_dir)])
    temp_dir = analysis_dir{i};
    
    [temp_before_sound,temp_all_sound,temp_choice] = ...
        HMM_ephys_20240522_regress_get_Data_shuffle(temp_dir, depth_def);

    prob_before_sound = [prob_before_sound; temp_before_sound];
    prob_all_sound = [prob_all_sound; temp_all_sound];
    prob_choice2 = [prob_choice2; temp_choice];
end
delete(gcp('nocreate'))

%check whether the STD = 1.96 is about 0.025
temp = find(prob_before_sound(:,1) > 0.01);
prob_before_sound(temp,:)
figure
plot(prob_before_sound(temp,1),prob_before_sound(temp,2),'k.')

figure
subplot(1,3,1)
plot(prob_before_sound(:,1), prob_before_sound(:,2), 'b.')
subplot(1,3,2)
plot(prob_all_sound(:,1), prob_all_sound(:,2), 'b.')
subplot(1,3,3)
plot(prob_choice2(:,1), prob_choice2(:,2), 'b.')

figure
subplot(1,3,1)
plot(-log10(prob_before_sound(:,1)), prob_before_sound(:,2), 'r.')
subplot(1,3,2)
plot(-log10(prob_all_sound(:,1)), prob_all_sound(:,2), 'r.')
subplot(1,3,3)
plot(-log10(prob_choice2(:,1)), prob_choice2(:,2), 'r.')
%convert to p_values
% sum(prob_before_sound)
% sum(prob_all_sound)
% sum(prob_choice2)

[N_before_sound,~] = size(prob_before_sound);
[N_all_sound,~] = size(prob_all_sound);
[N_choice2,~] = size(prob_choice2);

%We need to plot the -log10(:,1)
p_values(1).matrix = -log10(prob_before_sound(:,1));
p_values(2).matrix = -log10(prob_all_sound(:,1));
p_values(3).matrix = -log10(prob_choice2(:,1));

thre = 0.025;
sig_before_sound = length(find(p_values(1).matrix > -log10(thre)));
sig_all_sound = length(find(p_values(2).matrix > -log10(thre)));
sig_choice2 = length(find(p_values(3).matrix > -log10(thre)));

length_neuron = [sig_before_sound, N_before_sound, sig_before_sound/N_before_sound;
                sig_all_sound,     N_all_sound,    sig_all_sound/N_all_sound;
                sig_choice2,       N_choice2,      sig_choice2/N_choice2];


return


function [prob_before_sound,prob_sound,prob_choice] = ...
    HMM_ephys_20240522_regress_get_Data_shuffle(pathname, depth_def)

switch nargin
    case 0
        pathname = pwd;
    case 2
        disp('OK to analyze')
    otherwise
        hoge
end
cd(pathname)

% [filename1, pathname1,findex]=uigetfile('*.*','frame file');
% filename1 = [pathname1,filename1];
% load(filename1)
temp = dir('HMM_spike_count_neurons_20230316*');
if length(temp) ~= 1
    disp(temp)
    hoge
end
load(temp.name);
%neuron_index p_index

temp = dir('sig_HMM_neurons_20230310*');
if length(temp) ~= 1
    disp(temp)
    hoge
end
load(temp.name);
%p_task: around sound
%p_task2: around choice

temp = dir('HMM_ephys_20240520_regress_CV_analysis*');
if length(temp) ~= 1
    disp(temp)
    hoge
end
% all_MS ...
% no_choice_MS no_pre_choice_MS no_pre_reward_MS ...
% no_pre_sound_MS no_reward_MS no_run_MS no_sound_MS shuffle_MS
load(temp.name);

temp = dir('depth_spike_20230427*');
if length(temp) == 1
    load(temp.name);
    %spike_depth def_depth length_neuron
    if size(p_task,1) ~= length_neuron
        disp([length(p_task) length_neuron])
        hoge
    end
    
    if depth_def == 1
        depth_neuron = find(spike_depth <= def_depth(1));
    else
        depth_neuron = find(spike_depth > def_depth(1) & spike_depth <= def_depth(2));
    end
elseif isempty(temp)
    depth_neuron = 1:size(p_task,1); %Use all the neurons
    hoge
else
    disp(temp)
    hoge
end

% %Get sig neurons
new_p_thre = 10;

if ~exist('shuffle_MS', 'var')
    prob_before_sound = [];
    prob_sound = [];
    prob_choice = [];
else
  
    sig_before_sound = get_sig_neuron_time_window(p_task,new_p_thre,10:15); %600ms before sound
    sig_during_sound = get_sig_neuron_time_window(p_task,new_p_thre,16:21); %During sound
    sig_during_choice = get_sig_neuron_time_window(p_task2,new_p_thre,6:15); %During choice 1000ms

    %ADD depth definition
    sig_before_sound = intersect(sig_before_sound,depth_neuron);
    sig_during_sound = intersect(sig_during_sound,depth_neuron);
    sig_during_choice = intersect(sig_during_choice,depth_neuron);

    thre = 1.96;
    %Get the sig neurons and the regression analysis
    prob_before_sound = get_CV_std_shuffle(sig_before_sound, 1, thre, all_MS, shuffle_MS);

    prob_sound = get_CV_std_shuffle(sig_during_sound, 2, thre, all_MS, shuffle_MS);

    prob_choice = get_CV_std_shuffle(sig_during_choice, 3, thre, all_MS, shuffle_MS);
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function prob_gauss = get_CV_std_shuffle(sig_before_sound, use_number, thre, all_MS, shuffle_MS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(use_number)
[all_mean, ~] = get_ave_thre_shuffle(all_MS,sig_before_sound,use_number, thre);
[s_mean, s_std, ~] = get_ave_thre_shuffle(shuffle_MS,sig_before_sound,use_number, thre);

sabun = s_mean - all_mean;
sabun = sabun ./ s_std;

prob_gauss = normcdf(all_mean, s_mean, s_std);
prob_gauss(:,2) = sabun;

% figure
% plot(prob_gauss(:,1), prob_gauss(:,2), 'b.')
% %convert to p_values

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [all_mean, all_std, neuron_thre] = get_ave_thre_shuffle(all_MS,sig_before_sound,use_number, thre)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
all_mean = all_MS(1).matrix(sig_before_sound,use_number);
all_std  = all_MS(2).matrix(sig_before_sound,use_number);
neuron_thre = all_mean - thre * all_std;

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sig_neuron = get_sig_neuron_time_window(p_task,new_p_thre,window)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
use_p_sound = min(p_task(:,window),[],2);
use_p_sound = -log10(use_p_sound);
sig_neuron = find(use_p_sound > new_p_thre);

return


