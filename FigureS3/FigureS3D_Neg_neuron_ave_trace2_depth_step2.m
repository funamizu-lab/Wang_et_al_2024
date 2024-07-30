function FigureS3D_Neg_neuron_ave_trace2_depth_step2

[~,norm_neg_all_OFC_rep, propor_rep(1,:),num_rep(1,:),number_pos_OFC_rep,number_neg_OFC_rep] = Neg_neuron_ave_trace2_depth_step1_1('repeat_OFC_20240427');
[~,norm_neg_all_HPC_rep, propor_rep(2,:),num_rep(2,:),number_pos_HPC_rep,number_neg_HPC_rep] = Neg_neuron_ave_trace2_depth_step1_1('repeat_Hippo_20240427');
[~,norm_neg_all_AC_rep, propor_rep(3,:),num_rep(3,:),number_pos_AC_rep,number_neg_AC_rep] = Neg_neuron_ave_trace2_depth_step1_1('repeat_AC_20240427');
[~,norm_neg_all_PPC_rep, propor_rep(4,:),num_rep(4,:),number_pos_PPC_rep,number_neg_PPC_rep] = Neg_neuron_ave_trace2_depth_step1_1('repeat_PPC_20240427');

[~,norm_neg_all_OFC_zig, propor_zig(1,:),num_zig(1,:),number_pos_OFC_zig,number_neg_OFC_zig] = Neg_neuron_ave_trace2_depth_step1_1('altern_OFC_20240427');
[~,norm_neg_all_HPC_zig, propor_zig(2,:),num_zig(2,:),number_pos_HPC_zig,number_neg_HPC_zig] = Neg_neuron_ave_trace2_depth_step1_1('altern_Hippo_20240427');
[~,norm_neg_all_AC_zig, propor_zig(3,:),num_zig(3,:),number_pos_AC_zig,number_neg_AC_zig] = Neg_neuron_ave_trace2_depth_step1_1('altern_AC_20240427');
[~,norm_neg_all_PPC_zig, propor_zig(4,:),num_zig(4,:),number_pos_PPC_zig,number_neg_PPC_zig] = Neg_neuron_ave_trace2_depth_step1_1('altern_PPC_20240427');
[~,norm_neg_all_M1_zig, propor_zig(5,:),num_zig(5,:),number_pos_M1_zig,number_neg_M1_zig] = Neg_neuron_ave_trace2_depth_step1_1('altern_M1_20240427');
[~,norm_neg_all_STR_zig, propor_zig(6,:),num_zig(6,:),number_pos_STR_zig,number_neg_STR_zig] = Neg_neuron_ave_trace2_depth_step1_1('altern_STR_20240427');


OFC_rep = nan(length(number_pos_OFC_rep),2);
OFC_rep(:,1) = number_pos_OFC_rep(:,1)./ number_pos_OFC_rep(:,2);
OFC_rep(:,2) = number_neg_OFC_rep(:,1)./ number_neg_OFC_rep(:,2);


HPC_rep = nan(length(number_pos_HPC_rep),2);
HPC_rep(:,1) = number_pos_HPC_rep(:,1)./ number_pos_HPC_rep(:,2);
HPC_rep(:,2) = number_neg_HPC_rep(:,1)./ number_neg_HPC_rep(:,2);

AC_rep = nan(length(number_pos_AC_rep),2);
AC_rep(:,1) = number_pos_AC_rep(:,1)./ number_pos_AC_rep(:,2);
AC_rep(:,2) = number_neg_AC_rep(:,1)./ number_neg_AC_rep(:,2);

PPC_rep = nan(length(number_pos_PPC_rep),2);
PPC_rep(:,1) = number_pos_PPC_rep(:,1)./ number_pos_PPC_rep(:,2);
PPC_rep(:,2) = number_neg_PPC_rep(:,1)./ number_neg_PPC_rep(:,2);

OFC_zig = nan(length(number_pos_OFC_zig),2);
OFC_zig(:,1) = number_pos_OFC_zig(:,1)./ number_pos_OFC_zig(:,2);
OFC_zig(:,2) = number_neg_OFC_zig(:,1)./ number_neg_OFC_zig(:,2);

HPC_zig = nan(length(number_pos_HPC_zig),2);
HPC_zig(:,1) = number_pos_HPC_zig(:,1)./ number_pos_HPC_zig(:,2);
HPC_zig(:,2) = number_neg_HPC_zig(:,1)./ number_neg_HPC_zig(:,2);

AC_zig = nan(length(number_pos_AC_zig),2);
AC_zig(:,1) = number_pos_AC_zig(:,1)./ number_pos_AC_zig(:,2);
AC_zig(:,2) = number_neg_AC_zig(:,1)./ number_neg_AC_zig(:,2);


PPC_zig = nan(length(number_pos_PPC_zig),2);
PPC_zig(:,1) = number_pos_PPC_zig(:,1)./ number_pos_PPC_zig(:,2);
PPC_zig(:,2) = number_neg_PPC_zig(:,1)./ number_neg_PPC_zig(:,2);

M1_zig = nan(length(number_pos_M1_zig),2);
M1_zig(:,1) = number_pos_M1_zig(:,1)./ number_pos_M1_zig(:,2);
M1_zig(:,2) = number_neg_M1_zig(:,1)./ number_neg_M1_zig(:,2);

STR_zig = nan(length(number_pos_STR_zig),2);
STR_zig(:,1) = number_pos_STR_zig(:,1)./ number_pos_STR_zig(:,2);
STR_zig(:,2) = number_neg_STR_zig(:,1)./ number_neg_STR_zig(:,2);

rep_prop = [OFC_rep;HPC_rep;AC_rep;PPC_rep];
zig_prop = [OFC_zig;HPC_zig;AC_zig;PPC_zig;M1_zig;STR_zig];


p_rep = signrank(rep_prop(:,1),rep_prop(:,2));
p_zig = signrank(zig_prop(:,1),zig_prop(:,2));

p_rep_OFC = signrank(OFC_rep(:,1),OFC_rep(:,2));
p_rep_HPC = signrank(HPC_rep(:,1),HPC_rep(:,2));
p_rep_AC = signrank(AC_rep(:,1),AC_rep(:,2));
p_rep_PPC = signrank(PPC_rep(:,1),PPC_rep(:,2));

p_rep_all = [p_rep_OFC,p_rep_HPC,p_rep_AC,p_rep_PPC];

p_zig_OFC = signrank(OFC_zig(:,1),OFC_zig(:,2));
p_zig_HPC = signrank(HPC_zig(:,1),HPC_zig(:,2));
p_zig_AC = signrank(AC_zig(:,1),AC_zig(:,2));
p_zig_PPC = signrank(PPC_zig(:,1),PPC_zig(:,2));
p_zig_M1 = signrank(M1_zig(:,1),M1_zig(:,2));
p_zig_STR = signrank(STR_zig(:,1),STR_zig(:,2));

p_zig_all = [p_zig_OFC,p_zig_HPC,p_zig_AC,p_zig_PPC,p_zig_M1,p_zig_STR];


figure
plot_mean_se_moto(norm_neg_all_OFC_rep,'k',2)
hold on
plot_mean_se_moto(norm_neg_all_HPC_rep,[0.9290 0.6940 0.1250],2)
hold on
plot_mean_se_moto(norm_neg_all_AC_rep,[0.4660 0.6740 0.1880],2)
hold on
plot_mean_se_moto(norm_neg_all_PPC_rep,[0 0.4470 0.7410],2)
hold on

figure
plot_mean_se_moto(norm_neg_all_OFC_zig,'k',2)
hold on
plot_mean_se_moto(norm_neg_all_HPC_zig,[0.9290 0.6940 0.1250],2)
hold on
plot_mean_se_moto(norm_neg_all_AC_zig,[0.4660 0.6740 0.1880],2)
hold on
plot_mean_se_moto(norm_neg_all_PPC_zig,[0 0.4470 0.7410],2)
hold on
plot_mean_se_moto(norm_neg_all_M1_zig,[0.6350 0.0780 0.1840],2)
hold on
plot_mean_se_moto(norm_neg_all_STR_zig,[0.8500 0.3250 0.0980],2)
hold on

return


function [norm_pos_all,norm_neg_all, proportion,neuron_number,neuron_number_pos,neuron_number_neg] = ...
    Neg_neuron_ave_trace2_depth_step1_1(folders)

[analysis_dir,depth_def] = eval(folders);

norm_pos_all = [];
norm_neg_all = [];

neuron_number_pos = nan(length(analysis_dir),3);
neuron_number_neg = nan(length(analysis_dir),3);

for i = 1:length(analysis_dir)
    disp([i,length(analysis_dir)])
    temp_dir = analysis_dir{i};
    
    [sig_norm_spike_pos,~,~,neuron_number_pos(i,:),~,...
   sig_norm_spike_neg,~,~,neuron_number_neg(i,:),~] = ...
        Neg_neuron_ave_trace2_depth_step(temp_dir,depth_def);

    norm_pos_all = [norm_pos_all; sig_norm_spike_pos];

    norm_neg_all = [norm_neg_all; sig_norm_spike_neg];
end
delete(gcp('nocreate'))

porportion_pos = sum(neuron_number_pos(:,1))/sum(neuron_number_pos(:,2));
porportion_neg = sum(neuron_number_neg(:,1))/sum(neuron_number_neg(:,2));

proportion = [porportion_pos, porportion_neg];
neuron_number = [sum(neuron_number_pos(:,1)),sum(neuron_number_neg(:,1))];

b = ['Negative neuron: ' num2str(porportion_neg*100) '%'];

disp(b);


return



function [sig_norm_spike_pos,sig_spike01_pos,sig_max_time_pos,neuron_number_pos,per_sig_neuron_pos,...
   sig_norm_spike_neg,sig_spike01_neg,sig_max_time_neg,neuron_number_neg,per_sig_neuron_neg] = ...
    Neg_neuron_ave_trace2_depth_step(pathname,depth_def)

switch nargin
    case 0
        pathname = pwd;
    case 1
        hoge
    case 2
        disp('OK to analyze')
    otherwise
        hoge
end
cd(pathname)

% [filename1, pathname1,findex]=uigetfile('*.*','frame file');
% filename1 = [pathname1,filename1];
% load(filename1)
temp = dir('HMM_ephys1_20230315_make_ave_trace*');
if length(temp) ~= 1
    hoge
end
load(temp.name);
%norm_spike max_time spike01 mean_spike std_spike

temp = dir('sig_neg_task_neurons_20240417*');
if length(temp) ~= 1
    hoge
end
load(temp.name);
%p_task: around sound
%p_task2: around choice

temp = dir('depth_spike_20230427*');
if length(temp) == 1
    load(temp.name);
    %spike_depth def_depth length_neuron
    if size(p_task_neg,1) ~= length_neuron
        disp([length(p_task_neg) length_neuron])
        hoge
    end
    
    if depth_def == 1
        depth_neuron = find(spike_depth <= def_depth(1));
    else
        depth_neuron = find(spike_depth > def_depth(1) & spike_depth <= def_depth(2));
    end
elseif isempty(temp)
    depth_neuron = 1:size(p_task_neg,1); %Use all the neurons
    hoge
else
    hoge
end

%Get sig neurons
[size_neuron,~] = size(p_task_neg);
new_p_thre = 10;
sigsound_pos = get_sig_neuron_all(p_task_pos,new_p_thre);
sigchoice_pos = get_sig_neuron_all(p_task2_pos,new_p_thre);
sig_neuron_pos = union(sigsound_pos,sigchoice_pos);


sigsound_neg = get_sig_neuron_all(p_task_neg,new_p_thre);
sigchoice_neg = get_sig_neuron_all(p_task2_neg,new_p_thre);
sig_neuron_neg = union(sigsound_neg,sigchoice_neg);

sig_neuron_neg = setdiff(sig_neuron_neg, sig_neuron_pos);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ADD depth definition
sig_neuron_pos = intersect(sig_neuron_pos,depth_neuron);
sig_neuron_neg = intersect(sig_neuron_neg,depth_neuron);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
neuron_number_neg = [length(sig_neuron_neg), length(depth_neuron), size_neuron];
neuron_number_pos = [length(sig_neuron_pos), length(depth_neuron), size_neuron];
%neuron_number = [size_neuron, length(sig_neuron)];

sig_norm_spike_neg = norm_spike.all(sig_neuron_neg,:);

[~,time_trace] = size(sig_norm_spike_neg);
max_sig_trace = max(sig_norm_spike_neg,[],2);
min_sig_trace = min(sig_norm_spike_neg,[],2);

max_sig_trace = repmat(max_sig_trace,1,time_trace);
min_sig_trace = repmat(min_sig_trace,1,time_trace);
sig_spike01_neg = (sig_norm_spike_neg - min_sig_trace) ./ (max_sig_trace-min_sig_trace);

%sig_spike01 = spike01.all(sig_neuron,:);
sig_max_time_neg = max_time.all(sig_neuron_neg,:);


sig_norm_spike_pos = norm_spike.all(sig_neuron_pos,:);

[~,time_trace] = size(sig_norm_spike_pos);
max_sig_trace = max(sig_norm_spike_pos,[],2);
min_sig_trace = min(sig_norm_spike_pos,[],2);

max_sig_trace = repmat(max_sig_trace,1,time_trace);
min_sig_trace = repmat(min_sig_trace,1,time_trace);
sig_spike01_pos = (sig_norm_spike_pos - min_sig_trace) ./ (max_sig_trace-min_sig_trace);

%sig_spike01 = spike01.all(sig_neuron,:);
sig_max_time_pos = max_time.all(sig_neuron_pos,:);
per_sig_neuron_pos = neuron_number_pos(1)/neuron_number_pos(2);

per_sig_neuron_neg = neuron_number_neg(1)/neuron_number_neg(2);
 
 
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sig_neuron = get_sig_neuron_all(p_task,new_p_thre)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
use_p_sound = min(p_task,[],2);
use_p_sound = -log10(use_p_sound);
sig_neuron = find(use_p_sound > new_p_thre);

return








