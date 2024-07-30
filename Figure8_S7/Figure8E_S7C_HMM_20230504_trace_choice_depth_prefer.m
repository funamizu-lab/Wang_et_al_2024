

function Figure8E_S7C_HMM_20230504_trace_choice_depth_prefer

%repeat
process_HMM_20230504_trace_at_choice_ver2_depth_correct('repeat_OFC_20240427');
process_HMM_20230504_trace_at_choice_ver2_depth_correct('repeat_Hippo_20240427');
process_HMM_20230504_trace_at_choice_ver2_depth_correct('repeat_AC_20240427');
process_HMM_20230504_trace_at_choice_ver2_depth_correct('repeat_PPC_20240427');

%zigzag
process_HMM_20230504_trace_at_choice_ver2_depth_correct('altern_OFC_20240427');
process_HMM_20230504_trace_at_choice_ver2_depth_correct('altern_Hippo_20240427');
process_HMM_20230504_trace_at_choice_ver2_depth_correct('altern_AC_20240427');
process_HMM_20230504_trace_at_choice_ver2_depth_correct('altern_PPC_20240427');
process_HMM_20230504_trace_at_choice_ver2_depth_correct('altern_M1_20240427');
process_HMM_20230504_trace_at_choice_ver2_depth_correct('altern_STR_20240427');

return


function process_HMM_20230504_trace_at_choice_ver2_depth_correct(folders, kaiseki_number)

switch nargin
    case 0
        hoge
    case 1
        kaiseki_number = 3;
    case 2
        disp('OK to analyze')
    otherwise
        hoge
end

[analysis_dir,depth_def] = eval(folders);

all_during_choice = [];
all_p_during_choice = [];
all_norm_spike = [];
all_p_surprise = [];
all_sound_correct_error = [];
all_p_sound_correct_error = [];
all_neuron_choice_category = [];

for j = 1:14
    all_sound_trace(j).matrix = [];
    all_choice_trace(j).matrix = [];
end

for i = 1:length(analysis_dir)
    disp([i,length(analysis_dir)])
    temp_dir = analysis_dir{i};
    
    [sound_spike, choice_spike, choice_neuron, p_choice_neuron, norm_spike, p_surprise, ...
    sound_correct_error, p_sound_correct_error, neuron_choice_category] = ...
        HMM_ephys_20230504_trace_at_choice_depth(temp_dir, kaiseki_number,depth_def);

    for j = 1:14
        all_sound_trace(j).matrix = [all_sound_trace(j).matrix; sound_spike(j).matrix];
        all_choice_trace(j).matrix = [all_choice_trace(j).matrix; choice_spike(j).matrix];
    end
    
    all_during_choice = [all_during_choice; choice_neuron];
    all_p_during_choice = [all_p_during_choice; p_choice_neuron];
    all_norm_spike = [all_norm_spike; norm_spike];
    all_p_surprise = [all_p_surprise; p_surprise];

    all_sound_correct_error = [all_sound_correct_error; sound_correct_error];
    all_p_sound_correct_error = [all_p_sound_correct_error; p_sound_correct_error];

    all_neuron_choice_category = [all_neuron_choice_category; neuron_choice_category];
end
delete(gcp('nocreate'))

disp([length(all_during_choice),length(all_norm_spike),length(all_p_surprise)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot about low, high neurons

right_neuron = find(all_sound_correct_error(:,1) >= 0);
left_neuron = find(all_sound_correct_error(:,1) < 0);
sig_sound = find(all_p_during_choice(:,3) < 0.01); %Sig diff. between left and right

right_sig_neuron = intersect(right_neuron, sig_sound);
left_sig_neuron = intersect(left_neuron, sig_sound);
non_sig_neuron = setdiff(1:length(all_during_choice), sig_sound);
check_neuron = length(right_sig_neuron) + length(left_sig_neuron) + length(non_sig_neuron);
if check_neuron ~= length(all_during_choice)
    hoge
else
    disp([length(left_sig_neuron),length(right_sig_neuron),length(non_sig_neuron)])
end

length_neuron = length(left_sig_neuron)+length(right_sig_neuron);

make_prefer_trace_2(left_sig_neuron, right_sig_neuron, all_choice_trace,length_neuron);



return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [c_repeat_prefer,c_switch_prefer,c_switch_nonprefer,c_repeat_nonprefer, ...
          e_repeat_prefer,e_switch_prefer,e_switch_nonprefer,e_repeat_nonprefer] ...
           = get_trace_prefer(use_neuron, all_choice_trace, prefer_order)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c_repeat_prefer = all_choice_trace(prefer_order(1,1)).matrix(use_neuron,:);
c_switch_prefer = all_choice_trace(prefer_order(1,2)).matrix(use_neuron,:);
c_switch_nonprefer = all_choice_trace(prefer_order(1,3)).matrix(use_neuron,:);
c_repeat_nonprefer = all_choice_trace(prefer_order(1,4)).matrix(use_neuron,:);

e_repeat_prefer = all_choice_trace(prefer_order(2,1)).matrix(use_neuron,:);
e_switch_prefer = all_choice_trace(prefer_order(2,2)).matrix(use_neuron,:);
e_switch_nonprefer = all_choice_trace(prefer_order(2,3)).matrix(use_neuron,:);
e_repeat_nonprefer = all_choice_trace(prefer_order(2,4)).matrix(use_neuron,:);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function make_prefer_trace_2(left_sig_neuron, right_sig_neuron, all_sound_trace,length_neuron)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

left_order =  [7 9 11 13; 8 10 12 14]; %repeat_prefer, switch_prefer, switch_nonprefer, repeat_nonprefer
right_order = [13 11 9 7; 14 12 10 8]; %repeat_prefer, switch_prefer, switch_nonprefer, repeat_nonprefer

[c_repeat_prefer_l,c_switch_prefer_l,c_switch_nonprefer_l,c_repeat_nonprefer_l, ...
 ~,~,~,~] ...
  = get_trace_prefer(left_sig_neuron, all_sound_trace, left_order);
[c_repeat_prefer_r,c_switch_prefer_r,c_switch_nonprefer_r,c_repeat_nonprefer_r, ...
 ~,~,~,~] ...
  = get_trace_prefer(right_sig_neuron, all_sound_trace, right_order);

c_repeat_prefer = [c_repeat_prefer_l; c_repeat_prefer_r];
c_switch_prefer = [c_switch_prefer_l; c_switch_prefer_r];
c_switch_nonprefer = [c_switch_nonprefer_l; c_switch_nonprefer_r];
c_repeat_nonprefer = [c_repeat_nonprefer_l; c_repeat_nonprefer_r];


figure
plot_mean_se_moto(c_repeat_prefer,[1 0 0],0) %
hold on
plot_mean_se_moto(c_switch_prefer,[255 102 0]./255,0) %
hold on
plot_mean_se_moto(c_switch_nonprefer,[68 168 40]./255,0) %
hold on
plot_mean_se_moto(c_repeat_nonprefer,[0 0 1],0) %
hold on
plot([2000 2000], ylim, ':')
set(gcf,'Position',[200 200 325 319])
yl = ylim;
text(4000,yl(2)*0.9, strcat(string(length_neuron),' neurons'))
xticklabels({'-2',[],'0',[],'2',[],'4'})
xlabel('Time from choice (s)')
ylabel('Zscored activity')
set(gca,'fontname','Arial','Box','off')
set(gcf,'Position',[584,652,295,263])




return


