
%{
----------------------------------------------------------------------------
Determine the time window for task relevant neurons
%Start of task
%Sound on
%Sound off
%Before choice
%After choice (0sec)
%After choice (1sec)
%After choice (2sec)
%p = 0.001
%Each epoch, predict the prior, sensory and choice (integration)
----------------------------------------------------------------------------
%}
function Figure8C_S7A_HMM_20230504_index_choice_ver2_depth_correct


%repeat
 
process_HMM_20230504_index_at_choice_ver2_depth_correct('repeat_OFC_20230427');
process_HMM_20230504_index_at_choice_ver2_depth_correct('repeat_Hippo_20230427');
process_HMM_20230504_index_at_choice_ver2_depth_correct('repeat_AC_20230427');
process_HMM_20230504_index_at_choice_ver2_depth_correct('repeat_PPC_20230427');

%zigzag

process_HMM_20230504_index_at_choice_ver2_depth_correct('zigzag_OFC_20230427');
process_HMM_20230504_index_at_choice_ver2_depth_correct('zigzag_Hippo_20230427');
process_HMM_20230504_index_at_choice_ver2_depth_correct('zigzag_AC_20230427');
process_HMM_20230504_index_at_choice_ver2_depth_correct('zigzag_PPC_20230427');
process_HMM_20230504_index_at_choice_ver2_depth_correct('zigzag_M1_20230427');
process_HMM_20230504_index_at_choice_ver2_depth_correct('zigzag_STR_20230427');



function process_HMM_20230504_index_at_choice_ver2_depth_correct(folders, kaiseki_number)

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

for i = 1:length(analysis_dir)
    disp([i,length(analysis_dir)])
    temp_dir = analysis_dir{i};
    
    [choice_neuron, p_choice_neuron, norm_spike, p_surprise, ...
    sound_correct_error, p_sound_correct_error, neuron_choice_category] = ...
        HMM_ephys_20230504_surprise_at_choice_depth(temp_dir, kaiseki_number,depth_def);

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

size(all_sound_correct_error)
[~, ~,r,~] = get_corr(all_sound_correct_error(:,1), all_sound_correct_error(:,2));

length_neuron = size(all_sound_correct_error,1);

figure
plot_index_abs_index_sig_temp(all_sound_correct_error(:,1), all_sound_correct_error(:,2), right_sig_neuron,left_sig_neuron,non_sig_neuron);
text(0.2,-0.8,['r = ' string(r)])
text(0.2,-0.6,strcat(string(length_neuron), ' neurons'))
set(gca,'xlim',[-1, 1],'ylim',[-1, 1],'fontname','Arial','Box','off')
set(gcf,'Position',[584,652,295,263])

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p, median_value,rho,pval] = get_corr(index1, index2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[rho,pval] = corr(index1,index2);
p = signrank(abs(index1), abs(index2));
median_value = [median(abs(index1)), median(abs(index2))];

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_index_abs_index_sig_temp(index1, index2, right_neuron,left_neuron,non_sig_neuron)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ax2 = nexttile;
plot([-1 1],[0 0],'k:')
hold on
plot([0 0],[-1 1],'k:')
hold on
plot([-1 1],[-1 1],'k:')
hold on
%plot(index1, index2, 'b.')
plot(index1(non_sig_neuron), index2(non_sig_neuron), 'k.')
hold on
plot(index1(left_neuron), index2(left_neuron), 'b.')
hold on
plot(index1(right_neuron), index2(right_neuron), 'r.')
box(ax2,'off')
set(gcf,'Position',[200 200 260 256])
set(gca,'linewidth',1,'fontsize',10,'fontname','Arial');


