
%{
----------------------------------------------------------------------------
First_take number of frames in each tif files
Analyzing imaging data simply
At least for the correct rate
----------------------------------------------------------------------------
%}
function Figure8B_ChoiceTrace_20230919(root_dir)


neuronID_rep = 209;
neuronID_zig = 190;

pathname_rep = [root_dir '/Zigzag_repeat_Ephys/E_phys/repeat/W03/OFC/2021-09-30_09-59-21_done_W03_OFC_L_2nd'];
pathname_zig = [root_dir '/Zigzag_repeat_Ephys/E_phys/zigzag/W14/OFC/2022-03-16_13-50-15_W14_OFC_R_1st'];


figure
getTraceData(pathname_rep,neuronID_rep)

figure
getTraceData(pathname_zig,neuronID_zig)

return


function getTraceData(pathname,neuronID)
cd(pathname)
temp_frame_file = dir('task_frame_*');
frame_file = temp_frame_file.name;
load(frame_file);

cd(pathname)
temp_behavior_file = dir('Bpod_mat_*');
behavior_file = temp_behavior_file.name;
load(behavior_file);

trial_info = get_trial_basic_information(Outcome,Correct_side,Chosen_side,...
    EvidenceStrength,Intensity,StimDuration,InitBlock);
use_trial = trial_info.RemoveFirst;
use_trial = intersect(use_trial, trial_info.rewarded);


[~,~,~,use_block,...
    ~,~,~,~,~,...
    ~,~,~,~,~] ...
 = Wang_get_basic_task_structure_20221104(behavior_file);

spike_dir = dir('spike_ch*');
if length(spike_dir) ~= 1
    hoge
end
spike_dir = spike_dir.name;

left = find(Chosen_side == 0);
right = find(Chosen_side == 1);

left_correct = intersect(left, use_trial);
right_correct = intersect(right, use_trial);

%Switch based on previous low or high
post_left = left_correct + 1;
post_right = right_correct + 1;
post_left = intersect(use_block,post_left);
post_right = intersect(use_block,post_right);

pLcL = intersect(post_left, left); 
pLcR  = intersect(post_left, right); 
pRcL  = intersect(post_right, left); 
pRcR  = intersect(post_right, right); 


%Get the  frame
pre_frame = 500;
post_frame = 500;
pre_neuron = [500,2000];
post_neuron = [1000,1000];

cd(spike_dir);

temp_file = sprintf('task_spike_stripe20210520_%d',neuronID);

clear data
data = load(temp_file); %spike_mark
[~,~,pre_spike_trace,~] = get_sound_response_Wang(data.spike_mark, data.spike_filter, frame_choice, frame_spout,...
    left, right, post_left, post_right, pre_frame, post_frame, pre_neuron, post_neuron);

plot_ave_spike(pre_spike_trace, pLcL,pRcR,pLcR,pRcL)
xticks(0:500:2500)
xticklabels({[],0,[],'1',[],'2',[]})
xlabel('Time from choice (s)')
ylabel('Spike (Hz)')
set(gcf,'Position',[584,652,295,263])
title ('OFC')

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_ave_spike(spike_trace, pLcL, pRcR, pLcR,pRcL)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a1 = plot(mean(spike_trace(pLcL,:)*1000),'color',[1 0 1],'LineWidth',1.5); %Blue
hold on
a2 =plot(mean(spike_trace(pLcR,:)*1000),'color',[0 0 1],'LineWidth',1.5); %Green
hold on
a3 =plot(mean(spike_trace(pRcL,:)*1000),'color',[255 0 0]./255,'LineWidth',1.5); %Red
hold on
a4 =plot(mean(spike_trace(pRcR,:)*1000),'color',[0.3010 0.7450 0.9330],'LineWidth',1.5); %Orange

lgd = legend([a1 a2 a3 a4],{'Left--Left', 'Left--Right','Right--Left','Right--Right'},'Location','northwest');
legend('boxoff')
title(lgd,{'Previous--Current'})


return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [spike_count,p,pre_spike_trace,post_spike_trace] = get_sound_response_Wang(spike_mark, spike_filter, frame_choice, frame_spout,...
    left, right, post_left, post_high, pre_frame, post_frame, pre_neuron, post_neuron)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

choice_time = nan(length(frame_choice),1);
for j = 1 : length(frame_choice)
    if isnan(min(frame_choice(j,:)))
        choice_time(j,1) = frame_spout(j,4);
    else
        choice_time(j,1) = min(frame_choice(j,:));
    end
end 

spike_count = nan(length(choice_time),2);
pre_spike_trace = nan(length(choice_time),pre_neuron(1)+pre_neuron(2));
post_spike_trace = nan(length(choice_time),post_neuron(1)+post_neuron(2));

for i = 1:length(choice_time)
    temp_pre  = choice_time(i)-pre_frame : choice_time(i)-1;
    temp_post = choice_time(i) : choice_time(i)+post_frame-1;
    temp_pre_all = choice_time(i)-pre_neuron(1) : choice_time(i)+pre_neuron(2)-1;
    temp_post_all = choice_time(i)-post_neuron(1) : choice_time(i)+post_neuron(2)-1;
    
    temp_pre = spike_mark(temp_pre);
    temp_post = spike_mark(temp_post);
    spike_count(i,:) = [sum(temp_pre), sum(temp_post)];
    
    pre_spike_trace(i,:) = spike_filter(temp_pre_all);
    post_spike_trace(i,:) = spike_filter(temp_post_all);
end

%For pre activity, compare between pre_low, pre_high
p(1) = ranksum(spike_count(post_left,1),spike_count(post_high,2));
p(2) = ranksum(spike_count(left,2),spike_count(right,2));

return
