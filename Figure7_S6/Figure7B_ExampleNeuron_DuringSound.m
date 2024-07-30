
%{
---------------------------------------------------------------------------
First_take number of frames in each tif files
Analyzing imaging data simply
At least for the correct rate
----------------------------------------------------------------------------
%}
function Figure7B_ExampleNeuron_DuringSound(root_dir)


neuronID_rep = 168;
neuronID_zig = 174;

pathname_rep = [root_dir '/Zigzag_repeat_Ephys/E_phys/repeat/W03/OFC/2021-09-30_09-59-21_done_W03_OFC_L_2nd'];
pathname_zig = [root_dir '/Zigzag_repeat_Ephys/E_phys/zigzag/W28/OFC/2022-09-28_14-02-36_W28_OFC_L_1st'];


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

binary_tone = nan(length(Tone_cloud),1);
for i = 1:length(Tone_cloud)
    temp_tone = Tone_cloud(i).matrix;
    %Get the data in all sound
    %     temp0 = find(temp_tone <= 8);
    temp1 = find(temp_tone >= 9);
    binary_tone(i,1) = length(temp1) ./ length(temp_tone);
end

temptemp = Correct_side == 1; %Right is correct
if mean(binary_tone(temptemp)) < 0.5 %low right
    disp('flip tones')

    low_correct = intersect(trial_info.low, trial_info.rewarded);
    high_correct = intersect(trial_info.high, trial_info.rewarded);

    post_low_correct = low_correct+1;
    post_high_correct = high_correct+1;

    pLcR = intersect(trial_info.low, post_high_correct);
    pLcL = intersect(trial_info.high,post_high_correct);
    pRcL = intersect(trial_info.high,post_low_correct);
    pRcR = intersect(trial_info.low,post_low_correct);

else % low left
    low_correct = intersect(trial_info.low, trial_info.rewarded);
    high_correct = intersect(trial_info.high, trial_info.rewarded);

    post_low_correct = low_correct+1;
    post_high_correct = high_correct+1;

    pLcR = intersect(trial_info.high, post_low_correct);
    pLcL = intersect(trial_info.low,post_low_correct);
    pRcL = intersect(trial_info.low,post_high_correct);
    pRcR = intersect(trial_info.high,post_high_correct);

end

[~,~,~,use_trial,low,high,~,~,~,~,~,~,~,~] ...
    = Wang_get_basic_task_structure_20221104(behavior_file);

spike_dir = dir('spike_ch*');
if length(spike_dir) ~= 1
    hoge
end
spike_dir = spike_dir.name;

%Get task parameter
left = find(Chosen_side == 0);
right = find(Chosen_side == 1);

low = intersect(low, use_trial);
high = intersect(high, use_trial);
left = intersect(left, use_trial);
right = intersect(right, use_trial);

%Switch based on previous low or high
pre_low = low + 1;
pre_high = high + 1;

%Get the sound frame
pre_frame = 500;
post_frame = 500;
pre_neuron = [500,1500];
post_neuron = [1000,1000];

cd(spike_dir);

temp_file = sprintf('task_spike_stripe20210520_%d',neuronID);
data = load(temp_file); %spike_mark
[~,~,~, pre_spike_trace,~,~] = get_sound_response_Wang(data.spike_mark, data.spike_filter, frame_sound, ...
    left, right, pre_low, pre_high, pre_frame, post_frame, pre_neuron, post_neuron,pLcL,pLcR,pRcL,pRcR);

plot_ave_spike(pre_spike_trace, pLcL,pLcR,pRcL,pRcR)
xticks(0:500:2000)
xticklabels({-0.5,'0',0.5,'1',1.5})
xlim([0 1700])
xlabel('Time from sound onset (s)')
ylabel('Spike (Hz)')
set(gcf,'Position',[584,652,295,263])
title ('OFC')

return
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_ave_spike(spike_trace, pLcL,pLcR,pRcL,pRcR)
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
function [spike_count,p,p_sound,pre_spike_trace,post_spike_trace,IncreaseDecrease] = get_sound_response_Wang(spike_mark, spike_filter, frame_sound, ...
    ~, ~, ~, ~, pre_frame, post_frame, pre_neuron, post_neuron,pLcL,pLcR,pRcL,pRcR)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

spike_count = nan(length(frame_sound),2);
pre_spike_trace = nan(length(frame_sound),pre_neuron(1)+pre_neuron(2));
post_spike_trace = nan(length(frame_sound),post_neuron(1)+post_neuron(2));

for i = 1:length(frame_sound)
    temp_pre  = frame_sound(i)-pre_frame : frame_sound(i)-1;
    temp_post = frame_sound(i) : frame_sound(i)+post_frame-1;
    temp_pre_all = frame_sound(i)-pre_neuron(1) : frame_sound(i)+pre_neuron(2)-1;
    temp_post_all = frame_sound(i)-post_neuron(1) : frame_sound(i)+post_neuron(2)-1;

    temp_pre = spike_mark(temp_pre);
    temp_post = spike_mark(temp_post);
    spike_count(i,:) = [sum(temp_pre), sum(temp_post)];

    pre_spike_trace(i,:) = spike_filter(temp_pre_all);
    post_spike_trace(i,:) = spike_filter(temp_post_all);
end

%For pre activity, compare between pre_low, pre_high

p = ranksum(spike_count(:,1),spike_count(:,2));
p_sound(1) = ranksum(spike_count(pLcR,2),spike_count(pRcR,2));
p_sound(2) = ranksum(spike_count(pRcL,2),spike_count(pLcL,2));



% p = ranksum(spike_count(pre_low,2),spike_count(pre_high,2));
IncreaseDecrease = mean(spike_count(:,2))-mean(spike_count(:,1));

return

function [Choice_trial,tone_evidence,trial_evidence,use_block,...
    low,high,correct,error,flip_tone,...
    binary_tone,right_trial_all,number_trial_all,right_trial,number_trial] ...
 = Wang_get_basic_task_structure_20221104(filename1)

switch nargin
    case 0
        [filename1, pathname1]=uigetfile('*.mat','Block_mat');
        %[filename1, pathname1]=uigetfile('Bpod*.mat','Block_mat');
        filename1 = [pathname1, filename1];
        load(filename1)
    case 1
        load(filename1)
    otherwise
        hoge
end

check_duration = unique(StimDuration);
if length(check_duration) == 1

else
    disp('cannot use this program')
    length(check_duration)
    hoge
end

Choice_trial = find(Outcome == 1 | Outcome == 2);

% %Outcome
% outcome_EW     = 0; %early withdrawal
% outcome_IC     = 1; %incorrect choice
% outcome_reward = 2; %reward was dispensed (either automatically in early training, or after correct choice)
% outcome_NC     = 3; %no choice was made and time elapsed
% outcome_UN     = 4; %undefined or Free water:

temp_evi = unique(EvidenceStrength);
temp_evi_low  = 0.5 - temp_evi/2;
temp_evi_high = 0.5 + temp_evi/2;
temp_evi_all = [temp_evi_low', temp_evi_high'];
tone_evidence = sort(temp_evi_all);

%Put tone evidence in all trials;
trial_evidence = zeros(length(Outcome),1);
low  = find(Correct_side == 0);
high = find(Correct_side == 1);
temp = Chosen_side == Correct_side;
correct = find(temp == 1);
error   = find(temp == 0);
for i = 1:length(temp_evi)
    temp = find(EvidenceStrength == temp_evi(i));
    temp_left  = intersect(temp,low);
    temp_right = intersect(temp,high);
    trial_evidence(temp_left)  = temp_evi_low(i);
    trial_evidence(temp_right) = temp_evi_high(i);
end

%Make the true tone cloud value
for i = 1:length(Tone_cloud)
    temp_tone = Tone_cloud(i).matrix;
    %Get the data in all sound
    temp1 = find(temp_tone >= 9);
    binary_tone(i,1) = length(temp1) ./ length(temp_tone);
end


%TrialBlock
use_block = TrialCount(InitBlock+1:end);
use_block = intersect(use_block,Choice_trial);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Analyzing more data

%Based on the correct trial, flip the tone cloud
temp = Correct_side == 1; %Right is correct
if mean(binary_tone(temp)) < 0.5 %low for right correct
    disp('flip tones')
    binary_tone = 1 - binary_tone;
    flip_tone = 1;
else
    flip_tone = 0;
end

%Based on the binary tone decide the pseudo tone evidence
clear temp_tone
temp_tone(1).matrix = find(binary_tone == 0);
temp_tone(2).matrix = find(binary_tone > 0 & binary_tone <= 0.35);
temp_tone(3).matrix = find(binary_tone > 0.35 & binary_tone < 0.5);
temp_tone(4).matrix = find(binary_tone > 0.5 & binary_tone < 0.65);
temp_tone(5).matrix = find(binary_tone >= 0.65 & binary_tone < 1);
temp_tone(6).matrix = find(binary_tone == 1);

new_tone_evi = nan(length(binary_tone),1);
for i = 1:6
    new_tone_evi(temp_tone(i).matrix) = tone_evidence(i);
end
%Check nan
temp = unique(isnan(new_tone_evi));
temp = find(temp == 1, 1);
if ~isempty(temp)
    disp([binary_tone, new_tone_evi])
    disp('nan detected')
    hoge
end
trial_evidence = new_tone_evi;

%Switch based on previous low or high
pre_low = low + 1;
pre_high = high + 1;
pre_low = intersect(pre_low, use_block);
pre_high = intersect(pre_high, use_block);

%Get the number of right choice trials in each session
[right_trial.low, number_trial.low] = get_right_choice_trials(pre_low,Chosen_side,trial_evidence,tone_evidence);
[right_trial.high, number_trial.high] = get_right_choice_trials(pre_high,Chosen_side,trial_evidence,tone_evidence);
    
right_trial_all = right_trial.low + right_trial.high;
number_trial_all = number_trial.low + number_trial.high;

%right_trial_all ./ number_trial_all %Right choice prob

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [right_trial, number_trial] = get_right_choice_trials(use_trials,Chosen_side,trial_evidence,tone_evidence)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

use_choice = Chosen_side(use_trials);
use_tone = trial_evidence(use_trials);

for i = 1:length(tone_evidence)
    temp_trial = find(use_tone == tone_evidence(i));
    temp_choice = use_choice(temp_trial);
    number_trial(i) = length(temp_trial);
    right_trial(i) = sum(temp_choice);
end

return

