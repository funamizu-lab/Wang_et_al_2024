
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
function Figure7F_G_HMM_20231129_all_region_at_sound_depth

%repeat
[repeat_OFC_prefer, repeat_OFC_nonprefer, repeat_p(1,:), repeat_neuron(1)] = HMM_20230626_surprise_at_sound_depth('repeat_OFC_20230427',2);
[repeat_AC_prefer, repeat_AC_nonprefer, repeat_p(2,:), repeat_neuron(2)] = HMM_20230626_surprise_at_sound_depth('repeat_AC_20230427',2);
[repeat_PPC_prefer, repeat_PPC_nonprefer, repeat_p(3,:), repeat_neuron(3)] = HMM_20230626_surprise_at_sound_depth('repeat_PPC_20230427',2);

%zigzag

[zigzag_OFC_prefer, zigzag_OFC_nonprefer, zigzag_p(1,:), zigzag_neuron(1)] = HMM_20230626_surprise_at_sound_depth('zigzag_OFC_20230427',2);
[zigzag_Hippo_prefer, zigzag_Hippo_nonprefer, zigzag_p(2,:), zigzag_neuron(2)] = HMM_20230626_surprise_at_sound_depth('zigzag_Hippo_20230427',2);
[zigzag_AC_prefer, zigzag_AC_nonprefer, zigzag_p(3,:), zigzag_neuron(3)] = HMM_20230626_surprise_at_sound_depth('zigzag_AC_20230427',2);
[zigzag_PPC_prefer, zigzag_PPC_nonprefer, zigzag_p(4,:), zigzag_neuron(4)] = HMM_20230626_surprise_at_sound_depth('zigzag_PPC_20230427',2);
[zigzag_M1_prefer, zigzag_M1_nonprefer, zigzag_p(5,:), zigzag_neuron(5)] = HMM_20230626_surprise_at_sound_depth('zigzag_M1_20230427',2);
[zigzag_STR_prefer, zigzag_STR_nonprefer, zigzag_p(6,:), zigzag_neuron(6)] = HMM_20230626_surprise_at_sound_depth('zigzag_STR_20230427',2);


group_repeat_OFC = ones(length(repeat_OFC_prefer),1) * 1;
group_repeat_AC = ones(length(repeat_AC_prefer),1) * 2;
group_repeat_PPC = ones(length(repeat_PPC_prefer),1) * 3;

group_zigzag_OFC = ones(length(zigzag_OFC_prefer),1) * 4;
group_zigzag_Hippo = ones(length(zigzag_Hippo_prefer),1) * 5;
group_zigzag_AC = ones(length(zigzag_AC_prefer),1) * 6;
group_zigzag_PPC = ones(length(zigzag_PPC_prefer),1) * 7;
group_zigzag_M1 = ones(length(zigzag_M1_prefer),1) * 8;
group_zigzag_STR = ones(length(zigzag_STR_prefer),1) * 9;

repeat_prefer = [repeat_OFC_prefer; repeat_AC_prefer; repeat_PPC_prefer];
repeat_nonprefer = [repeat_OFC_nonprefer; repeat_AC_nonprefer; repeat_PPC_nonprefer];
group_repeat = [group_repeat_OFC; group_repeat_AC; group_repeat_PPC];
group_repeat2 = group_repeat + 3;
zigzag_prefer = [zigzag_OFC_prefer; zigzag_Hippo_prefer; zigzag_AC_prefer; zigzag_PPC_prefer; zigzag_M1_prefer; zigzag_STR_prefer];
zigzag_nonprefer = [zigzag_OFC_nonprefer; zigzag_Hippo_nonprefer; zigzag_AC_nonprefer; zigzag_PPC_nonprefer; zigzag_M1_nonprefer; zigzag_STR_nonprefer];
group_zigzag = [group_zigzag_OFC; group_zigzag_Hippo; group_zigzag_AC; group_zigzag_PPC; group_zigzag_M1; group_zigzag_STR];
group_zigzag2 = group_zigzag + 6;

figure(1)
boxplot(cat(1, repeat_prefer, repeat_nonprefer), cat(1, group_repeat, group_repeat2),'symbol', '',...
    'Labels',{'OFC','AC','PPC','OFC','AC','PPC'},'Colors', [0 0.5 0])
yline(0,':k')
xt = 1:6;
xt = xt - 0.3;
yt = [2 2 2 2 2 2];
str = {num2str(repeat_p(1,1)),num2str(repeat_p(2,1)),num2str(repeat_p(3,1)),num2str(repeat_p(1,2)),num2str(repeat_p(2,2)),num2str(repeat_p(3,2))};
text(xt,yt,str)
text(1, 2.3, 'Preferred Sound')
text(3, 2.3, 'Non-preferred Sound')
ylim([-2.5 2.5])
ylabel('Repeat - Alterante')
title('During sound (Repeating condition)')
set(gca,'FontName', 'Arial')
set(gcf,'Position',[441,351,616,360])



figure(2)
boxplot(cat(1, zigzag_prefer, zigzag_nonprefer), cat(1, group_zigzag, group_zigzag2),'symbol', '',...
      'Labels',{'OFC','HPC','AC','PPC','M1','STR','OFC','HPC','AC','PPC','M1','STR'},'Colors', 'm')
xt = 1:12;
xt = xt - 0.3;
yt = [2 2 2 2 2 2 2 2 2 2 2 2];
str = {num2str(zigzag_p(1,1)),num2str(zigzag_p(2,1)),num2str(zigzag_p(3,1)),num2str(zigzag_p(4,1)),num2str(zigzag_p(5,1)),num2str(zigzag_p(6,1))...
    num2str(zigzag_p(1,2)),num2str(zigzag_p(2,2)),num2str(zigzag_p(3,2)),num2str(zigzag_p(4,2)),num2str(zigzag_p(5,2)),num2str(zigzag_p(6,2))};
text(xt,yt,str)
text(2.5, 2.3, 'Preferred Sound')
text(7.5, 2.3, 'Non-preferred Sound')
yline(0,':k')
ylim([-2.5 2.5])
ylabel('Repeat - Alterante')
title('During sound (Alternating condition)')
set(gca,'FontName', 'Arial')
set(gcf,'Position',[389,433,1044,360])


return


function [prefer_sabun, nonprefer_sabun, p_sabun, length_neuron] = HMM_20230626_surprise_at_sound_depth(folders, kaiseki_number)


[analysis_dir,depth_def] = eval(folders);

all_pre_sound = [];
all_p_pre_sound = [];
all_norm_spike = [];
all_p_surprise = [];
all_sound_correct_error = [];
all_tone_correct_error = [];
all_p_sound_correct_error = [];
all_p_tone_correct_error = [];
all_neuron_choice_category = [];

for i = 1:length(analysis_dir)
    disp([i,length(analysis_dir)])
    temp_dir = analysis_dir{i};
    
    [sound_neuron, p_sound_neuron, norm_spike, p_surprise, ...
    sound_correct_error, p_sound_correct_error, neuron_choice_category,tone_correct_error,p_tone_correct_error] = ...
        HMM_20240710_surprise_at_sound_depth(temp_dir, kaiseki_number,depth_def);
    
    all_pre_sound = [all_pre_sound; sound_neuron];
    all_p_pre_sound = [all_p_pre_sound; p_sound_neuron];
    all_norm_spike = [all_norm_spike; norm_spike];
    all_p_surprise = [all_p_surprise; p_surprise];

    all_sound_correct_error = [all_sound_correct_error; sound_correct_error];
    all_p_sound_correct_error = [all_p_sound_correct_error; p_sound_correct_error];

    all_tone_correct_error = [all_tone_correct_error;tone_correct_error];
    all_p_tone_correct_error = [all_p_tone_correct_error;p_tone_correct_error];

    all_neuron_choice_category = [all_neuron_choice_category; neuron_choice_category];
end
delete(gcp('nocreate'))



target_analysis = 1;
%target_analysis = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Based on focusing on the sound or choice, should change the analysis way
if target_analysis == 1 %Focus on choice
    %Focus on choice
    target_sig = 1;
    target_choice = [1,2]; %left, right
    target_surprise = [25,26,27,28];
elseif target_analysis == 2 %Focus on sound
    %Focus on sound
    target_sig = 2;
    target_choice = [3,4]; %left, right
    target_surprise = [29,30,31,32];
end

[prefer_sabun, nonprefer_sabun, p_sabun, length_neuron] = ...
    get_analysis_sound_or_choice(target_sig, target_choice, target_surprise, ...
    all_neuron_choice_category,all_pre_sound,all_p_pre_sound,all_sound_correct_error,all_tone_correct_error,all_norm_spike);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [prefer_sabun, nonprefer_sabun, p_sabun, length_neuron] = ...
    get_analysis_sound_or_choice(target_sig, ~, target_surprise, ...
    ~,all_pre_sound,all_p_pre_sound,all_sound_correct_error,~,all_norm_spike)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot about low, high neurons

right_neuron = find(all_sound_correct_error(:,1) > 0);
left_neuron = find(all_sound_correct_error(:,1) < 0);
%sig_sound = find(all_p_pre_sound(:,1) < 0.01); %Sig diff. between left and right
sig_sound = find(all_p_pre_sound(:,target_sig) < 0.01); %Sig diff. between left and right

right_sig_neuron = intersect(right_neuron, sig_sound);
left_sig_neuron = intersect(left_neuron, sig_sound);
non_sig_neuron = setdiff(1:length(all_pre_sound), sig_sound);
check_neuron = length(right_sig_neuron) + length(left_sig_neuron) + length(non_sig_neuron);
if check_neuron ~= length(all_pre_sound)
    hoge
else
    disp([length(left_sig_neuron),length(right_sig_neuron),length(non_sig_neuron)])
end

%Use only on low sig neuron
[Lsig_PostLow_left, Lsig_PostHigh_left, Lsig_PostLow_right, Lsig_PostHigh_right] = ...
    test_surprise_sound_choice(all_norm_spike, left_sig_neuron, target_surprise);
%Use only on high sig neuron
[Rsig_PostLow_left, Rsig_PostHigh_left, Rsig_PostLow_right, Rsig_PostHigh_right] = ...
    test_surprise_sound_choice(all_norm_spike, right_sig_neuron, target_surprise);

%Activity is already flipped
Repeat_prefer = [Lsig_PostLow_left; Rsig_PostHigh_right];
Repeat_nonprefer = [Lsig_PostHigh_right; Rsig_PostLow_left];
Switch_prefer = [Lsig_PostHigh_left; Rsig_PostLow_right];
Switch_nonprefer = [Lsig_PostLow_right; Rsig_PostHigh_left];

%plot_surprise_choice(Repeat_prefer, Switch_prefer, Switch_nonprefer, Repeat_nonprefer)
plot_surprise_choice(Lsig_PostLow_left, Lsig_PostHigh_left, Lsig_PostLow_right, Lsig_PostHigh_right);
plot_surprise_choice(Rsig_PostHigh_right, Rsig_PostLow_right, Rsig_PostHigh_left, Rsig_PostLow_left);

[prefer_sabun, nonprefer_sabun, p_sabun] = plot_surprise_choice(Repeat_prefer, Switch_prefer, Switch_nonprefer, Repeat_nonprefer);
temp1 = length(prefer_sabun);
temp2 = length(nonprefer_sabun);

if temp1 ~= temp2
    hoge
else
    length_neuron = temp1;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [PostLow_left, PostHigh_left, PostLow_right, PostHigh_right] = ...
    test_surprise_sound_choice(all_norm_spike, low_sig_neuron, target_surprise)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PostLow_left = all_norm_spike(low_sig_neuron,target_surprise(1)); %PostLowCorrect_left
PostHigh_left = all_norm_spike(low_sig_neuron,target_surprise(2)); %PostHighCorrect_left
PostLow_right = all_norm_spike(low_sig_neuron,target_surprise(3)); %PostLowCorrect_right
PostHigh_right = all_norm_spike(low_sig_neuron,target_surprise(4)); %PostHighCorrect_right
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function plot_surprise_choice(PostLow_left, PostHigh_left, PostLow_right, PostHigh_right)
function [prefer_sabun, nonprefer_sabun, p_sabun] = plot_surprise_choice(Repeat_prefer, Switch_prefer, Switch_nonprefer, Repeat_nonprefer)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


prefer_sabun = Repeat_prefer-Switch_prefer;
nonprefer_sabun = Repeat_nonprefer - Switch_nonprefer;
if ~isempty(prefer_sabun)
    p_sabun(1) = signrank(prefer_sabun);
    p_sabun(2) = signrank(nonprefer_sabun);

    signrank(Repeat_prefer, Switch_prefer)
    signrank(Repeat_nonprefer, Switch_nonprefer)
else
    p_sabun = nan(1,2);
end

return

function [sound_neuron, p_sound_neuron, norm_spike, p_surprise, ...
    sound_correct_error, p_sound_correct_error, neuron_choice_category,tone_correct_error,p_tone_correct_error] = ...
    HMM_20240710_surprise_at_sound_depth(pathname, kaiseki_number,depth_def)

switch nargin
    case 0
        pathname = pwd;
    case 3
        disp('OK to analyze')
    otherwise
        hoge
end
cd(pathname)

% [filename1, pathname1,findex]=uigetfile('*.*','frame file');
% filename1 = [pathname1,filename1];
% load(filename1)
temp = dir('HMM_spike_count_neurons_20230626*');
if length(temp) ~= 1
    hoge
end
load(temp.name);
%neuron_index p_index

temp = dir('sig_HMM_neurons_20230310*');
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
    
    hoge
end

%Get sig neurons
new_p_thre = 10;
%sig_before_sound = get_sig_neuron_time_window(p_task,new_p_thre,[6:15]); %1000ms before sound
sig_before_sound = get_sig_neuron_time_window(p_task,new_p_thre,10:15); %600ms before sound
sig_during_sound = get_sig_neuron_time_window(p_task,new_p_thre,16:21); %During sound
%sig_during_choice = get_sig_neuron_time_window(p_task2,new_p_thre,[6:8]); %During choice
sig_during_choice = get_sig_neuron_time_window(p_task2,new_p_thre,6:15); %During choice 1000ms

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if kaiseki_number == 1 %Before sound
    use_sig_neuron = sig_before_sound;
    use_neuron_index = neuron_index.before_sound;
    use_p_index = p_index.before_sound;
    use_norm_spike = norm_spike_group.before_sound;
    use_p_surprise = p_surprise.before_sound;

elseif kaiseki_number == 2 %during sound
    use_sig_neuron = sig_during_sound;
    use_neuron_index = neuron_index.all_sound;
    use_p_index = p_index.all_sound;
    use_norm_spike = norm_spike_group.all_sound;
    use_p_surprise = p_surprise.all_sound;
    
elseif kaiseki_number == 3 %during choice %1000ms
    use_sig_neuron = sig_during_choice;
    use_neuron_index = neuron_index.choice2;
    use_p_index = p_index.choice2;
    use_norm_spike = norm_spike_group.choice2;
    use_p_surprise = p_surprise.choice2;
    
else
    hopge
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ADD depth definition
use_sig_neuron = intersect(use_sig_neuron,depth_neuron);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%The activity depended on the current choice or sound:
%How are they modulated?
%And do they have surprised activity?

%First, focus on the choice or sound index
%Neuron index for low or high sounds
sound_neuron(:,1) = use_neuron_index(use_sig_neuron,1); %Left right: choice
sound_neuron(:,2) = use_neuron_index(use_sig_neuron,2); %Low high: sound
p_sound_neuron(:,1) = use_p_index(use_sig_neuron,1); %Left right: choice
p_sound_neuron(:,2) = use_p_index(use_sig_neuron,2); %Low high

sound_correct_error(:,1) = use_neuron_index(use_sig_neuron,4); %Correct Low high
sound_correct_error(:,2) = use_neuron_index(use_sig_neuron,5); %Error Low high
p_sound_correct_error(:,1) = use_p_index(use_sig_neuron,4);
p_sound_correct_error(:,2) = use_p_index(use_sig_neuron,5);

tone_correct_error(:,1) = use_neuron_index(use_sig_neuron,12); %Correct Low high
tone_correct_error(:,2) = use_neuron_index(use_sig_neuron,13); %Error Low high
p_tone_correct_error(:,1) = use_p_index(use_sig_neuron,12);
p_tone_correct_error(:,2) = use_p_index(use_sig_neuron,13);

norm_spike = use_norm_spike(use_sig_neuron,:);
p_surprise = use_p_surprise(use_sig_neuron,:);

temp = norm_spike(:,2) - norm_spike(:,1); %right-left
neuron_choice_category = zeros(length(temp),1);
temp = temp > 0;
neuron_choice_category(temp) = 1; %right=1, left=0

%%%
%Remove NAN
test = [sound_neuron, norm_spike, sound_correct_error];
test = mean(test,2);
test = isnan(test);
test_nan = find(test == 1);
sound_neuron(test_nan,:) = [];
p_sound_neuron(test_nan,:) = [];
norm_spike(test_nan,:) = [];
p_surprise(test_nan,:) = [];
sound_correct_error(test_nan,:) = [];
p_sound_correct_error(test_nan,:) = [];
neuron_choice_category(test_nan,:) = [];

return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sig_neuron = get_sig_neuron_time_window(p_task,new_p_thre,window)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
use_p_sound = min(p_task(:,window),[],2);
use_p_sound = -log10(use_p_sound);
sig_neuron = find(use_p_sound > new_p_thre);

return


