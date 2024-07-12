
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
function Figure7C_S5A_HMM_20231129_index_at_sound_depth

%repeat
process_HMM_20230626_index_surprise_at_sound_depth('repeat_OFC_20230427');
process_HMM_20230626_index_surprise_at_sound_depth('repeat_AC_20230427');
process_HMM_20230626_index_surprise_at_sound_depth('repeat_PPC_20230427');

%zigzag
process_HMM_20230626_index_surprise_at_sound_depth('zigzag_OFC_20230427');
process_HMM_20230626_index_surprise_at_sound_depth('zigzag_Hippo_20230427');
process_HMM_20230626_index_surprise_at_sound_depth('zigzag_AC_20230427');
process_HMM_20230626_index_surprise_at_sound_depth('zigzag_PPC_20230427');
process_HMM_20230626_index_surprise_at_sound_depth('zigzag_M1_20230427');
process_HMM_20230626_index_surprise_at_sound_depth('zigzag_STR_20230427');


return

function process_HMM_20230626_index_surprise_at_sound_depth(folders, kaiseki_number)

switch nargin
    case 0
        hoge
    case 1
        kaiseki_number = 2;
    case 2
        disp('OK to analyze')
    otherwise
        hoge
end

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
        HMM_ephys_20230710_surprise_at_sound_depth(temp_dir, kaiseki_number,depth_def);

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
get_analysis_sound_or_choice(target_sig, target_choice, target_surprise, ...
    all_neuron_choice_category,all_pre_sound,all_p_pre_sound,all_sound_correct_error,all_tone_correct_error,all_norm_spike);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function get_analysis_sound_or_choice(target_sig, ~, ~, ...
    ~,all_pre_sound,all_p_pre_sound,all_sound_correct_error,all_tone_correct_error,~)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

right_neuron = find(all_sound_correct_error(:,1) >= 0);
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

length_neuron = size(all_tone_correct_error,1);
[~, ~,r,~] = plot_index_abs_index_sig(all_sound_correct_error(:,1), all_sound_correct_error(:,2));



figure
plot_index_abs_index_sig_temp(all_sound_correct_error(:,1), all_sound_correct_error(:,2), right_sig_neuron,left_sig_neuron,non_sig_neuron);
text(0.2,-0.8,['r = ' string(r)])
text(0.2,-0.6,strcat(string(length_neuron), ' neurons'))
set(gca,'xlim',[-1, 1],'ylim',[-1, 1],'fontname','Arial','Box','off')
set(gcf,'Position',[584,652,295,263])



return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p, median_value,rho,pval] = plot_index_abs_index_sig(index1, index2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[rho,pval] = corr(index1,index2);
p = signrank(abs(index1), abs(index2));
median_value = [median(abs(index1)), median(abs(index2))];

return


function [rho,pval] = plot_index_abs_index_sig_temp(index1, index2, right_neuron,left_neuron,non_sig_neuron)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

[rho,pval] = corr(index1,index2);

function [sound_neuron, p_sound_neuron, norm_spike, p_surprise, ...
    sound_correct_error, p_sound_correct_error, neuron_choice_category,tone_correct_error,p_tone_correct_error] = ...
    HMM_ephys_20230710_surprise_at_sound_depth(pathname, kaiseki_number,depth_def)

switch nargin
    case 0
        pathname = pwd;
    case 3
        disp('OK to analyze')
    otherwise
        hoge
end
cd(pathname)

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
sig_before_sound = get_sig_neuron_time_window(p_task,new_p_thre,10:15); %600ms before sound
sig_during_sound = get_sig_neuron_time_window(p_task,new_p_thre,16:21); %During sound
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
