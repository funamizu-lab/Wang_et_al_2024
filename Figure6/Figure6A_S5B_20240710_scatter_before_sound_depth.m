
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
function Figure6A_S5B_20240710_scatter_before_sound_depth


% % repeat
process_HMM_20240710_surprise_before_sound_depth('repeat_OFC_20230427');
process_HMM_20240710_surprise_before_sound_depth('repeat_AC_20230427');

% % zigzag
process_HMM_20240710_surprise_before_sound_depth('zigzag_OFC_20230427');
process_HMM_20240710_surprise_before_sound_depth('zigzag_AC_20230427');
process_HMM_20240710_surprise_before_sound_depth('zigzag_PPC_20230427');
process_HMM_20240710_surprise_before_sound_depth('zigzag_M1_20230427');
process_HMM_20240710_surprise_before_sound_depth('zigzag_STR_20230427');


return


function process_HMM_20240710_surprise_before_sound_depth(folders, kaiseki_number)

switch nargin
    case 0
        hoge
    case 1
        kaiseki_number = 1;
    case 2
        disp('OK to analyze')
    otherwise
        hoge
end

brainarea = folders(1:10);


[analysis_dir,depth_def] = eval(folders);

all_pre_sound = [];
all_p_pre_sound = [];
all_norm_spike = [];
all_p_surprise = [];
all_choice_correct_error = [];
all_p_choice_correct_error = [];
all_neuron_choice_category = [];

all_tone_correct_error = [];
all_p_tone_correct_error = [];

all_pre_tone_correct_error = [];
all_p_pre_tone_correct_error = [];

all_pre_choice_correct_error = [];
all_p_pre_choice_correct_error = [];


for i = 1:length(analysis_dir)
    disp([i,length(analysis_dir)])
    temp_dir = analysis_dir{i};
    
    [pre_sound, p_pre_sound, norm_spike, p_surprise, ...
    choice_correct_error, p_choice_correct_error, neuron_choice_category,...
    tone_correct_error, p_tone_correct_error, pre_tone_correct_error, p_pre_tone_correct_error,...
    pre_choice_correct_error, p_pre_choice_correct_error, ~,~,~,~,~,~,~,~,~,~,~,~,~,~] = ...
        HMM_ephys_20240709_surprise_before_sound_depth(temp_dir, kaiseki_number,depth_def);

%     p_pre_sound(:,1) = use_p_index(use_sig_neuron,1); %current choice index
%     p_pre_sound(:,2) = use_p_index(use_sig_neuron,6); %pre tone index
%     p_pre_sound(:,3) = use_p_index(use_sig_neuron,7); %pre choice index
%     p_pre_sound(:,4) = use_p_index(use_sig_neuron,10); %pre SoundCorrect, current left
%     p_pre_sound(:,5) = use_p_index(use_sig_neuron,11); %pre SoundCorrect, current right
%     p_pre_sound(:,6) = use_p_index(use_sig_neuron,11); % current tone index

    all_pre_sound = [all_pre_sound; pre_sound];
    all_p_pre_sound = [all_p_pre_sound; p_pre_sound];
    all_norm_spike = [all_norm_spike; norm_spike];
    all_p_surprise = [all_p_surprise; p_surprise];

    all_choice_correct_error = [all_choice_correct_error; choice_correct_error];
    all_p_choice_correct_error = [all_p_choice_correct_error; p_choice_correct_error];

    all_tone_correct_error = [all_tone_correct_error; tone_correct_error]; % current tone correct index
    all_p_tone_correct_error = [all_p_tone_correct_error; p_tone_correct_error]; % current tone error index
       
    all_pre_tone_correct_error = [all_pre_tone_correct_error; pre_tone_correct_error]; % pre tone correct index
    all_p_pre_tone_correct_error = [all_p_pre_tone_correct_error; p_pre_tone_correct_error]; % pre tone error index
        
    all_pre_choice_correct_error = [all_pre_choice_correct_error; pre_choice_correct_error]; % pre choice correct index
    all_p_pre_choice_correct_error = [all_p_pre_choice_correct_error; p_pre_choice_correct_error]; % pre choice error index

    all_neuron_choice_category = [all_neuron_choice_category; neuron_choice_category];


end
delete(gcp('nocreate'))

disp([length(all_pre_sound),length(all_norm_spike),length(all_p_surprise)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if contains(brainarea,'repeat')
    right_neuron = find(all_pre_choice_correct_error(:,1) >= 0);
    left_neuron = find(all_pre_choice_correct_error(:,1) < 0);
    sig_sound = find(all_p_pre_sound(:,3) < 0.01); %Sig diff. between left and right
else 
    left_neuron = find(all_pre_choice_correct_error(:,1) >= 0);
    right_neuron = find(all_pre_choice_correct_error(:,1) < 0);
    sig_sound = find(all_p_pre_sound(:,3) < 0.01); %Sig diff. between left and right
end

right_sig_neuron = intersect(right_neuron, sig_sound);
left_sig_neuron = intersect(left_neuron, sig_sound);
non_sig_neuron = setdiff(1:length(all_p_pre_sound), sig_sound);
check_neuron = length(right_sig_neuron) + length(left_sig_neuron) + length(non_sig_neuron);
if check_neuron ~= length(all_p_pre_sound)
    hoge
else
    disp([length(left_sig_neuron),length(right_sig_neuron),length(non_sig_neuron)])
end

length_neuron = length(left_sig_neuron)+length(right_sig_neuron);
%Use only on low sig neuron
[L_sig_PostLow_left, L_sig_PostHigh_left, L_sig_PostLow_right, L_sig_PostHigh_right] = ...
    test_surprise_choice(all_norm_spike, left_sig_neuron, [0 0 1]);
%Use only on high sig neuron
[R_sig_PostLow_left, R_sig_PostHigh_left, R_sig_PostLow_right, R_sig_PostHigh_right] = ...
    test_surprise_choice(all_norm_spike, right_sig_neuron, [1 0 0]);

%Activity is already flipped

if contains(brainarea,'repeat')
    Repeat_prefer = [L_sig_PostLow_left; R_sig_PostHigh_right];
    Repeat_nonprefer = [L_sig_PostHigh_right; R_sig_PostLow_left];
    Switch_prefer = [L_sig_PostLow_right; R_sig_PostHigh_left];
    Switch_nonprefer = [L_sig_PostHigh_left; R_sig_PostLow_right];
else
    
    Repeat_prefer = [L_sig_PostHigh_right; R_sig_PostLow_left];
    Repeat_nonprefer = [L_sig_PostLow_left; R_sig_PostHigh_right];
    Switch_prefer = [L_sig_PostHigh_left; R_sig_PostLow_right];
    Switch_nonprefer = [L_sig_PostLow_right; R_sig_PostHigh_left];  
end

length_L_sig = length(left_sig_neuron);

prefer_sabun = Repeat_prefer-Switch_prefer;
nonprefer_sabun = Repeat_nonprefer-Switch_nonprefer;

p_sabun(1) = signrank(prefer_sabun);
p_sabun(2) = signrank(nonprefer_sabun);

figure
plot_surprise_choice_prefer(Repeat_prefer(1:length_L_sig),Switch_prefer(1:length_L_sig) ,[0 0 1]);
hold on
plot_surprise_choice_prefer(Repeat_prefer(length_L_sig+1:end),Switch_prefer(length_L_sig+1:end), [1 0 0]);
text(6,0,['p = ' string(p_sabun(1))])
text(6,2,strcat(string(length_neuron),' neurons'))
set(gca,'fontname','Arial','Box','off')
set(gcf,'Position',[584,652,295,263])
title({'Preferred choice';brainarea(8:10)})


figure
plot_surprise_choice_nonprefer(Repeat_nonprefer(1:length_L_sig),Switch_nonprefer(1:length_L_sig) ,[0 0 1]);
hold on
plot_surprise_choice_nonprefer(Repeat_nonprefer(length_L_sig+1:end),Switch_nonprefer(length_L_sig+1:end) ,[1 0 0]);
text(6,0,['p = ' string(p_sabun(2))])
text(6,2,strcat(string(length_neuron),' neurons'))
set(gca,'fontname','Arial','Box','off')
set(gcf,'Position',[584,652,295,263])
title({'Non-Preferred choice';brainarea(8:10)})

return






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [PostLow_left, PostHigh_left, PostLow_right, PostHigh_right] = ...
    test_surprise_choice(all_norm_spike, low_sig_neuron, ~)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PostLow_left = all_norm_spike(low_sig_neuron,25); %PostLowCorrect_left
PostHigh_left = all_norm_spike(low_sig_neuron,26); %PostHighCorrect_left
PostLow_right = all_norm_spike(low_sig_neuron,27); %PostLowCorrect_right
PostHigh_right = all_norm_spike(low_sig_neuron,28); %PostHighCorrect_right


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_surprise_choice_prefer(Repeat_prefer, Switch_prefer, plot_color)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% figure
% plot(predict_ave, choice_ave, 'k.')
plot(Repeat_prefer, Switch_prefer, '.','color', plot_color) %same choice, X_repeat or Y_zigzag
hold on
plot([-1, 10],[-1, 10],'k')
set(gca,'xlim',[-1, 10],'ylim',[-1, 10])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_surprise_choice_nonprefer(Switch_nonprefer, Repeat_nonprefer, plot_color)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% figure
% plot(predict_ave, choice_ave, 'k.')

plot(Switch_nonprefer, Repeat_nonprefer, '.','color', plot_color)  %same choice, X_zigzag or Y_repeat
hold on
plot([-1, 10],[-1, 10],'k')
set(gca,'xlim',[-1, 10],'ylim',[-1, 10])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pre_sound, p_pre_sound, norm_spike, p_surprise, ...
    choice_correct_error, p_choice_correct_error, neuron_choice_category,...
    tone_correct_error, p_tone_correct_error, pre_tone_correct_error, p_pre_tone_correct_error,...
    pre_choice_correct_error, p_pre_choice_correct_error,...
    SameDif_preTone_correct_error,p_SameDif_preTone_correct_error,...
    preCorrectTone_index_correct_correct_error,p_preCorrectTone_index_correct_correct_error,...
    preCorrectChoice_index_correct_correct_error,p_preCorrectChoice_index_correct_correct_error,...
    pre_same_pLcL_pLcH_index,p_pre_same_pLcL_pLcH_index,...
    pre_same_pHcL_pHcH_index,p_pre_same_pHcL_pHcH_index,...
    cur_same_pLcL_pHcL_index,p_cur_same_pLcL_pHcL_index,...
    cur_same_pLcH_pHcH_index,p_cur_same_pLcH_pHcH_index] = ...
    HMM_ephys_20240709_surprise_before_sound_depth(pathname, kaiseki_number,depth_def)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


switch nargin
    case 0
        pathname = pwd;
    case 3
        disp('OK to analyze')
    otherwise
        hoge
end
cd(pathname)

temp = dir('HMM_spike_count_neurons_20240120*');
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


temp = dir('depth_spike_20230427*');
if length(temp) == 1
    load(temp.name);
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

% %Get sig neurons
new_p_thre = 10;

sig_before_sound = get_sig_neuron_time_window(p_task,new_p_thre,10:15); %600ms before sound
sig_during_sound = get_sig_neuron_time_window(p_task,new_p_thre,16:21); %During sound
sig_during_choice = get_sig_neuron_time_window(p_task2,new_p_thre,6:15); %During choice 1000ms

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%neuron_index.before_sound
%neuron_index.init_sound
%neuron_index.late_sound
%neuron_index.all_sound
%neuron_index.before_choice
%neuron_index.choice1
%neuron_index.choice2
%neuron_index.choice3
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

pre_sound(:,1) = use_neuron_index(use_sig_neuron,1); %current_choice
pre_sound(:,2) = use_neuron_index(use_sig_neuron,6); %pre sound
pre_sound(:,3) = use_neuron_index(use_sig_neuron,7); %pre choice
pre_sound(:,4) = use_neuron_index(use_sig_neuron,10); %pre SoundCorrect, current left
pre_sound(:,5) = use_neuron_index(use_sig_neuron,11); %pre SoundCorrect, current right

pre_sound(:,6) = use_neuron_index(use_sig_neuron,2); %current tone index


p_pre_sound(:,1) = use_p_index(use_sig_neuron,1); %current choice index
p_pre_sound(:,2) = use_p_index(use_sig_neuron,6); %pre tone index
p_pre_sound(:,3) = use_p_index(use_sig_neuron,7); %pre choice index
p_pre_sound(:,4) = use_p_index(use_sig_neuron,10); %pre SoundCorrect, current left
p_pre_sound(:,5) = use_p_index(use_sig_neuron,11); %pre SoundCorrect, current right

p_pre_sound(:,6) = use_p_index(use_sig_neuron,2); %current tone index

choice_correct_error(:,1) = use_neuron_index(use_sig_neuron,4); %Correct choice index
choice_correct_error(:,2) = use_neuron_index(use_sig_neuron,5); %Error choice index
p_choice_correct_error(:,1) = use_p_index(use_sig_neuron,4);
p_choice_correct_error(:,2) = use_p_index(use_sig_neuron,5);

tone_correct_error(:,1) = use_neuron_index(use_sig_neuron,12); %Correct choice index
tone_correct_error(:,2) = use_neuron_index(use_sig_neuron,13); %Error choice index
p_tone_correct_error(:,1) = use_p_index(use_sig_neuron,12);
p_tone_correct_error(:,2) = use_p_index(use_sig_neuron,13);

pre_tone_correct_error(:,1) = use_neuron_index(use_sig_neuron,14); %Correct choice index
pre_tone_correct_error(:,2) = use_neuron_index(use_sig_neuron,15); %Error choice index
p_pre_tone_correct_error(:,1) = use_p_index(use_sig_neuron,14);
p_pre_tone_correct_error(:,2) = use_p_index(use_sig_neuron,15);

pre_choice_correct_error(:,1) = use_neuron_index(use_sig_neuron,16); %Correct choice index
pre_choice_correct_error(:,2) = use_neuron_index(use_sig_neuron,17); %Error choice index
p_pre_choice_correct_error(:,1) = use_p_index(use_sig_neuron,16);
p_pre_choice_correct_error(:,2) = use_p_index(use_sig_neuron,17);

SameDif_preTone_correct_error(:,1) = use_neuron_index(use_sig_neuron,18);
SameDif_preTone_correct_error(:,2) = use_neuron_index(use_sig_neuron,19);
p_SameDif_preTone_correct_error(:,1) = use_p_index(use_sig_neuron,18);
p_SameDif_preTone_correct_error(:,2) = use_p_index(use_sig_neuron,19);


preCorrectTone_index_correct_correct_error(:,1) = use_neuron_index(use_sig_neuron,20);
preCorrectTone_index_correct_correct_error(:,2) = use_neuron_index(use_sig_neuron,21);
p_preCorrectTone_index_correct_correct_error(:,1) = use_p_index(use_sig_neuron,20);
p_preCorrectTone_index_correct_correct_error(:,2) = use_p_index(use_sig_neuron,21);

preCorrectChoice_index_correct_correct_error(:,1) = use_neuron_index(use_sig_neuron,22);
preCorrectChoice_index_correct_correct_error(:,2) = use_neuron_index(use_sig_neuron,23);
p_preCorrectChoice_index_correct_correct_error(:,1) = use_p_index(use_sig_neuron,22);
p_preCorrectChoice_index_correct_correct_error(:,2) = use_p_index(use_sig_neuron,23);


pre_same_pLcL_pLcH_index(:,1) = use_neuron_index(use_sig_neuron,24);
pre_same_pHcL_pHcH_index(:,1) = use_neuron_index(use_sig_neuron,25);
cur_same_pLcL_pHcL_index(:,1) = use_neuron_index(use_sig_neuron,26);
cur_same_pLcH_pHcH_index(:,1) = use_neuron_index(use_sig_neuron,27);

p_pre_same_pLcL_pLcH_index(:,1) = use_p_index(use_sig_neuron,24);
p_pre_same_pHcL_pHcH_index(:,1) = use_p_index(use_sig_neuron,25);
p_cur_same_pLcL_pHcL_index(:,1) = use_p_index(use_sig_neuron,26);
p_cur_same_pLcH_pHcH_index(:,1) = use_p_index(use_sig_neuron,27);


norm_spike = use_norm_spike(use_sig_neuron,:);
p_surprise = use_p_surprise(use_sig_neuron,:);

temp = norm_spike(:,2) - norm_spike(:,1); %right-left
neuron_choice_category = zeros(length(temp),1);
temp = find(temp > 0);
neuron_choice_category(temp) = 1; %right=1, left=0

%Remove NAN
test = [pre_sound, norm_spike, choice_correct_error,tone_correct_error,...
    pre_tone_correct_error,pre_choice_correct_error,SameDif_preTone_correct_error, p_SameDif_preTone_correct_error,...
    preCorrectTone_index_correct_correct_error,p_preCorrectTone_index_correct_correct_error,...
    preCorrectChoice_index_correct_correct_error,p_preCorrectChoice_index_correct_correct_error,...
    pre_same_pLcL_pLcH_index, p_pre_same_pLcL_pLcH_index,...
    pre_same_pHcL_pHcH_index,p_pre_same_pHcL_pHcH_index,...
    cur_same_pLcL_pHcL_index, p_cur_same_pLcL_pHcL_index,...
    cur_same_pLcH_pHcH_index,p_cur_same_pLcH_pHcH_index ];
test = mean(test,2);
test = isnan(test);
test_nan = find(test == 1);
pre_sound(test_nan,:) = [];
p_pre_sound(test_nan,:) = [];
norm_spike(test_nan,:) = [];
p_surprise(test_nan,:) = [];
choice_correct_error(test_nan,:) = [];
p_choice_correct_error(test_nan,:) = [];
neuron_choice_category(test_nan,:) = [];

tone_correct_error(test_nan,:) = []; %Correct choice index
p_tone_correct_error(test_nan,:) = [];

pre_tone_correct_error(test_nan,:) = []; %Correct choice index
p_pre_tone_correct_error(test_nan,:) = [];

pre_choice_correct_error(test_nan,:) = [];
p_pre_choice_correct_error(test_nan,:) = [];

SameDif_preTone_correct_error(test_nan,:) = [];
p_SameDif_preTone_correct_error(test_nan,:) = [];

preCorrectTone_index_correct_correct_error(test_nan,:) = [];
p_preCorrectTone_index_correct_correct_error(test_nan,:) = [];

preCorrectChoice_index_correct_correct_error(test_nan,:) = [];
p_preCorrectChoice_index_correct_correct_error(test_nan,:) = [];

pre_same_pLcL_pLcH_index(test_nan,:) = [];
p_pre_same_pLcL_pLcH_index(test_nan,:) = [];
pre_same_pHcL_pHcH_index(test_nan,:) = [];
p_pre_same_pHcL_pHcH_index(test_nan,:) = [];
cur_same_pLcL_pHcL_index(test_nan,:) = [];
p_cur_same_pLcL_pHcL_index(test_nan,:) = [];
cur_same_pLcH_pHcH_index(test_nan,:) = [];
p_cur_same_pLcH_pHcH_index(test_nan,:) = [];


return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sig_neuron = get_sig_neuron_time_window(p_task,new_p_thre,window)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
use_p_sound = min(p_task(:,window),[],2);
use_p_sound = -log10(use_p_sound);
sig_neuron = find(use_p_sound > new_p_thre);

return





