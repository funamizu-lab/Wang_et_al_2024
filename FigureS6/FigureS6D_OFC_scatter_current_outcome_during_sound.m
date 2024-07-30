function FigureS6D_OFC_scatter_current_outcome_during_sound


process_CurrentOutcome_20240705_surprise_at_sound_OFC('repeat_OFC_20240427',1);
process_CurrentOutcome_20240705_surprise_at_sound_OFC('altern_OFC_20240427',1);

function process_CurrentOutcome_20240705_surprise_at_sound_OFC(folders, kaiseki_number)

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


Brainarea = [folders(1:6),' ', folders(8:10)];

all_norm_spike = [];
all_p_surprise = [];
all_sound_correct_error = [];
all_tone_correct_error = [];
all_neuron_choice_category = [];

all_pre_sound = [];
all_p_pre_sound = [];
for i = 1:length(analysis_dir)
    disp([i,length(analysis_dir)])
    temp_dir = analysis_dir{i};
    
    [sound_neuron, p_sound_neuron, norm_spike, p_surprise, ...
    sound_correct_error, ~, neuron_choice_category,tone_correct_error,~] = ...
        HMM_sound_trans20240705_surprise_at_sound_depth(temp_dir, kaiseki_number,depth_def);
    
    all_pre_sound = [all_pre_sound; sound_neuron];
    all_p_pre_sound = [all_p_pre_sound; p_sound_neuron];
    all_norm_spike = [all_norm_spike; norm_spike];
    all_p_surprise = [all_p_surprise; p_surprise];
    all_sound_correct_error = [all_sound_correct_error; sound_correct_error];
    all_tone_correct_error = [all_tone_correct_error;tone_correct_error];
    all_neuron_choice_category = [all_neuron_choice_category; neuron_choice_category];
end
delete(gcp('nocreate'))

disp([length(all_pre_sound),length(all_norm_spike),length(all_p_surprise)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Based on focusing on the sound or choice, should change the analysis way

%Focus on choice
target_sig = 1;
target_choice = [1,2]; %left, right
target_surprise_both = [25,26,27,28];

target_surprise_correct = [33,37,35,39];
target_surprise_error = [34,38,36,40];


get_analysis_sound_or_choice(target_sig, target_choice, target_surprise_both,target_surprise_correct, target_surprise_error,...
    all_neuron_choice_category,all_pre_sound,all_p_pre_sound,all_sound_correct_error,all_tone_correct_error,all_norm_spike,Brainarea);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function get_analysis_sound_or_choice(target_sig, ~, ~,target_surprise_correct, target_surprise_error,...
    ~,all_pre_sound,all_p_pre_sound,all_sound_correct_error,~,all_norm_spike,Brainarea)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

right_neuron = find(all_sound_correct_error(:,1) > 0);
left_neuron = find(all_sound_correct_error(:,1) < 0);
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

[L_sig_PostLow_leftcorrect, L_sig_PostHigh_leftcorrect, L_sig_PostLow_rightcorrect, L_sig_PostHigh_rightcorrect] = ...
    test_surprise_sound_choice(all_norm_spike, left_sig_neuron, target_surprise_correct);
%Use only on high sig neuron
[R_sig_PostLow_leftcorrect, R_sig_PostHigh_leftcorrect, R_sig_PostLow_rightcorrect, R_sig_PostHigh_rightcorrect] = ...
    test_surprise_sound_choice(all_norm_spike, right_sig_neuron, target_surprise_correct);

[L_sig_PostLow_lefterror, L_sig_PostHigh_lefterror, L_sig_PostLow_righterror, L_sig_PostHigh_righterror] = ...
    test_surprise_sound_choice(all_norm_spike, left_sig_neuron, target_surprise_error);
%Use only on high sig neuron
[R_sig_PostLow_lefterror, R_sig_PostHigh_lefterror, R_sig_PostLow_righterror, R_sig_PostHigh_righterror] = ...
    test_surprise_sound_choice(all_norm_spike, right_sig_neuron, target_surprise_error);

%Activity is already flipped

Repeat_prefer_correct = [L_sig_PostLow_leftcorrect; R_sig_PostHigh_rightcorrect];
Repeat_nonprefer_correct = [L_sig_PostHigh_rightcorrect; R_sig_PostLow_leftcorrect];
Switch_prefer_correct = [L_sig_PostHigh_leftcorrect; R_sig_PostLow_rightcorrect];
Switch_nonprefer_correct = [L_sig_PostLow_rightcorrect; R_sig_PostHigh_leftcorrect];

Repeat_prefer_error = [L_sig_PostLow_lefterror; R_sig_PostHigh_righterror];
Repeat_nonprefer_error = [L_sig_PostHigh_righterror; R_sig_PostLow_lefterror];
Switch_prefer_error = [L_sig_PostHigh_lefterror; R_sig_PostLow_righterror];
Switch_nonprefer_error = [L_sig_PostLow_righterror; R_sig_PostHigh_lefterror];



p_prefer_correct = signrank(Repeat_prefer_correct,Switch_prefer_correct);
p_nonprefer_correct = signrank(Repeat_nonprefer_correct,Switch_nonprefer_correct);



p_prefer_error = signrank(Repeat_prefer_error,Switch_prefer_error);
p_nonprefer_error = signrank(Repeat_nonprefer_error,Switch_nonprefer_error);


p_sabun = [p_prefer_correct, p_nonprefer_correct, p_prefer_error, p_nonprefer_error];



temp1 = length(Repeat_prefer_correct);
temp2 = length(Repeat_nonprefer_correct);

if temp1 ~= temp2
    hoge
end

figure
subplot(2,2,1)
plot_surprise_choice_prefer(L_sig_PostLow_leftcorrect, L_sig_PostHigh_leftcorrect,[0 0 1]);
hold on
plot_surprise_choice_prefer(R_sig_PostHigh_rightcorrect, R_sig_PostLow_rightcorrect, [1 0 0]);
title('Preferred choice -- Current Correct')

subplot(2,2,2)

plot_surprise_choice_prefer(L_sig_PostLow_lefterror, L_sig_PostHigh_lefterror,[0 0 1]);
hold on
plot_surprise_choice_prefer(R_sig_PostHigh_righterror, R_sig_PostLow_righterror, [1 0 0]);
title('Preferred choice -- Current Error')


subplot(2,2,3)
plot_surprise_choice_nonprefer(L_sig_PostHigh_rightcorrect, L_sig_PostLow_rightcorrect,[0 0 1]);
hold on
plot_surprise_choice_prefer(R_sig_PostLow_leftcorrect, R_sig_PostHigh_leftcorrect, [1 0 0]);
title('Nonpreferred choice -- Current Correct')

subplot(2,2,4)
plot_surprise_choice_nonprefer(L_sig_PostHigh_righterror, L_sig_PostLow_righterror,[0 0 1]);
hold on
plot_surprise_choice_prefer(R_sig_PostLow_lefterror, R_sig_PostHigh_lefterror, [1 0 0]);
title('Nonpreferred choice -- Current Error')
sgtitle(Brainarea)
set(gca,'FontName', 'Arial')

sgtitle(Brainarea)


return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [PostLow_left, PostHigh_left, PostLow_right, PostHigh_right] = ...
    test_surprise_sound_choice(all_norm_spike, low_sig_neuron, target_surprise)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PostLow_left = all_norm_spike(low_sig_neuron,target_surprise(1)); %PostLowCorrect_left
PostHigh_left = all_norm_spike(low_sig_neuron,target_surprise(2)); %PostHighCorrect_left
PostLow_right = all_norm_spike(low_sig_neuron,target_surprise(3)); %PostLowCorrect_right
PostHigh_right = all_norm_spike(low_sig_neuron,target_surprise(4)); %PostHighCorrect_right

function plot_surprise_choice_prefer(Repeat_prefer, Switch_prefer, plot_color)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% figure
% plot(predict_ave, choice_ave, 'k.')
plot(Repeat_prefer, Switch_prefer, '.','color', plot_color) %same choice, X_repeat or Y_zigzag
hold on
plot([-1, 10],[-1, 10],'k')
set(gca,'xlim',[-1, 10],'ylim',[-1, 10])


function plot_surprise_choice_nonprefer(Switch_nonprefer, Repeat_nonprefer, plot_color)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% figure
% plot(predict_ave, choice_ave, 'k.')

plot(Switch_nonprefer, Repeat_nonprefer, '.','color', plot_color)  %same choice, X_zigzag or Y_repeat
hold on
plot([-1, 10],[-1, 10],'k')
set(gca,'xlim',[-1, 10],'ylim',[-1, 10])

