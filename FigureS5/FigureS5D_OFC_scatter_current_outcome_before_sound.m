function FigureS5D_OFC_scatter_current_outcome_before_sound


process_CurrentOutcome_20240703_surprise_before_sound_OFC('repeat_OFC_20240427',1);
process_CurrentOutcome_20240703_surprise_before_sound_OFC('altern_OFC_20240427',1);

function process_CurrentOutcome_20240703_surprise_before_sound_OFC(folders, kaiseki_number)


Brainarea = [folders(1:6),' ', folders(8:10)];
[analysis_dir,depth_def] = eval(folders);

all_pre_sound = [];
all_p_pre_sound = [];
all_norm_spike = [];
all_p_surprise = [];

all_pre_choice_correct_error = [];
all_p_pre_choice_correct_error = [];

for i = 1:length(analysis_dir)
    disp([i,length(analysis_dir)])
    temp_dir = analysis_dir{i};
    
    [pre_sound, p_pre_sound, norm_spike, ~, ~, ~, ~, ~, ~, ~, ~,...
    pre_choice_correct_error, p_pre_choice_correct_error,~,~,~,~,~,~,~,~,~,~,~,~,~,~] = ...
        HMM_ephys_20240401_surprise_before_sound_depth(temp_dir, kaiseki_number,depth_def);

%     p_pre_sound(:,1) = use_p_index(use_sig_neuron,1); %current choice index
%     p_pre_sound(:,2) = use_p_index(use_sig_neuron,6); %pre tone index
%     p_pre_sound(:,3) = use_p_index(use_sig_neuron,7); %pre choice index
%     p_pre_sound(:,4) = use_p_index(use_sig_neuron,10); %pre SoundCorrect, current left
%     p_pre_sound(:,5) = use_p_index(use_sig_neuron,11); %pre SoundCorrect, current right
%     p_pre_sound(:,6) = use_p_index(use_sig_neuron,11); % current tone index

    all_pre_sound = [all_pre_sound; pre_sound];
    all_p_pre_sound = [all_p_pre_sound; p_pre_sound];
    all_norm_spike = [all_norm_spike; norm_spike];

    all_pre_choice_correct_error = [all_pre_choice_correct_error; pre_choice_correct_error]; % pre choice correct index
    all_p_pre_choice_correct_error = [all_p_pre_choice_correct_error; p_pre_choice_correct_error]; % pre choice error index
end
delete(gcp('nocreate'))

disp([length(all_pre_sound),length(all_norm_spike),length(all_p_surprise)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if contains(Brainarea,'repeat')
    right_neuron = find(all_pre_choice_correct_error(:,1) >= 0);
    left_neuron = find(all_pre_choice_correct_error(:,1) < 0);
    sig_sound = find(all_p_pre_sound(:,3) < 0.01); %Sig diff. between previous low and high
else
    left_neuron = find(all_pre_choice_correct_error(:,1) >= 0);
    right_neuron = find(all_pre_choice_correct_error(:,1) < 0);
    sig_sound = find(all_p_pre_sound(:,3) < 0.01); %Sig diff. between previous low and high

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

%Use only on low sig neuron
[L_sig_PostLow_leftcorrect, L_sig_PostHigh_leftcorrect, L_sig_PostLow_rightcorrect, L_sig_PostHigh_rightcorrect,...
    L_sig_PostLow_lefterror, L_sig_PostHigh_lefterror, L_sig_PostLow_righterror, L_sig_PostHigh_righterror] = ...
    test_surprise_choice(all_norm_spike, left_sig_neuron);
%Use only on high sig neuron
[R_sig_PostLow_leftcorrect, R_sig_PostHigh_leftcorrect, R_sig_PostLow_rightcorrect, R_sig_PostHigh_rightcorrect,....
    R_sig_PostLow_lefterror, R_sig_PostHigh_lefterror, R_sig_PostLow_righterror, R_sig_PostHigh_righterror] = ...
    test_surprise_choice(all_norm_spike, right_sig_neuron);

%Activity is already flipped

if contains(Brainarea,'repeat')
    Repeat_prefer_correct = [L_sig_PostLow_leftcorrect; R_sig_PostHigh_rightcorrect];
    Repeat_nonprefer_correct = [L_sig_PostHigh_rightcorrect; R_sig_PostLow_leftcorrect];
    Switch_prefer_correct = [L_sig_PostLow_rightcorrect; R_sig_PostHigh_leftcorrect];
    Switch_nonprefer_correct = [L_sig_PostHigh_leftcorrect; R_sig_PostLow_rightcorrect];

    Repeat_prefer_error = [L_sig_PostLow_lefterror; R_sig_PostHigh_righterror];
    Repeat_nonprefer_error = [L_sig_PostHigh_righterror; R_sig_PostLow_lefterror];
    Switch_prefer_error = [L_sig_PostLow_righterror; R_sig_PostHigh_lefterror];
    Switch_nonprefer_error = [L_sig_PostHigh_lefterror; R_sig_PostLow_righterror];
else
    
    Repeat_prefer_correct = [L_sig_PostHigh_rightcorrect; R_sig_PostLow_leftcorrect];
    Repeat_nonprefer_correct = [L_sig_PostLow_leftcorrect; R_sig_PostHigh_rightcorrect];
    Switch_prefer_correct = [L_sig_PostHigh_leftcorrect; R_sig_PostLow_rightcorrect];
    Switch_nonprefer_correct = [L_sig_PostLow_rightcorrect; R_sig_PostHigh_leftcorrect];


    Repeat_prefer_error = [L_sig_PostHigh_righterror; R_sig_PostLow_lefterror];
    Repeat_nonprefer_error = [L_sig_PostLow_lefterror; R_sig_PostHigh_righterror];
    Switch_prefer_error = [L_sig_PostHigh_lefterror; R_sig_PostLow_righterror];
    Switch_nonprefer_error = [L_sig_PostLow_righterror; R_sig_PostHigh_lefterror];


   
end

p_prefer_correct = signrank(Repeat_prefer_correct,Switch_prefer_correct);
p_nonprefer_correct = signrank(Repeat_nonprefer_correct,Switch_nonprefer_correct);
p_prefer_error = signrank(Repeat_prefer_error,Switch_prefer_error);
p_nonprefer_error = signrank(Repeat_nonprefer_error,Switch_nonprefer_error);

p_sabun = [p_prefer_correct, p_nonprefer_correct, p_prefer_error, p_nonprefer_error];

temp1 = length((Repeat_prefer_correct - Switch_prefer_correct));
temp2 = length((Repeat_nonprefer_correct - Switch_nonprefer_correct));

if temp1 ~= temp2
    hoge
end


length_L_sig = length(left_sig_neuron);

figure
subplot(2,2,1)
plot_surprise_choice_prefer(Repeat_prefer_correct(1:length_L_sig),Switch_prefer_correct(1:length_L_sig),[0 0 1]);
hold on
plot_surprise_choice_prefer(Repeat_prefer_correct(length_L_sig+1:end),Switch_prefer_correct(length_L_sig+1:end), [1 0 0]);
title('Preferred choice -- Current Correct')
box off

subplot(2,2,2)

plot_surprise_choice_prefer(Repeat_prefer_error(1:length_L_sig),Switch_prefer_error(1:length_L_sig),[0 0 1]);
hold on
plot_surprise_choice_prefer(Repeat_prefer_error(length_L_sig+1:end),Switch_prefer_error(length_L_sig+1:end),[1 0 0]);
title('Preferred choice -- Current Error')
box off


subplot(2,2,3)
plot_surprise_choice_nonprefer(Repeat_nonprefer_correct(1:length_L_sig),Switch_nonprefer_correct(1:length_L_sig),[0 0 1]);
hold on
plot_surprise_choice_prefer(Repeat_nonprefer_correct(length_L_sig+1:end),Switch_nonprefer_correct(length_L_sig+1:end), [1 0 0]);
title('Nonpreferred choice -- Current Correct')
box off

subplot(2,2,4)
plot_surprise_choice_nonprefer(Repeat_nonprefer_error(1:length_L_sig),Switch_nonprefer_error(1:length_L_sig),[0 0 1]);
hold on
plot_surprise_choice_nonprefer(Repeat_nonprefer_error(length_L_sig+1:end),Switch_nonprefer_error(length_L_sig+1:end),[1 0 0]);
title('Nonpreferred choice -- Current Error')
sgtitle(Brainarea)
set(gca,'FontName', 'Arial')

box off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [PostLow_Leftcorrect, PostHigh_Leftcorrect, PostLow_Rightcorrect, PostHigh_Rightcorrect,...
    PostLow_Lefterror,PostHigh_Lefterror,PostLow_Righterror,PostHigh_Righterror] = ...
    test_surprise_choice(all_norm_spike, low_sig_neuron)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PostLow_Leftcorrect = all_norm_spike(low_sig_neuron,33); %PostLowCorrect_left
PostHigh_Leftcorrect = all_norm_spike(low_sig_neuron,37); %PostHighCorrect_left
PostLow_Rightcorrect = all_norm_spike(low_sig_neuron,35); %PostLowCorrect_right
PostHigh_Rightcorrect = all_norm_spike(low_sig_neuron,39); %PostHighCorrect_right

PostLow_Lefterror = all_norm_spike(low_sig_neuron,34); %PostLowCorrect_left
PostHigh_Lefterror = all_norm_spike(low_sig_neuron,38); %PostHighCorrect_left
PostLow_Righterror = all_norm_spike(low_sig_neuron,36); %PostLowCorrect_right
PostHigh_Righterror = all_norm_spike(low_sig_neuron,40); %PostHighCorrect_right

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

