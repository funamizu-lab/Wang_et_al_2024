%{
----------------------------------------------------------------------------
%Stimulus bias changes by block
%Just model analysis
%Based on Rao 2010 Front Comp Neurosci
%Value updating + Prior updating
----------------------------------------------------------------------------
%}

function Figure2D_top_ExampleSessionRL_20240115(root_dir)

root_dir = '/Volumes/Extreme SSD';

path1 = [root_dir '/Zigzag_repeat_Ephys/Behavior_during_learning/Repeat_first/W14/repeat/W14L_ToneClouds_HMM3_Rep_20210330_Aug11_2021_Session2'];
path2 = [root_dir '/Zigzag_repeat_Ephys/Behavior_during_learning/Repeat_first/W14/repeat'];

cd(path1)
filename1 = dir('Bpod_mat*');
filename1 = filename1.name;
cd(path2)
filename2 = dir('HMM_RL_20231128_indivi_session*');
data = load(filename2.name);

use_session_number= 7;

tic

[Choice_trial,tone_evidence,~,~,...
    low,high,~,~,~,~,~,...
    binary_tone,~,~,~,~] ...
 = HMM_get_basic_task_structure_20240714(filename1,path1);

load([path1 '/' filename1])

trial_use = 2:length(Choice_trial);
PostLow  = low+1;
PostHigh = high+1;
PostLow = intersect(PostLow,trial_use);
PostHigh = intersect(PostHigh,trial_use);

% %Take out chose success trials:

%Use all trials for analysis
N_trial = length(Choice_trial);
Correct_trial = Correct_side(Choice_trial);
Chosen_trial = Chosen_side(Choice_trial);
binary_tone = binary_tone(Choice_trial);
Reward_LCR = Reward_LCR(Choice_trial,[1,3]); %choose chosen trials
TrialBlock = TrialBlock(Choice_trial);

clear init_Q
init_Q = BlockReward(1,:); %2.4 2.4
    
%reward amount for each trial

reward_trial = zeros(N_trial,1);
for i = 1:N_trial
    %Reward trials
    if Correct_trial(i) == Chosen_trial(i)
        reward_trial(i) = Reward_LCR(i,Chosen_trial(i)+1);
    end
end

%Parameter
% bin_x = [0:0.001:1];
% stim_LR   = [0,0.25,0.45,0.55,0.75,1]; %3stim
% stim_prob_LR   = [0.5,0.25,0.25,0.25,0.25,0.5]; %3stim
% stim_prob_kitei = [0.25,0.125,0.125,0.125,0.125,0.25]; %3stim
%Sense_std = [0.01:0.05:1];
%reward = [3,1]; %left-right
%stim_prob = [0.5, 0.5];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%
%Based on the binary tone decide the pseudo tone evidence
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

mean_BIC = mean(data.BIC_all);
use_parameter = mean_BIC == min(mean_BIC); 
% use_parameter = 5;
temp_use_para = [1 2 3 4];
temp = data.para_max(use_session_number).matrix;
para_max = temp(use_parameter,:);


%%
disp('start simulation')
% Simulate choice and make psychometric function

parfor i = 1:100
    disp([i,100])

    [simu_action(i,:),~,~,~] = ...
        HMM_Thre_update_211030_simulate(para_max, Reward_LCR, binary_tone, init_Q, temp_use_para);


    [opt_L_max(i,:),opt_R_max(i,:),~,~] = ...
    simulated_choice_full_psycho2(simu_action(i,:),trial_evidence, tone_evidence,binary_tone,PostLow,PostHigh);


end
delete(gcp('nocreate'))

%% Plot psychometric
evi_x = 0:0.01:1;
figure
hold on;%Max trials
plot_mean_se_moto_x_axis(opt_L_max,evi_x,[0 0 1],1);
plot_mean_se_moto_x_axis(opt_R_max,evi_x,[1 0 0],1);
set(gca,'xlim',[-0.1 1.1],'ylim',[0 1])

toc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [opt_L_max,opt_R_max,right_maxL,right_maxR] = ...
    simulated_choice_full_psycho2(simu_choice,trial_evidence, tone_evidence,binary_tone,max_L,max_R)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Make artificial choice based on the Choice_prob
simu_choice = simu_choice - 1;

[right_maxL, ~] = get_right_prob(trial_evidence,simu_choice,tone_evidence,max_L);
[right_maxR, ~] = get_right_prob(trial_evidence,simu_choice,tone_evidence,max_R);
%%%%%%%%%%%%%%%%%

[opt_L_max,opt_R_max] = get_Gauss_standard_fit2(max_L, max_R, ...
    simu_choice, binary_tone, right_maxL, right_maxR);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [min_right,min_right05] = get_right_prob(trial_evidence,simu_choice,tone_evidence,minD_trial)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Min trials
min_evi = trial_evidence(minD_trial);
min_choice = simu_choice(minD_trial);
min_right = nan(length(tone_evidence),1);
for i = 1:length(tone_evidence)
    temp_choice = min_choice(min_evi == tone_evidence(i));    
    [min_right(i),min_right05(i,:)] = binofit(sum(temp_choice),length(temp_choice));
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sim_action,Q,Choice_prob,para] = ...
    HMM_Thre_update_211030_simulate(para, reward, Stim_freq, init_Q, use_para)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% para(1) = 0.1; % forgetting value
% para(2) = 0.2;  % std
% para(3) = 0.3683;  % beta
% para(4) = 0;  % bias
% %para(5) = 0.1;  %inverse temperature for prior
% use_para(1).matrix = [2];
% use_para(2).matrix = [1 2];
% use_para(3).matrix = [1 2 3];
% use_para(4).matrix = [1 2 4];
% use_para(5).matrix = [1 2 3 4];

temp = find(use_para == 1, 1);
if isempty(temp)
    para(1) = 0;
end
temp = find(use_para == 3, 1);
if isempty(temp)
    softmax_check = 0;
else
    softmax_check = 1;
end
temp = find(use_para == 4, 1);
if isempty(temp)
    para(4) = 0;
end

%para shusei
if para(1) < 0
    para(1) = 0;
elseif para(1) > 1
    para(1) = 1;
end
if para(2) < 0
    para(2) = 0;
end
if para(3) < 0
    para(3) = 0;
end

%Plot the gaussian based class for left and right
%Write the figures for bayes computation
N_trial = length(reward);

clear Prior Q_left Q_right Choice_prob
%Belief_basis for left choice reward probability
Q = zeros(N_trial,2);
Q(1,:)  = init_Q;

Choice_prob = zeros(N_trial,2);

sim_action = nan(N_trial,1);
for i = 1:N_trial
    clear Posterior_LR Decision_LR
    
    if softmax_check
        choice_prior = softmax(Q(i,:),para(3));
    else
        choice_prior = Q(i,:)./sum(Q(i,:));
    end
    decision_x = choice_prior(1);

    left_choice  = normcdf(decision_x + para(4),Stim_freq(i),para(2)) - normcdf(0,Stim_freq(i),para(2));
    right_choice = normcdf(1,Stim_freq(i),para(2)) - normcdf(decision_x + para(4),Stim_freq(i),para(2));
    left_choice = left_choice / (left_choice + right_choice);
    right_choice = 1 - left_choice;

    if left_choice > 1
        left_choice = 1;
        right_choice = 0;
    elseif right_choice > 1
        left_choice = 0;
        right_choice = 1;
    end

    if left_choice > rand
        sim_action(i) = 1; %left
    else
        sim_action(i) = 2; %right
    end
    
    Choice_prob(i,:) = [left_choice,right_choice];
    %Get probability from action and reward;
%     likelihood(i) = Choice_prob(i,action(i));
    
    %Forgetting Q learning
    Q(i+1,:) = Qlearning_ori(Q(i,:), sim_action(i), reward(i), para(1));
end

%likelihood keisan
% log_likeli = sum(log(likelihood));
% log_likeli = log_likeli / N_trial;
% ave_likeli = exp(log_likeli);
% ave_likeli = -ave_likeli;

return

function post_prob = Qlearning_ori(Q, choice, reward, alpha)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Control alpha
temp_Q = Q(:,choice) + alpha .* (reward - Q(:,choice)); 

if choice == 1
    post_prob = [temp_Q, (1-alpha) * Q(:,2)];
else
    post_prob = [(1-alpha) * Q(:,1), temp_Q];
end
return

function stim_prob = softmax(temp_stim,beta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear stime_prob
stim_prob(1) = exp(beta * temp_stim(1)) / (exp(beta * temp_stim(1)) + exp(beta * temp_stim(2)));
stim_prob(2) = exp(beta * temp_stim(2)) / (exp(beta * temp_stim(1)) + exp(beta * temp_stim(2)));

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [opt_L,opt_R,X_L,X_R] = get_Gauss_standard_fit2(block_L, block_R, ...
    Chosen_side, trial_evidence, right_prob_L, right_prob_R)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

evi_x = 0:0.01:1;

lapse_L = [right_prob_L(1), 1-right_prob_L(6)]; %limit for lapse
lapse_R = [right_prob_R(1), 1-right_prob_R(6)]; %limit for lapse

opt = optimset('Display','off');
para = [0.5 0.01 0 0
        0.5 0.1  0 0
        0.5 0.2  0 0
        0.5 0.5  0 0
        0.5 10    0 0
        0.5 0.01 0.1 0.1
        0.5 0.1  0.1 0.1
        0.5 0.2  0.1 0.1
        0.5 0.5  0.1 0.1
        0.5 10    0.1 0.1];

for i = 1:10
    [X_L(i,:),FCAL_L(i)] = fminsearch(@Opt_psychometric_Gauss,para(i,:),opt,Chosen_side(block_L), trial_evidence(block_L), evi_x, lapse_L);
    [X_R(i,:),FCAL_R(i)] = fminsearch(@Opt_psychometric_Gauss,para(i,:),opt,Chosen_side(block_R), trial_evidence(block_R), evi_x, lapse_R);
    %data fitting is not least square but with maximum likelihood
end

min_L = find(FCAL_L == min(FCAL_L),1);
min_R = find(FCAL_R == min(FCAL_R),1);
X_L = X_L(min_L,:);
X_R = X_R(min_R,:);

[~,opt_L,X_L] = Opt_psychometric_Gauss_max(X_L, Chosen_side(block_L), trial_evidence(block_L), evi_x, lapse_L);
[~,opt_R,X_R] = Opt_psychometric_Gauss_max(X_R, Chosen_side(block_R), trial_evidence(block_R), evi_x, lapse_R);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function log_likeli = Opt_psychometric_Gauss(para, chosen_side, tone_evi, evi_x, lapse_limit)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[log_likeli,~] = Opt_psychometric_Gauss_max(para, chosen_side, tone_evi, evi_x, lapse_limit);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [log_likeli,neurometric,para] = Opt_psychometric_Gauss_max(para, chosen_side, tone_evi, evi_x, lapse_limit)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%para(1): threthold: should be around 0.5
%para(2): standard deviation of Gaussian: from 0.01 to 1 
%para(3): lambda1
%para(4): lambda2

%4 parameters
%Yn = lambda1 + (1-lambda1-lambda2) .* norminv(x,std)

if para(1) < 0 %para(2) shoud be positive
    para(1) = 0;
elseif para(1) > 1
    para(1) = 1;
end
if para(2) < 0 %para(2) shoud be positive
    para(2) = eps;
end
if para(3) < 0
    para(3) = 0;
elseif para(3) > lapse_limit(1)
    para(3) = lapse_limit(1);
end
if para(4) < 0
    para(4) = 0;
elseif para(4) > lapse_limit(2)
    para(4) = lapse_limit(2);
end

N_trial = size(tone_evi,1);
likelihood = zeros(1,N_trial);

%Trial by trial, get the likelihood
clear temp_p temp_exp
temp_right = normcdf(1, tone_evi, para(2)) - normcdf(para(1), tone_evi, para(2)); %Right choice probablity
temp_left  = normcdf(para(1), tone_evi, para(2)) - normcdf(0, tone_evi, para(2)); %Left choice probablity
temp_gauss = temp_right ./ (temp_right + temp_left); %normalized with truncated part
temp_p = para(3) + (1-para(3)-para(4)) .* temp_gauss;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
temp1 = find(chosen_side == 0);
temp2 = find(chosen_side == 1);
likelihood(temp1) = 1-temp_p(temp1); %left choice
likelihood(temp2) = temp_p(temp2);

if length(temp1)+length(temp2) ~= N_trial
    disp([length(temp1),length(temp2),N_trial])
    hoge
end

%likelihood keisan
log_likeli = sum(log(likelihood));
log_likeli = -log_likeli;

%get the tuning function with evi_x
temp_right = normcdf(1, evi_x, para(2)) - normcdf(para(1), evi_x, para(2)); %Right choice probablity
temp_left  = normcdf(para(1), evi_x, para(2)) - normcdf(0, evi_x, para(2)); %Left choice probablity
temp_gauss = temp_right ./ (temp_right + temp_left); %normalized with truncated part
neurometric = para(3) + (1-para(3)-para(4)) .* temp_gauss;

return


%{
----------------------------------------------------------------------------
Analyzing behavioral data
At least for the correct rate
----------------------------------------------------------------------------
%}

function [Choice_trial,tone_evidence,trial_evidence,use_trial_remove_first,...
    low,high,correct,error,flip_tone,number_use_trial,number_use_trial_remove_first,...
    binary_tone,right_trial_all,number_trial_all,right_trial,number_trial] ...
 = HMM_get_basic_task_structure_20240714(filename1,path1)

load([path1 '/' filename1])


check_duration = unique(StimDuration);
if length(check_duration) ~= 1
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

temp_evi = unique(EvidenceStrength);%
temp_evi_low  = 0.5 - temp_evi/2;%
temp_evi_high = 0.5 + temp_evi/2;%
temp_evi_all = [temp_evi_low', temp_evi_high'];%
tone_evidence = sort(temp_evi_all);%

%Put tone evidence in all trials;
trial_evidence = zeros(length(Outcome),1);%
low  = find(Correct_side == 0);
high = find(Correct_side == 1);
temp = Chosen_side == Correct_side; 
correct = find(temp == 1);
error   = find(temp == 0);
for i = 1:length(temp_evi)
    temp = find(EvidenceStrength == temp_evi(i));% 
    temp_left  = intersect(temp,low); %
    temp_right = intersect(temp,high);%
    trial_evidence(temp_left)  = temp_evi_low(i);%
    trial_evidence(temp_right) = temp_evi_high(i); % 
end

% left  = 0;
% right = 1;
% undefined = 2;

% low   = 0;
% high  = 1;

% lowLeft = 1;
% lowRight = 2;
% centerValveCode = 2;

%Make the true tone cloud value
binary_tone = nan(length(Tone_cloud),1);
for i = 1:length(Tone_cloud)
    temp_tone = Tone_cloud(i).matrix;
    %Get the data in all sound
    temp1 = find(temp_tone >= 9);
    binary_tone(i,1) = length(temp1) ./ length(temp_tone);  
end


% if max_block ~= 1
%     hoge
% end
use_trial = Choice_trial;
use_trial_remove_first = use_trial(InitBlock+1:length(use_trial));
number_use_trial = length(use_trial);
number_use_trial_remove_first = length(use_trial_remove_first);
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
% disp(temp_tone(2).matrix)
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
    disp([binary_tone, new_tone_evi]);
    disp('nan detected')
    hoge
end
trial_evidence = new_tone_evi;
PostLow  = low+1;
PostHigh = high+1;

low_trial = intersect(PostLow, use_trial_remove_first);
high_trial = intersect(PostHigh, use_trial_remove_first);

%Get the number of right choice trials in each session
if BlockProb(1) > 0.5 % Zig
    [right_trial.zig, number_trial.zig] = get_right_choice_trials(use_trial_remove_first,Chosen_side,trial_evidence,tone_evidence);
    [right_trial.Z_low, number_trial.Z_low] = get_right_choice_trials(low_trial,Chosen_side,trial_evidence,tone_evidence);
    [right_trial.Z_high, number_trial.Z_high] = get_right_choice_trials(high_trial,Chosen_side,trial_evidence,tone_evidence);
    right_trial.stay = [];
    number_trial.stay = [];
    right_trial.S_low = [];
    number_trial.S_low = [];
    right_trial.S_high = [];
    number_trial.S_high = [];
    
    right_trial_all = right_trial.zig;
    number_trial_all = number_trial.zig;
else %Stay
    [right_trial.stay, number_trial.stay] = get_right_choice_trials(use_trial_remove_first,Chosen_side,trial_evidence,tone_evidence);
    [right_trial.S_low, number_trial.S_low] = get_right_choice_trials(low_trial,Chosen_side,trial_evidence,tone_evidence);
    [right_trial.S_high, number_trial.S_high] = get_right_choice_trials(high_trial,Chosen_side,trial_evidence,tone_evidence);
    right_trial.zig = [];
    number_trial.zig = [];
    right_trial.Z_low = [];
    number_trial.Z_low = [];
    right_trial.Z_high = [];
    number_trial.Z_high = [];
    
    right_trial_all = right_trial.stay;
    number_trial_all = number_trial.stay;
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [right_trial, number_trial] = get_right_choice_trials(use_trials,Chosen_side,trial_evidence,tone_evidence)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

use_choice = Chosen_side(use_trials);
use_tone = trial_evidence(use_trials);
number_trial = nan(length(tone_evidence),1);
right_trial = nan(length(tone_evidence),1);
for i = 1:length(tone_evidence)
    temp_trial = find(use_tone == tone_evidence(i));
    temp_choice = use_choice(temp_trial);
    number_trial(i) = length(temp_trial);
    right_trial(i) = sum(temp_choice);
end

return



