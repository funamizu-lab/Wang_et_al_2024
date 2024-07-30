%{
----------------------------------------------------------------------------
%Stimulus bias changes by block
%Just model analysis
%Based on Rao 2010 Front Comp Neurosci
%Value updating + Prior updating
----------------------------------------------------------------------------
%}

function [ave_likeli,BIC_all,log_likeli,para_max,N_trial] = ...
    HMM_RL_model_20211029_model_only_each_session(Correct_trial, Chosen_trial, binary_tone, Reward_LCR)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Use all trials for analysis
N_trial = length(Correct_trial);
% Correct_trial = Correct_side(Choice_trial);
% Chosen_trial = Chosen_side(Choice_trial);
% binary_tone = binary_tone(Choice_trial);
% Reward_LCR = Reward_LCR(Choice_trial,[1,3]); %choose chosen trials

clear init_Q
init_Q = Reward_LCR(1,:); %2.4 2.4
    
%reward amount for each trial

reward_trial = zeros(N_trial,1);
for i = 1:N_trial
    %Reward trials
    if Correct_trial(i) == Chosen_trial(i)
        reward_trial(i) = Reward_LCR(i,Chosen_trial(i)+1);
    end
end

%Parameter
bin_x = [0:0.001:1];
stim_LR   = [0,0.25,0.45,0.55,0.75,1]; %3stim
stim_prob_LR   = [0.5,0.25,0.25,0.25,0.25,0.5]; %3stim
stim_prob_kitei = [0.25,0.125,0.125,0.125,0.125,0.25]; %3stim
%Sense_std = [0.01:0.05:1];
%reward = [3,1]; %left-right
%stim_prob = [0.5, 0.5];

% %Parameter
% para(1) = 0.1; % forgetting value
% para(2) = 0.2;  % std
% para(3) = 0.3683;  % beta
% para(4) = 0;  % bias
% %para(5) = 0.1;  %inverse temperature for prior

use_para(1).matrix = [2];
use_para(2).matrix = [1 2];
use_para(3).matrix = [1 2 3];
use_para(4).matrix = [1 2 4];
use_para(5).matrix = [1 2 3 4];

para_temp = [0   0.2 0.5 0;
             0   0.5 1   0;
             0   1   0.5 0;
             0.1 0.2 1   0;
             0.1 0.5 0.5 0;
             0.3 0.2 0.5 0;
             0.3 0.5 1 0;
             ];
[size_y,size_x] = size(para_temp);
         
%Get maximum prediction
opt = optimset('Display','off');
clear ave_likeli BIC_all log_likeli
for j = 1:length(use_para)
    j
    temp_use_para = use_para(j).matrix;
    parfor i = 1:size_y
        i
        [X1(i,:),FCAL1(i,:),EXITFLAG,OUTPUT] = ...
            fminsearch(@HMM_Thre_update_211030_all,para_temp(i,:), opt, reward_trial, Chosen_trial, binary_tone, bin_x, stim_LR, stim_prob_LR, init_Q, temp_use_para);
    end
    min_FCAL = find(FCAL1 == min(FCAL1),1);
    temp_para_max = X1(min_FCAL,:);

    [ave_likeli(1,j),likelihood,Q,Choice_prob,para_max(j,:)] = ...
        HMM_Thre_update_211030_max_all(temp_para_max, reward_trial, Chosen_trial, binary_tone, bin_x, stim_LR, stim_prob_LR, init_Q, temp_use_para);
    
    sumlog_likeli = sum(log(likelihood));
    log_likeli(1,j) = sumlog_likeli;
    BIC_all(1,j) = -2 * sum(log(likelihood)) + length(temp_use_para) * log(N_trial);
        %log_likelihood(filecount).matrix = log(likelihood);
end
ave_likeli
BIC_all
log_likeli
para_max
N_trial

% figure
% plot(Q(:,1),'b')
% hold on
% plot(Q(:,2),'r')
% set(gca,'xlim',[1 N_trial])
% set(gca,'xtick',[0:40:N_trial])
% set(gca,'ylim',[0 4])
% 
% save_file = 'HMM_all_session_RL_20211029.mat';
% save(save_file, 'ave_likeli' ,'BIC_all','log_likeli','para_max','N_trial')

% delete(gcp('nocreate'))
%hoge

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ave_likeli,likelihood,Q,Choice_prob,para] = ...
    HMM_Thre_update_211030_max_all(para, reward, action, Stim_freq, bin_x, stim_LR, stim_prob_LR, init_Q, use_para)
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
% %long stimulus
% for i = 1:length(stim_LR)
%     temp_x = normpdf(bin_x,stim_LR(i),para(2));
%     temp_x = temp_x ./ sum(temp_x);
%     LR_stim(i,:) = stim_prob_LR(i) .* temp_x;
%     norm_pdf_long(i,:) = temp_x;
% end
% %Combination of left and right
% L_long = sum(LR_stim([1:3],:));
% R_long = sum(LR_stim([4:6],:));

clear Prior Q_left Q_right Choice_prob
%Belief_basis for left choice reward probability
Q = zeros(N_trial,2);
Q(1,:)  = init_Q;
likelihood = zeros(N_trial,1);
Choice_prob = zeros(N_trial,2);
action = action + 1; %0,1 -> 1,2

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
    
    Choice_prob(i,:) = [left_choice,right_choice];
    %Get probability from action and reward;
    likelihood(i) = Choice_prob(i,action(i));
    
    %Forgetting Q learning
    Q(i+1,:) = Qlearning_ori(Q(i,:), action(i), reward(i), para(1));
end

%likelihood keisan
log_likeli = sum(log(likelihood));
log_likeli = log_likeli / N_trial;
ave_likeli = exp(log_likeli);
ave_likeli = -ave_likeli;

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ave_likeli_mean = HMM_Thre_update_211030_all(para, reward_trial, Chosen_trial, binary_tone, bin_x, stim_LR, stim_prob_LR, init_Q, use_para)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ave_likeli_mean,~] = ...
    HMM_Thre_update_211030_max_all(para, reward_trial, Chosen_trial, binary_tone, bin_x, stim_LR, stim_prob_LR, init_Q, use_para);
ave_likeli_mean
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stim_prob = softmax(temp_stim,beta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear stime_prob
stim_prob(1) = exp(beta * temp_stim(1)) / (exp(beta * temp_stim(1)) + exp(beta * temp_stim(2)));
stim_prob(2) = exp(beta * temp_stim(2)) / (exp(beta * temp_stim(1)) + exp(beta * temp_stim(2)));

return


