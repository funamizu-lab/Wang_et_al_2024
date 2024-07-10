

function Figure1D_HMM_psycho_full_model_Gauss_20240710


[filename1, pathname1]=uigetfile('*.mat','Block_mat');
filename1 = [pathname1, filename1];
load(filename1)


[~,tone_evidence,trial_evidence,trial_remove_first,...
    low,high,~,~,~,~,~,...
    ~,~,~,right_trial,number_trial] ...
 = HMM_get_basic_task_structure_20210514(filename1);
    

PostLow  = low+1;
PostHigh = high+1;
block_use = trial_remove_first;
low_trial = intersect(PostLow, trial_remove_first);
high_trial = intersect(PostHigh, trial_remove_first);
conf_low = nan(6,2);
conf_high = nan(6,2);

if BlockProb(1) > BlockProb(2)
    right_prob = right_trial.zig ./ number_trial.zig;
    right_low = right_trial.Z_low ./ number_trial.Z_low;
    right_high = right_trial.Z_high ./ number_trial.Z_high;
    
    for i = 1:6
        % [~,conf(i,:)] = binofit(right_trial.zig(i),number_trial.zig(i));
        [~,conf_low(i,:)] = binofit(right_trial.Z_low(i),number_trial.Z_low(i));
        [~,conf_high(i,:)] = binofit(right_trial.Z_high(i),number_trial.Z_high(i));
    end
else 
    right_prob = right_trial.stay ./ number_trial.stay;
    right_low = right_trial.S_low ./ number_trial.S_low;
    right_high = right_trial.S_high ./ number_trial.S_high;
    
    for i = 1:6
        % [p(i,1),conf(i,:)] = binofit(right_trial.stay(i),number_trial.stay(i));
        [~,conf_low(i,:)] = binofit(right_trial.S_low(i),number_trial.S_low(i));
        [~,conf_high(i,:)] = binofit(right_trial.S_high(i),number_trial.S_high(i));
    end
end

if length(low_trial) + length(high_trial) ~= length(block_use)
    hoge
end

[~,~] = get_Gauss_standard_fit_20210514(block_use, Correct_side, Chosen_side, trial_evidence, right_prob);
[opt_low,~] = get_Gauss_standard_fit_20210514(low_trial, Correct_side, Chosen_side, trial_evidence, right_low);
[opt_high,~] = get_Gauss_standard_fit_20210514(high_trial, Correct_side, Chosen_side, trial_evidence, right_high);

evi_x = 0:0.01:1;
tick_sabun = 0.005;

figure
plot(evi_x,opt_low,'b','LineWidth',1)
hold on
plot(evi_x,opt_high,'r','LineWidth',1)
hold on

for i = 1:length(tone_evidence)
    plot([tone_evidence(i)-tick_sabun,tone_evidence(i)-tick_sabun],conf_low(i,:),'b','LineWidth',1)
    hold on
    plot([tone_evidence(i)+tick_sabun,tone_evidence(i)+tick_sabun],conf_high(i,:),'r','LineWidth',1)
    hold on
end
plot(tone_evidence-tick_sabun,right_low,'b.')
hold on
plot(tone_evidence+tick_sabun,right_high,'r.')

if BlockProb(1) > BlockProb(2) 
    title('ZigZag block')
else
    title('Repeat block')
end
set(gca,'xlim',[-0.1 1.1],'ylim',[0 1])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [opt_fit,X] = get_Gauss_standard_fit_20210514(block_use, ~, Chosen_side, trial_evidence, right_prob)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

evi_x = 0:0.01:1;

lapse = [right_prob(1), 1-right_prob(6)]; %limit for lapse

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

X = nan(10,4);
FCAL = nan(1,10);
for i = 1:10
    [X(i,:),FCAL(i),~,~] = fminsearch(@Opt_psychometric_Gauss,...
        para(i,:),opt,Chosen_side(block_use), trial_evidence(block_use), evi_x, lapse);
    %data fitting is not least square but with maximum likelihood
end

min_FCAL = find(FCAL == min(FCAL),1);
X = X(min_FCAL,:);

[~,opt_fit,X] = Opt_psychometric_Gauss_max(X, Chosen_side(block_use), trial_evidence(block_use), evi_x, lapse);

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

[N_trial,~] = size(tone_evi);
likelihood = zeros(1,N_trial);

%Trial by trial, get the likelihood
clear temp_p temp_exp
temp_right = normcdf(1, tone_evi, para(2)) - normcdf(para(1), tone_evi, para(2)); %Right choice probablity
temp_left  = normcdf(para(1), tone_evi, para(2))          - normcdf(0, tone_evi, para(2)); %Left choice probablity
temp_gauss = temp_right ./ (temp_right + temp_left); %normalized with truncated part
temp_p = para(3) + (1-para(3)-para(4)) .* temp_gauss;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
temp1 = find(chosen_side == 0);
temp2 = find(chosen_side == 1);
likelihood(temp1) = 1-temp_p(temp1); %left choice
likelihood(temp2) = temp_p(temp2);

if length(temp1)+length(temp2) ~= N_trial
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Choice_trial,tone_evidence,trial_evidence,use_trial_remove_first,...
    low,high,correct,error,flip_tone,number_use_trial,number_use_trial_remove_first,...
    binary_tone,right_trial_all,number_trial_all,right_trial,number_trial] ...
 = HMM_get_basic_task_structure_20210514(filename1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


load(filename1)


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
for i = 1:length(temp_evi) %
    temp = find(EvidenceStrength == temp_evi(i));% 
    temp_left  = intersect(temp,low); %
    temp_right = intersect(temp,high);%
    trial_evidence(temp_left)  = temp_evi_low(i);%
    trial_evidence(temp_right) = temp_evi_high(i); % 
end

%Make the true tone cloud value
binary_tone = nan(length(Tone_cloud),1);
for i = 1:length(Tone_cloud)
    temp_tone = Tone_cloud(i).matrix;
    %Get the data in all sound
    temp1 = find(temp_tone >= 9);
    binary_tone(i,1) = length(temp1) ./ length(temp_tone);  
end

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
    % [binary_tone, new_tone_evi];
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

number_trial = nan(length(tone_evidence));
right_trial = nan(length(tone_evidence));
for i = 1:length(tone_evidence)
    temp_trial = find(use_tone == tone_evidence(i));
    temp_choice = use_choice(temp_trial);
    number_trial(i) = length(temp_trial);
    right_trial(i) = sum(temp_choice);
end

return

