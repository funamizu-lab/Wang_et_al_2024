function Figure1F_G_H_nBackTrial_learningphase_20240420_aki


root_dir = '/Volumes/Extreme SSD/';

path_file_repeat = repeat_learning_20240420(root_dir);
path_file_alternate = alternate_learning_20240420(root_dir);


[sabun_likeli_rep,ave_sabun_likeli_rep,use_para_rep,mean_para_correct_rep,mean_para_error_rep]...
    = getData(path_file_repeat);

[sabun_likeli_alt,ave_sabun_likeli_alt,use_para_alt,mean_para_correct_alt,mean_para_error_alt]...
    = getData(path_file_alternate);


plot_para_correct_error(sabun_likeli_rep,ave_sabun_likeli_rep,use_para_rep, mean_para_correct_rep,mean_para_error_rep,...
    sabun_likeli_alt,ave_sabun_likeli_alt,use_para_alt,mean_para_correct_alt,mean_para_error_alt)



function [sabun_likeli,ave_sabun_likeli,use_para, mean_para_correct,mean_para_error] = getData(path_file)

folders = path_file;
all_directory = eval('folders');
use_session = 1:7;
n_back = 7;

BIC_all = [];
for i = 1:length(all_directory)
    disp(['mouse' string(i)])
    cd(all_directory{i}); %mouse directory
    
    use_directory = dir;
    clear dir_name
    count = 0;
    for j = 1:length(use_directory)
        temp_dir = use_directory(j).name;
        if length(temp_dir) > 50 
            count = count + 1;
            dir_name{count} = temp_dir;
        end
    end
    clear BIC loglikeli
    for j = 1:length(use_session)
        cd(dir_name{use_session(j)})
        bpod_file = dir('Bpod*.mat');
        if length(bpod_file) ~= 1
            hoge
        else
            [xcorrect,xerror,Choice,TrialEvidence] = getChoiceData(bpod_file.name);
            [BIC(j,:), loglikeli(j,:), Coef_rep(i,j).matrix,Pval_rep(i,j).matrix]= nBackRegression_aki(Choice,TrialEvidence,xcorrect,xerror,n_back);
        end
        
        cd ../
    end
    
    %check BIC
    BIC_all = [BIC_all; BIC];
    ave_BIC(i,:) = mean(BIC);
    ave_likeli(i,:) = mean(loglikeli);
end
mean(ave_BIC)

mean_likeli = mean(ave_likeli);
%compute the likelihood ratio test
p_chi = nan(length(mean_likeli)-1,2);
ave_sabun_likeli = nan(1,length(mean_likeli)-1);
% sabun_likeli = nan(11,length(mean_likeli)-1);

for i = 1:length(mean_likeli)-1
    temp1 = mean_likeli(i+1) - mean_likeli(i);
    temp1 = 2 * temp1;
    p_chi(i,1) = chi2cdf(temp1,2,'upper');
    %p_chi(i) = chi2cdf(temp1,1,'upper');
    
    [~,p_chi(i,2)] = lratiotest(mean_likeli(i+1),mean_likeli(i),2);
    
    %This is for plotting
    sabun_likeli(:,i) = ave_likeli(:,i+1) - ave_likeli(:,i);
    ave_sabun_likeli(i) = mean_likeli(i+1) - mean_likeli(i);
end
% p_chi
p_chi(1,1)
p_chi(2,1)

use_para = 5;
session_all = [];
for i = 1:length(all_directory)
    temp = length(use_session);
    session_all = [session_all; ones(temp,1)*i];
end

%Analyze the regression coefficient
%Coef(i+1).matrix
all_para_sound = [];
all_para_correct = [];
all_para_error = [];
for i = 1:length(all_directory)
    clear para_sound para_correct para_error
    for j = 1:length(use_session)
        temp = Coef_rep(i,j).matrix;
        temp_para = temp(use_para+1).matrix';
        para_sound(j,1) = temp_para(1);
        para_correct(j,:) = temp_para(2:1+use_para);
        para_error(j,:) = temp_para(2+use_para:1+use_para*2);
    end
    all_para_sound = [all_para_sound; para_sound];
    all_para_correct = [all_para_correct; para_correct];
    all_para_error = [all_para_error; para_error];
    
    mean_para_sound(i,:) = mean(para_sound);
    mean_para_correct(i,:) = mean(para_correct);
    mean_para_error(i,:) = mean(para_error);
end

p_lme_correct = nan(use_para,2);
p_lme_error = nan(use_para,2);
p_correct = nan(i);
p_error = nan(i);


data_name{1} = 'value';
data_name{2} = 'random1';
for i = 1:use_para
    p_correct(i) = signrank(mean_para_correct(:,i));
    p_error(i) = signrank(mean_para_error(:,i));
    
    %linear mixed effect model
    data = [all_para_correct(:,i),session_all];
    tbl = table(data(:,1),data(:,2),'VariableNames',data_name);
    lme = fitlme(tbl,'value ~ 1+(1|random1)');
    p_lme_correct(i,1) = lme.Coefficients.pValue;
    p_lme_correct(i,2) = signrank(all_para_correct(:,i));
    
    data = [all_para_error(:,i),session_all];
    tbl = table(data(:,1),data(:,2),'VariableNames',data_name);
    lme = fitlme(tbl,'value ~ 1+(1|random1)');
    p_lme_error(i,1) = lme.Coefficients.pValue;
    p_lme_error(i,2) = signrank(all_para_error(:,i));
end
% p_correct
% p_error
signrank(abs(mean_para_correct(:,1)),abs(mean_para_correct(:,2)))
signrank(abs(mean_para_error(:,1)),abs(mean_para_error(:,2)))

% p_lme_correct
% p_lme_error

signrank(abs(mean_para_correct(:,1)),abs(mean_para_error(:,1)))




function plot_para_correct_error(sabun_likeli_rep,ave_sabun_likeli_rep,use_para_rep, mean_para_correct_rep,mean_para_error_rep,...
    sabun_likeli_alt,ave_sabun_likeli_alt,use_para_alt,mean_para_correct_alt,mean_para_error_alt)

figure
plot(sabun_likeli_rep(:,1:use_para_rep)','k')
hold on
plot(ave_sabun_likeli_rep(1:use_para_rep),'r')
set(gca,'xlim',[0.5 use_para_rep+0.5])
hold on
yline(0,":k")
xlabel('Number of back trials')
ylabel('loglikelihood')
title('Repeat (8 mice)')

figure
plot(sabun_likeli_alt(:,1:use_para_alt)','k')
hold on
plot(ave_sabun_likeli_alt(1:use_para_alt),'r')
set(gca,'xlim',[0.5 use_para_rep+0.5])
hold on
yline(0,":k")
xlabel('Number of back trials')
ylabel('loglikelihood')
title('Alternate (11 mice)')


figure
plot_mean_se_moto(mean_para_correct_rep,[1 0 0],2)
hold on
plot_mean_se_moto(mean_para_correct_alt,[0 0 1],2)
hold on
yline(0,":k")
set(gca,'xlim',[0.5 use_para_rep+0.5])
xlabel('Number of back trials')
ylabel('regression coefficient')
title('Correct trials')


figure
plot_mean_se_moto(mean_para_error_rep,[1 0 0],2)
hold on
plot_mean_se_moto(mean_para_error_alt,[0 0 1],2)
hold on
yline(0,":k")
set(gca,'xlim',[0.5 use_para_alt+0.5])
xlabel('Number of back trials')
ylabel('regression coefficient')
title('Error trials')


rand_x_rep = (rand(length(mean_para_correct_rep(:,1)),1)-0.5) * 0.2;
figure
subplot(1,2,1)
boxplot([abs(mean_para_correct_rep(:,1)),abs(mean_para_correct_rep(:,2))])
hold on
plot(1+rand_x_rep,abs(mean_para_correct_rep(:,1)),'k.')
hold on
plot(2+rand_x_rep,abs(mean_para_correct_rep(:,2)),'k.')
title('correct coefficient')
subplot(1,2,2)
boxplot([abs(mean_para_error_rep(:,1)),abs(mean_para_error_rep(:,2))])
hold on
plot(1+rand_x_rep,abs(mean_para_error_rep(:,1)),'k.')
hold on
plot(2+rand_x_rep,abs(mean_para_error_rep(:,2)),'k.')
title('error coefficient')
sgtitle('Repeating condition')


rand_x_alt = (rand(length(mean_para_correct_alt(:,1)),1)-0.5) * 0.2;
figure
subplot(1,2,1)
boxplot([abs(mean_para_correct_alt(:,1)),abs(mean_para_correct_alt(:,2))])
hold on
plot(1+rand_x_alt,abs(mean_para_correct_alt(:,1)),'k.')
hold on
plot(2+rand_x_alt,abs(mean_para_correct_alt(:,2)),'k.')
title('correct coefficient')
subplot(1,2,2)
boxplot([abs(mean_para_error_alt(:,1)),abs(mean_para_error_alt(:,2))])
hold on
plot(1+rand_x_alt,abs(mean_para_error_alt(:,1)),'k.')
hold on
plot(2+rand_x_alt,abs(mean_para_error_alt(:,2)),'k.')
title('error coefficient')
sgtitle('Alternating condition')



return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xcorrect,xerror,Choice,TrialEvidence] = getChoiceData(filename1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~,~,trial_evidence,use_trial_remove_first,...
    ~,~,correct,error,~,~,~,...
    ~,~,~,~,~] ...
    = HMM_get_basic_task_structure_20210514(filename1);
load(filename1);
Choice = Chosen_side(use_trial_remove_first,:);
TrialEvidence = trial_evidence(use_trial_remove_first,:);
xcorrect = ones(length(Chosen_side),1);
xerror =  ones(length(Chosen_side),1);

left = find(Chosen_side == 0);

left = intersect(left, use_trial_remove_first);
correct = intersect(correct, use_trial_remove_first);
error = intersect(error, use_trial_remove_first);

xcorrect(left) = -1;
xcorrect(error) = 0;
xerror(left) = -1;
xerror(correct) = 0;

xcorrect = xcorrect(use_trial_remove_first,:);
xerror = xerror(use_trial_remove_first,:);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [BIC, loglikeli, Coef,Pval]= nBackRegression_aki(Choice,TrialEvidence,xcorrect,xerror,N)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PredicName=cell(1,2*N+2);
CorrectSide = ones(length(Choice),1);
for t = 1 : length(TrialEvidence)
    if TrialEvidence(t)<0.5
        CorrectSide(t) = 0;
    end
end

% predictor 
CorrectPred = zeros(length(xcorrect)-N,N);
ErrorPred = zeros(length(xerror)-N,N);
for i=1:N
    CorrectPred(:,i) = xcorrect(N+1-i:end-i);
    ErrorPred(:,i) = xerror(N+1-i:end-i);
    PredicName{i+1}=['correct',num2str(i)];
    PredicName{i+1+N}=['error',num2str(i)];
end
RightTone = TrialEvidence(N+1:end);


% data y (explain) %

y = Choice(N+1:end);
choice0 = find(y == 0);

%start with 0 back
X=RightTone;
[~,size_para] = size(X);
mdl = fitglm(X,y,'Distribution','Binomial','Link','logit');
Coef(1).matrix = mdl.Coefficients.Estimate(2:end);
Pval(1).matrix = mdl.Coefficients.pValue(2:end);

likeli = mdl.Fitted.Probability;
likeli(choice0) = 1-likeli(choice0);
sum_loglikeli = sum(log(likeli));
%mdl.LogLikelihood
    
%Compute BIC
clear BIC loglikeli
loglikeli(1) = sum_loglikeli;
BIC(1) = -2*sum_loglikeli + (size_para+1) * log(length(y));


for i = 1:N
    X=[RightTone,CorrectPred(:,1:i),ErrorPred(:,1:i)];
    [~,size_para] = size(X);

    % linear regression %
    mdl = fitglm(X,y,'Distribution','Binomial','Link','logit');
    Coef(i+1).matrix = mdl.Coefficients.Estimate(2:end);
    Pval(i+1).matrix = mdl.Coefficients.pValue(2:end);

    likeli = mdl.Fitted.Probability;
    likeli(choice0) = 1-likeli(choice0);
    sum_loglikeli = sum(log(likeli));
    %mdl.LogLikelihood
    
    %Compute BIC
    loglikeli(i+1) = sum_loglikeli;
    BIC(i+1) = -2*sum_loglikeli + (size_para+1) * log(length(y));
end

return


function analysis_dir = repeat_learning_20240420(root_dir)

analysis_dir = {
    [root_dir '/Zigzag_repeat_Ephys/Behavior_during_learning/Repeat_first/a16/repeat']
    [root_dir '/Zigzag_repeat_Ephys/Behavior_during_learning/Repeat_first/W01/repeat']
    [root_dir '/Zigzag_repeat_Ephys/Behavior_during_learning/Repeat_first/W02/repeat']
    [root_dir '/Zigzag_repeat_Ephys/Behavior_during_learning/Repeat_first/W08/repeat']
    [root_dir '/Zigzag_repeat_Ephys/Behavior_during_learning/Repeat_first/W13/repeat']
    [root_dir '/Zigzag_repeat_Ephys/Behavior_during_learning/Repeat_first/W14/repeat']
    [root_dir '/Zigzag_repeat_Ephys/Behavior_during_learning/Repeat_first/W22/repeat']
    [root_dir '/Zigzag_repeat_Ephys/Behavior_during_learning/Repeat_first/W31/repeat']
    };
return

%}
function analysis_dir = alternate_learning_20240420(root_dir)

analysis_dir = {
    [root_dir '/Zigzag_repeat_Ephys/Behavior_during_learning/Zigzag_first/a17/zigzag']
    [root_dir '/Zigzag_repeat_Ephys/Behavior_during_learning/Zigzag_first/W03/zigzag']
    [root_dir '/Zigzag_repeat_Ephys/Behavior_during_learning/Zigzag_first/W04/zigzag']
    [root_dir '/Zigzag_repeat_Ephys/Behavior_during_learning/Zigzag_first/W05/zigzag']
    [root_dir '/Zigzag_repeat_Ephys/Behavior_during_learning/Zigzag_first/W06/zigzag']
    [root_dir '/Zigzag_repeat_Ephys/Behavior_during_learning/Zigzag_first/W21/zigzag']
    [root_dir '/Zigzag_repeat_Ephys/Behavior_during_learning/Zigzag_first/W23/zigzag']
    [root_dir '/Zigzag_repeat_Ephys/Behavior_during_learning/Zigzag_first/W26/zigzag']
    [root_dir '/Zigzag_repeat_Ephys/Behavior_during_learning/Zigzag_first/W27/zigzag']
    [root_dir '/Zigzag_repeat_Ephys/Behavior_during_learning/Zigzag_first/W29/zigzag']
    [root_dir '/Zigzag_repeat_Ephys/Behavior_during_learning/Zigzag_first/W30/zigzag']
    };
