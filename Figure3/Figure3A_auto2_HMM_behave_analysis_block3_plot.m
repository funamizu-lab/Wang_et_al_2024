

function Figure3A_auto2_HMM_behave_analysis_block3_plot


data_rep = load('repeat_sabun_all_each_session');
ylimit_rep = [0.1, 0.7];
plot_sabun_each_session(data_rep,ylimit_rep)

ylimit_zig = [-0.7, 0.1];
data_zig = load('zigzag_sabun_all_each_session');
plot_sabun_each_session(data_zig,ylimit_zig)

return


function plot_sabun_each_session(data, ylimit)


sabun_all = data.sabun_all;

[size_mouse,~] = size(sabun_all);
temp_x1 = [1;2;3];
temp_x2 = [1;1;1];
temp_x1 = [temp_x1,temp_x2];

mouse_b = nan(size_mouse,2);

for i = 1:size_mouse
    b = regress(sabun_all(i,1:3)',temp_x1);
    mouse_b(i,:) = b;
end


figure

plot_mean_se_moto_x_axis(sabun_all(:,1:3), [1,2,3], [0 0 0],2)
set(gca,'xlim',[0.5 3.5],'ylim',ylimit)

signrank(mouse_b(:,1))
signrank(mouse_b(:,2))

output = getCoeffPval_overtrain(sabun_all);
output.CI
output.Pval

clear sabun_all

return

function output = getCoeffPval_overtrain(Bias)
mouseval=Bias(:,5);

Pval = nan(4,1);
CIs_lower = nan(4,1);
CIs_upper = nan(4,1);

for i = 1:4
    data = Bias(:,i);

    tbl = table(data,mouseval,'VariableNames',{'value','mouse'});
    lme = fitlme(tbl,'value ~ 1 + (1|mouse)'); %linear regression

    estVal= lme.Coefficients.Estimate;
    Pval(i) = lme.Coefficients.pValue;
    CIs_lower(i,:)= lme.Coefficients.Lower;
    CIs_upper(i,:)= lme.Coefficients.Upper;

end
CI = [CIs_upper,CIs_lower];
output.estVal= estVal;
output.CI = CI(1:3,:);
output.Pval  = Pval(4);
output.CI
output.Pval
return


