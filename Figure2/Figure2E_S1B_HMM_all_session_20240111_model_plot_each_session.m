%{
----------------------------------------------------------------------------
Analyzing behavioral data
At least for the correct rate
----------------------------------------------------------------------------
%}

function Figure2E_S1B_HMM_all_session_20240111_model_plot_each_session

mouseval_allsession = [1 1 2 3 3 4 4 4 4 5 6 7 7 7 7 8 9 10 11 11 12 13 14 15 15 15 15 15 16 17 17 18 19 19 19];
mouseval_goodsession = [1 2	3 3	4 4	4 4	5 6	6 6	6 7	8 9	10 10 11 12	12 12 12 13	14 14 15 16];

sabun_all = [
    18.29560121	3.804150731	16.5656278	7.246587	2.766478	-0.190671464	0.070073072	0.52569862	11.33123958	1.474455361	12.33872217	18.40549051	10.26928872	23.12186285	-0.245144748	-0.143225688	-0.174629641	16.00827204	7.655632198	-1.200961628	0.546338487	0.145428473	18.05529587	37.68366365	2.831355713	-0.209063233	0.317259961	0.172767192	2.247189142	-0.082955269	103.5435175	4.463242704	38.84353986	10.41078231	3.320376793];
sabun_good = [
    3.804150731	16.5656278	7.246587	2.766478	-0.190671464	0.070073072	0.52569862	11.33123958	12.33872217	18.40549051	10.26928872	23.12186285	-0.245144748	-0.143225688	-0.174629641	16.00827204	7.655632198	-1.200961628	0.145428473	2.831355713	-0.209063233	0.317259961	0.172767192	2.247189142	-0.082955269	103.5435175	4.463242704	3.320376793];

Tokyo3_sabun2(1).value(1,:) = sabun_all;
Tokyo3_sabun2(1).value(2,:) = mouseval_allsession;

Tokyo3_sabun2(2).value(1,:) = sabun_good;
Tokyo3_sabun2(2).value(2,:) = mouseval_goodsession;

p = nan(2,1);p_signrank = nan(2,1);

for i = 1 : 2
    size_x = [zeros(length(Tokyo3_sabun2(i).value),1)];

    plot_x1 = 0.2 * (rand(1,length(Tokyo3_sabun2(i).value)) - 0.5);
    figure('Position', [600 600 300 400])
    boxplot(Tokyo3_sabun2(i).value(1,:),size_x,'Color', [0 0 0])
    hold on
    plot(plot_x1+1, Tokyo3_sabun2(i).value(1,:),'k.')
    hold on
    yline(0,':k')
    ylim([-20,120])
    if i == 1
        title('All Sessions')
    else
        title('Good Sessions')
    end
    
    output = getCoeffPval(Tokyo3_sabun2(i).value);
    p(i) = output.Pval;
    p_signrank(i) = signrank(Tokyo3_sabun2(i).value(1,:));
    % 1 all session 2 good session

end

return


function output = getCoeffPval(Bias)
data = Bias(1,:)';
mouseval=Bias(2,:)';

tbl = table(data,mouseval,'VariableNames',{'value','mouse'});
lme = fitlme(tbl,'value ~ 1 + (1|mouse)'); %linear regression

estVal= lme.Coefficients.Estimate;
Pval = lme.Coefficients.pValue;
CIs_lower= lme.Coefficients.Lower;
CIs_upper= lme.Coefficients.Upper;


CI = [CIs_upper,CIs_lower];
output.estVal= estVal;
output.CI = CI;
output.Pval  = Pval;
output.CI;
output.Pval;
return
