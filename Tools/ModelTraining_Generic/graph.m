function [outputArg1,outputArg2] = graph(inputArg1,inputArg2)
%GRAPH Summary of this function goes here
%   Detailed explanation goes here


T = load('C:\Users\BL\Documents\Data_synthetic_bias\TrueInspector.mat');
D = load('C:\Users\BL\Documents\Data_synthetic_bias\Generated Data.mat');
D.GeneratedData.TrueState(MdataEngy.RemovedData(:,1))=[];
D.GeneratedData.ShortObserved(MdataEngy.RemovedData(:,1))=[];



T = load('C:\Users\BL\Documents\Synthetique_Bias_Sept\TrueInspector.mat');
%graph
PARAM
% %histogram of estimated parameters
figure;
subplot(1,2,1);
histogram(UpdatedInspectorsData{1}(:,3));
xlabel('Estimated \sigma_V ^2(I_i)');
subplot(1,2,2);
hold on
histogram(UpdatedInspectorsData{1}(:,2));
xlabel('Estimated \mu_V (I_i)');
xline(mean(UpdatedInspectorsData{1}(:,2)));
hold off
%
%
% %comparaison mu/sigma for each inspector
figure;
scatter(mu_v,mu_v2hat,'filled');
xlabel('Estimated \mu_V (I_i)');
ylabel('Estimated \sigma_V^2(I_i)');
%
%% %comparaison mu/sigma for each inspector
figure;
scatter(mu_v,mu_v2hat,'filled');
xlabel('Estimated \mu_V (I_i)');
ylabel('Estimated \sigma_V^2(I_i)');


T = load('C:\Users\BL\Documents\Synthetique_Bias_Sept\TrueInspector.mat');
graph
       % comparaison ecart

         T = load('C:\Users\BL\Documents\Synthetique_Bias_Nov_alt_Q\TrueInspector.mat');
%         
%         ecart_sigma = sum((T.TrueInsp(:,1) - UpdatedInspectorsData{1}(:,3)).^2)
%         
%         if OperationIndex>=3
%             ecart_mu = sum((T.TrueInsp(:,2) - UpdatedInspectorsData{1}(:,2)).^2)   
%         end
%         
        figure;
        subplot(1,2,1);
        %comparaison to the true values \sigm_V  for each inspecteur
        scatter(T.TrueInsp(:,1),UpdatedInspectorsData{1}(:,3),'o');
        hold on
        plot([0 36],[0,36],'--');
        xlabel('\fontsize{20}True \sigma_V (I_i)');
        ylabel('\fontsize{20}Estimated \sigma_V  (I_i)');
        hold off
        axis equal;
        ylim([1 7]);
        xlim([1 7]);
        legend('\fontsize{20}Online \sigma_{V}');
        title('\fontsize{20}Standard deviation estimation');
        %comparaison to the true values \mu_V  for each inspecteur
        subplot(1,2,2);
        scatter(T.TrueInsp(:,2),UpdatedInspectorsData{1}(:,2),'o');
        hold on
        plot([-20 20],[-20,20],'--');
        xlabel('\fontsize{20}True \mu_V (I_i)');
        ylabel('\fontsize{20}Estimated \mu_V  (I_i)');
        hold off
        axis equal;
        ylim([-5 5]);
        xlim([-5 5]);
        title('\fontsize{20}Bias estimation');
%scatter plot
figure;
subplot(1,2,1);
%comparaison to the true values \sigm_V  for each inspecteur
scatter(T.TrueInsp(:,1),sqrt(variance)','o');
hold on
plot([0 36],[0,36],'--');
xlabel('True \sigma_V (I_i)');
ylabel('Estimated \sigma_V  (I_i)');
hold off
axis equal;
ylim([1 7]);
xlim([1 7]);
legend('Online \sigma_{V}');
title('Standard deviation estimation, initial [' + string(sd_v0) + ',' + string(sd_sd_v0) + ']');
%comparaison to the true values \mu_V  for each inspecteur
subplot(1,2,2);
scatter(T.TrueInsp(:,2),mu_v','o');
hold on
plot([-20 20],[-20,20],'--');
xlabel('True \mu_V (I_i)');
ylabel('Estimated \mu_V  (I_i)');
hold off
axis equal;
ylim([-5 5]);
xlim([-5 5]);
title('Bias estimation, initial [' + string(mu_v0) + ',' + string(mu_sd_V0) + ']');

%
%
%
%
% %evolution of the estimation for inspector ii
ii = 86;
figure;
Mu_Mu = unique(save_mu_mu(:,ii)',"stable");
Mu_Sig = unique(save_mu_sig(:,ii)',"stable");
Var_Mu = unique(save_mu_mu(:,ii)',"stable");
Var_Sig = unique(save_var_sig(:,ii)',"stable");

nbo = count(ii);
x = linspace(1,nbo,nbo);
subplot(1,2,1);
hold on
plot(x, save_mu_mu(1:nbo,ii)','b');
plot(x,  T.TrueInsp(ii,2)*ones(size(x)),'--');
s1 = patch([x fliplr(x)], [(save_mu_mu(1:nbo,ii) + sqrt(save_mu_sig(1:nbo,ii)))' fliplr((save_mu_mu(1:nbo,ii) - sqrt(save_mu_sig(1:nbo,ii)))')], 'b');
s2 = patch([x fliplr(x)], [(save_mu_mu(1:nbo,ii) + 2*sqrt(save_mu_sig(1:nbo,ii)))' fliplr((save_mu_mu(1:nbo,ii) - 2*sqrt(save_mu_sig(1:nbo,ii)))')], 'b');
hold off
alpha(s1,.3);
alpha(s2,.2);
xlabel('number of observations');
ylabel('\mu_v');
legend('estimated \mu_v','true \mu_v','\mu_v \pm \sigma_{\muV}','\mu_v \pm 2*\sigma_{\muV}');
ylim([-5 5]);
xlim([0 nbo]);
subplot(1,2,2);
hold on
plot(x, save_var_mu(1:nbo,ii)','b');
plot(x, T.TrueInsp(ii,1).^2*ones(size(x)),'--');
s1 = patch([x fliplr(x)], [(save_var_mu(1:nbo,ii) +  save_var_sig(1:nbo,ii))' fliplr((save_var_mu(1:nbo,ii) - save_var_sig(1:nbo,ii))')], 'b');
s2 = patch([x fliplr(x)], [(save_var_mu(1:nbo,ii) + 2*save_var_sig(1:nbo,ii))' fliplr((save_var_mu(1:nbo,ii) - 2*save_var_sig(1:nbo,ii))')], 'b');
hold off
alpha(s1,.3);
alpha(s2,.2);
xlabel('number of observations');
ylabel('\sigma_v');
legend('estimated \sigma_v^2','true \sigma_v^2','\sigma_v^2 \pm \sigma_{\sigmaV}','\sigma_v^2 \pm 2*\sigma_{\sigmaV}');
ylim([0 100]);
xlim([0 nbo]);
suptitle('estimation for I' + string(99000 + ii) + ', number of observation : ' + string(count(ii)));

%distribution error
ii = 86;
Y_real = Data_inspect;
[row, column, dataValues] = find(Inspectorlabel == 99000 + ii);
[rowSrt, ind] = sort(row);
columnSrt = column(ind);
est_mu = save_mu_mu(end,ii);
est_sig = sqrt(save_var_mu(end,ii));
x = -45:0.5:45;
figure;
hold on
plot(x,normpdf(x,est_mu,est_sig), 'b','Linewidth',2);
plot(x,normpdf(x,0,T.TrueInsp(ii,1)), 'r','Linewidth',2);
b1 = NaN(size(rowSrt,1),1) ;
b2 = NaN(size(rowSrt,1),1) ;
for a = 1:size(rowSrt,1)
    if ~ isnan(Y_real(rowSrt(a),columnSrt(a)))
        Obe1 = Y_real(rowSrt(a),columnSrt(a));
        Tru = D.GeneratedData.TrueState{rowSrt(a), 1}(1,D.GeneratedData.ShortObserved{rowSrt(a), 1}(1,3)-1964+columnSrt(a));
        Diff1 = Obe1 - Tru;
        b1(a,1) = Diff1;
        Obe2 = test(rowSrt(a),columnSrt(a));
        Diff2 = Obe2 - Tru;
        b2(a,1) = Diff2;
    end
end
histogram(b1,15,'Normalization','pdf');
histogram(b2,15,'Normalization','pdf');
xline(0,'--r');
xline(est_mu,'--b');
hold off
legend('Estimated distribution', 'True distribution', 'y(t) - x(t)', '\mu(t) - x(t)');
title(sprintf('Distribution of observation errors for Inspector ' + string(99000 + ii)));
xlim([-15 15]);
%

T = load('C:\Users\BL\Documents\Data_synthetic_bias\TrueInspector.mat');
D = load('C:\Users\BL\Documents\Synthetique_Bias\Generated Data.mat');
D.GeneratedData.TrueState(MdataEngy.RemovedData(:,1))=[];
D.GeneratedData.ShortObserved(MdataEngy.RemovedData(:,1))=[];

% %distribution best worst
figure;
[~, min_sig] = min(abs(T.TrueInsp(:,1) - UpdatedInspectorsData{1}(:,3)));
[~,min_bias] = min(abs(T.TrueInsp(:,2) - UpdatedInspectorsData{1}(:,2)));
[~, max_sig] = max(abs(T.TrueInsp(:,1) - UpdatedInspectorsData{1}(:,3)));
[~,max_bias] = max(abs(T.TrueInsp(:,2) - UpdatedInspectorsData{1}(:,2)));
l = [min_sig min_bias max_sig max_bias];
m = ["Best", "Best", "Worst", "Worst"];
p = ["standard deviation", "bias", "standard deviation", "bias"];
for c = 1:4
    ii = l(c);
    [row, column, ~] = find(Inspectorlabel == 99000 + ii);
    [rowSrt, ind] = sort(row);
    columnSrt = column(ind);
    est_mu = UpdatedInspectorsData{1}(ii,2);
    est_sig = UpdatedInspectorsData{1}(ii,3);
    x = -45:0.5:45;
    subplot(2,2,c)
    hold on
    plot(x,normpdf(x,est_mu,est_sig), 'b','Linewidth',2);
    plot(x,normpdf(x,T.TrueInsp(ii,2),T.TrueInsp(ii,1)), 'r','Linewidth',2);
    b1 = NaN(size(rowSrt,1),1) ;
    b2 = NaN(size(rowSrt,1),1) ;
    for a = 1:size(rowSrt,1)
        if ~ isnan(Data_inspect(rowSrt(a),columnSrt(a)))
            Obe1 = Data_inspect(rowSrt(a),columnSrt(a));
            Tru = D.GeneratedData.TrueState{rowSrt(a), 1}(1,D.GeneratedData.ShortObserved{rowSrt(a), 1}(1,3)-1964+columnSrt(a));
            Diff1 = Obe1 - Tru;
            b1(a,1) = Diff1;
            Obe2 = test(rowSrt(a),columnSrt(a));
            Diff2 = Obe2 - Tru;
            b2(a,1) = Diff2;
        end
    end
    histogram(b1,30,'Normalization','pdf','FaceColor','#EDB120');
    histogram(b2,30,'Normalization','pdf','FaceColor','#7E2F8E');
    xline(T.TrueInsp(ii,2),'--r');
    xline(est_mu,'--b');
    hold off
    xlim([-15 15]);
    legend('Estimated distribution', 'True distribution', 'y(t) - x(t)', '\mu(t) - x(t)');
    legend('Estimated distribution', 'True distribution', '\mu(t) - x(t)');
    title(sprintf('Distribution of observation errors for Inspector ' + string(99000 + ii)+  ',\n' + string(m(c)) +' estimation of the ' + p(c)));
    %%%%%
end
%



figure;
subplot(1,2,2);
plot(1:n_e , LL);
xlabel('number of epoch');
title('Loglikelihood');
subplot(1,2,1);
em = zeros(n_e,1);
for i = 1:n_e
    em(i,1) = mean((T.TrueInsp(:,1)-sqrt(save_epoch(i,:)')).^2);
end
plot(1:n_e, em);
xlabel('number of epoch');
title('Mean of the square difference between true and estimated standard deviation');

figure;
hold on
scatter(T.TrueInsp(:,1),sqrt(save_epoch(1,:))','ob');
scatter(T.TrueInsp(:,1),sqrt(save_epoch(5,:))','og');
scatter(T.TrueInsp(:,1),sqrt(save_epoch(10,:))','om');
sz = 25;
scatter(T.TrueInsp(:,1),UpdatedInspectorsData{1, 1}(:,3), sz, 'r')
plot([0 6],[0,6],'--');
hold off
xlabel('True \sigma_V (I_i)');
ylabel('Estimated \sigma_V  (I_i)');
legend('1 epoch', '5 epoch', '10 epoch', 'NR')




for i = 1:n_e
    figure;
    scatter(T.TrueInsp(:,1),sqrt(save_epoch(i,:))','o');
    hold on
    scatter(T.TrueInsp(:,1),UpdatedInspectorsData{1, 1}(:,3),'o')
    plot([0 36],[0,36],'--');
    xlabel('True \sigma_V (I_i)');
    ylabel('Estimated \sigma_V  (I_i)');
    hold off
    axis equal;
    ylim([1 7]);
    xlim([1 7]);
    legend('Online \sigma_{V}', 'NR \sigma_{V}');
    title('Standard deviation estimation, initial [' + string(sd_v0) + ',' + string(sd_sd_v0) + ']' + 'epochs = ' + string(i));
end
end

