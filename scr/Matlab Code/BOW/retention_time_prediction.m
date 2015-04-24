close all
clear all

addpath('~/libsvm-3.20/matlab');

tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculating features
%[retention_time,peptide,orginal] = xlsread('retention_time_peptide.xlsx');
%[data,peptide_2,orginal_2] = xlsread('hydrophobicity.xlsx');
%length_of_peptide = data(:,1);
%hydrophobicity = data(:,2);


%[row,~] = size(peptide);


%counting_mat = extract_aal(peptide); % A R N D C Q E G H I L K M F P S T W Y V O U 
%[idx_x,idx_y] = find(counting_mat>0);
%idx_list = [idx_x,idx_y];
%feature_mat = forming_feature(peptide,idx_list); 

% method 2 using 3 words as feature

% counting_mat = extract_aal_3(peptide); % A R N D C Q E G H I L K M F P S T W Y V O U 
% idx = find(counting_mat>80);
% feature_mat = forming_feature_3(peptide,idx); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% lumi feature benchmark

[retentiontime, peptide, orignal] = xlsread('lumi_fea_final.xlsx');
retention_time = retentiontime(:,1);
feature_mat = retentiontime(:,2:end);

[row,~] = size(peptide);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%weight1 = 0.4;
%weight2 = 0.1;
%weight3 = 0.5;
%feature_mat = [weight1.*feature_mat,weight2.*hydrophobicity,weight3.*length_of_peptide];

for i= 1:size(feature_mat,1)
    feature_mat(i,:) = feature_mat(i,:)/sum(feature_mat(i,:)); 
end

% set up train test set using bow as features
ratio = 0.8;
train_row = round(ratio*row);

train_set = feature_mat(1:train_row,:);
train_targ = retention_time(1:train_row,:);
test_set = feature_mat(train_row+1:end,:);
test_targ = retention_time(train_row+1:end,:);

% SVR
%disp('building the model.....');
%model = svmtrain(train_targ,train_set,'-s 3 -t 2 -g 40 -p 0.01 -c 100 -h 0 '); 
%disp('predicting the label....');
%[predicted_label, accuracy, decision_values] = svmpredict(test_targ, test_set, model);



%RVM

%kernel_= 'gauss';
%width = 0.2;
%maxIts = 500;

%disp('running the mvrvm');

%PHI     = sbl_kernelFunction(train_set,train_set,kernel_,width);
%[weights, used, alpha, beta] = mvrvm(PHI,train_targ,maxIts);

%PHI     = sbl_kernelFunction(test_set,train_set(used,:),kernel_,width);
%predicted_label = PHI*weights;


%%%%%%%%%GPML

x = train_set;
y = train_targ;
z = test_set;


fitgen = zeros(10,1);
timewin = zeros(10,1);
i = 1;


likfunc = @likGauss; sn = 0.1; hyp.lik = log(sn);
covfunc = @covSEiso; hyp2.cov = [0.5 ; 0.5]; hyp2.lik = log(1);


for t = 10:10:100
disp('minize the hyperparmeter');
hyp2 = minimize(hyp2, @gp, -t, @infExact, [], covfunc, likfunc, x, y);
%[hyp2 , fitgen, timewin, counter_time, counter_pre] = PSO_hypa(hyp2, @gp, test_targ, @infExact, [], covfunc, likfunc, x, y, z);
disp('predicting the label');
[predicted_label , s2] = gp(hyp2, @infExact, [], covfunc, likfunc, x, y, z);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

diff_mat = abs(predicted_label - test_targ);
max_t = max(test_targ);
min_t = min(test_targ);


%figure;
%hold on;
%plot(test_targ,predicted_label,'r*');
%plot([0,max_t],[0,max_t]);
%xlabel('observed retention time');
%ylabel('predicted retention time');
%title('BOW predictor of retention time');

%figure(2);
%hold on;
%plot(test_targ,'r*');
%plot(predicted_label,'b.');
%legend('real value','predicted value','Location','NorthEast');

step=100;


%figure(3);
%hist(diff_mat,step);
my_hist = hist(diff_mat,step);
time_interval = time_95_diff(my_hist,max_t,min_t,step);
timewin(i,1) = time_interval/max_t;
gamma = corrcoef(test_targ,predicted_label);

fitgen(i,1) = gamma(1,2);

i = i+1;
end

x =1:size(fitgen,1);
figure(1)
hold on;
grid on;
plot(x,fitgen,'r*-','LineWidth',3);
plot(x,timewin,'bx-','LineWidth',3);
legend('correlation','minimal time window','Location','NorthEast');
title('Result for each generation of PSO');

%figure(2)
%hold on;
%grid on;
%plot(counter_pre','*-','LineWidth',2);
%xlabel('Number of generation');
%ylabel('Predicted Precision of each Particle');

%figure(3)
%hold on;
%grid on;
%plot(counter_time','o-','LineWidth',2);
%xlabel('Number of generation');
%ylabel('Minimal Time window of each Particle');


toc;














