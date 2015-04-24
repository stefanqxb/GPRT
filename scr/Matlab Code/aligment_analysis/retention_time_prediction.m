close all
clear all

addpath('~/libsvm-3.20/matlab');

tic;
% calculating features
[retention_time,peptide,orginal] = xlsread('retention_time_peptide.xlsx');
[data,peptide_2,orginal_2] = xlsread('hydrophobicity.xlsx');
length_of_peptide = data(:,1);
hydrophobicity = data(:,2);

% split data reduction



% method 1 using  3 words as feature

data_set = extract_words(peptide);
counting_mat = extract_aal(peptide,data_set); % A R N D C Q E G H I L K M F P S T W Y V O U 

% demension reduction

feature_mat = cell2mat(counting_mat);

row = size(feature_mat,1);

for i= 1:row
    feature_mat(i,:) = feature_mat(i,:)/sum(feature_mat(i,:)); 
end



% set up train test set using bow as features
ratio = 0.8;
train_row = round(ratio*row);

train_set = feature_mat(1:train_row,:);
train_targ = retention_time(1:train_row,:);
test_set = feature_mat(train_row+1:end,:);
test_targ = retention_time(train_row+1:row,:);

% RVM
%kernel_ = 'gauss';
%width = 0.2;
%maxIts = 500;


%PHI     = sbl_kernelFunction(train_set,train_set,kernel_,width);
%[weights, used, alpha, beta] = mvrvm(PHI,train_targ,maxIts);
%PHI     = sbl_kernelFunction(test_set,train_set(used,:),kernel_,width);
%predicted_label = PHI*weights;

%SVR
%model = svmtrain(train_targ,train_set,'-s 3 -t 2 -g 50 -p 0.1 -c 1000 -h 0 '); 
%[predicted_label, accuracy, decision_values] = svmpredict(test_targ, test_set, model);

%%%%%%%%%%%%%%% GPML
x = train_set;
y = train_targ;
z = test_set;
gamma_t = zeros(10:1);
timewin_t = zeros(10,1);

for t= 10:10:100
likfunc = @likGauss; sn = 0.1; hyp.lik = log(sn);
covfunc = @covSEiso; hyp2.cov = [0; 0]; hyp2.lik = log(0.1);
disp('minize the hyperparmeter');
hyp2 = minimize(hyp2, @gp, -t, @infExact, [], covfunc, likfunc, x, y);
disp('predicting the label');
[predicted_label , s2] = gp(hyp2, @infExact, [], covfunc, likfunc, x, y, z);

%%%%%%%%%%%%%%%

%figure;
%hold on;
%plot(test_targ,predicted_label,'r*');
%plot([0,max_v],[0,max_v],'b','LineWidth',3);
%xlabel('observed retention time');
%ylabel('predicted retention time');
%title('Alignment_predictor of retention time');

gamma = corrcoef(test_targ,predicted_label);
gamma_t(i,1) = gamma(1,2);

diff_mat = abs(predicted_label - test_targ);
max_t = max(test_targ);
min_t = min(test_targ);
step=100;

%figure(2);
%hist(diff_mat,step);
my_hist = hist(diff_mat,step);
time_interval = time_95_diff(my_hist,max_t,min_t,step);
timewin_t(i,1) = time_interval/max_t;

i = i+1;
end

x =1:size(gamma_t,1);
figure(1)
hold on;
grid on;
plot(x,gamma_t,'r*-','LineWidth',3);
plot(x,timewin_t,'bx-','LineWidth',3);
legend('correlation','minimal time window','Location','NorthEast');
title('Result for conjugate gradient of different number of function evaluation');


toc;










