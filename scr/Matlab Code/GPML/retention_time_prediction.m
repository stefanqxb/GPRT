close all
clear all

% addpath('D:\Matlab\R2013a\mvrvm');

tic;
% calculating features
[aa_index,aa] =xlsread('rt_index.xls');
%[aa_index,aa] =xlsread('rt_index2.xls');
%[aa_index,aa] =xlsread('rt_index3.xls');
[retention_time,peptide,orginal] = xlsread('retention_time_peptide.xlsx');
% [data,peptide_2,orginal_2] = xlsread('hydrophobicity.xlsx');
% length_of_peptide = data(:,1);
% hydrophobicity = data(:,2);

% split data reduction

peptide = peptide(1:3000);

% method 1 using  3 words as feature


data_set = extract_words(peptide);
% changable 
% [~,q] = size(data_set);
% data_set = data_set(1:round(q/2));
%
% score_mat1 =make_score_mat(aa_index); % order of aa is different
score_mat2 = blosum(50);
score_mat =  score_mat2;
counting_mat = extract_aal(peptide,data_set,score_mat,aa_index); % A R N D C Q E G H I L K M F P S T W Y V

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

%%%%%%%%%%%%%%% GPML

x = train_set;
y = train_targ;
z = test_set;

likfunc = @likGauss; sn = 0.5; hyp.lik = log(sn);
covfunc = @covSEiso; hyp2.cov = [0; 0]; hyp2.lik = log(0.1);
hyp2 = minimize(hyp2, @gp, -50, @infExact, [], covfunc, likfunc, x, y);
nlml2 = gp(hyp2, @infExact, [], covfunc, likfunc, x, y);
[predicted_label , s2] = gp(hyp2, @infExact, [], covfunc, likfunc, x, y, z);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

max_v = max(test_targ);

figure;
hold on;
plot(test_targ,predicted_label,'r*');
plot([0,max_v],[0,max_v],'b','LineWidth',3);
xlabel('observed retention time');
ylabel('predicted retention time');
title('predictor of retention time');


gamma = corrcoef(test_targ,predicted_label);

diff_mat = abs(predicted_label - test_targ);
max_t = max(test_targ);
min_t = min(test_targ);
step=100;

figure(2);
hist(diff_mat,step);
my_hist = hist(diff_mat,step);
time_interval = time_95_diff(my_hist,max_t,min_t,step);

toc;


















