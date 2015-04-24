close all
clear all

tic;
% calculating features
[retention_time,peptide,orginal] = xlsread('retention_time_peptide.xlsx');
[data,peptide_2,orginal_2] = xlsread('hydrophobicity.xlsx');
length_of_peptide = data(:,1);
hydrophobicity = data(:,2);

% split data reduction

peptide = peptide(1:1000);

[row,~] = size(peptide);

% method 1 using  4 words as feature

data_set = extract_words(peptide);
counting_mat = extract_aal(peptide,data_set); % A R N D C Q E G H I L K M F P S T W Y V O U 

% demension reduction

feature_mat = cell2mat(counting_mat);



xlswrite('chosen_peptide.xlsx',feature_mat);
toc;

