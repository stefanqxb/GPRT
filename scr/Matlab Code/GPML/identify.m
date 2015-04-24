function [ score] = identify(a,data_set,score_mat,aa_index)
  % output : score : compare with 3-grams and the peptide 
  % input : a: one peptide sequence
  %         data_set: subsequence for compare
  
  [~,low_d] = size(data_set);
  score = zeros(1,low_d);
  
   % A R N D C Q E G H I L K M F P S T W Y V 
   
  for i =1: low_d
       score(1,i) = score_sys(a,data_set(i),score_mat,aa_index);      
  end
  
%   threshold = mean(score); 
%   ind = find(score<threshold);
%   score(ind) = 0;
  
end


































