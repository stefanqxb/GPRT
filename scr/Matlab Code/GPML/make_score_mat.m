function [ score_mat ] = make_score_mat(aa_index)
      % input percent : percentage of each amoid acid 
      %       peptide : peptide
      % output score_mat: score_mat similar to blosum
      
      p = 0.04;
      scoremat_1 = repmat(aa_index,[1,20]);
      scoremat_2 = scoremat_1';
      score_mat = exp(real((scoremat_1 -scoremat_2).^p));

end

