function [count_mat ] = extract_aal(orign,data_set)
  [row,~] = size(orign);
  count_mat =cell(row,1);
  
  parfor i =1: row
      temp = char(orign{i});
      t = strfind(temp,'.');
      peptide = [ temp(1:t(1)-1),temp(t(1)+1:t(2)-1),temp(t(2)+1:end)];
      temp_mat = identify(peptide,data_set);
      count_mat{i} = sparse(temp_mat);
      i
  end    
  
  
  
end






