function [ feature_mat ] = forming_feature_3(orgin,idx_list)
  [row,~] = size(orgin);
  [row2,~] = size(idx_list);
  feature_mat =zeros(row,row2);
  
  refer = forming_string_3(idx_list);
  
  for i =1: row
      temp = char(orgin{i});
      t = strfind(temp,'.');
      peptide = [ temp(1:t(1)-1),temp(t(1)+1:t(2)-1),temp(t(2)+1:end)];
      pep_len = size(peptide,2)-2;  
         for j =1: pep_len
             a = peptide(j:j+2);
             result = strcmp(refer,a);
             loc = find(result==1);
             if size(loc,2)~=0                
                feature_mat(i,loc) = feature_mat(i,loc)+1;
             end
         end  
  end


end

