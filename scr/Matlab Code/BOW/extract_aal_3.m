function [ count_mat ] = extract_aal_3(orign)
  [row,~] = size(orign);
  count_mat =zeros(20,20,20);
  
  for i =1: row
      temp = char(orign{i});
      t = strfind(temp,'.');
      peptide = [ temp(1:t(1)-1),temp(t(1)+1:t(2)-1),temp(t(2)+1:end)];
      pep_len = size(peptide,2)-2;  
         for j =1: pep_len
             a = peptide(j:j+2);
             idx =zeros(3,1);
                for k =1:3
                    temp_rec = identify(a(k));
                    if ~isempty(temp_rec) % avoid .-(gap)
                        idx(k) = temp_rec;
                    end
                end
             if idx(1)~=0 && idx(2)~=0 && idx(3)~=0
                count_mat(idx(1),idx(2),idx(3)) = count_mat(idx(1),idx(2),idx(3))+1;
             end
         end  
  end



end

