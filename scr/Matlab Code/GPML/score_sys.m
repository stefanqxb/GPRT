function [ final_score ] = score_sys(a,word,score_mat,aa_index)
    % a: peptide
    % word: compare subsequence
    % score : local alignment
    refer = 'ACDEFGHIKLMNPQRSTVWY';
    word = char(word);
    len = size(a,2);
    len2 = size(word,2);
    temp = zeros(len2+1,len+1); 
    d = 1;
    
    for i =2:len2+1
        for j =2:len+1
            ind_x = strfind(refer,word(i-1));         
            ind_y = strfind(refer,a(j-1));
%             if a(j-1) == '-'
%                if j-2> 0 
%                   ind_b = strfind(refer,a(j-2));
%                else
%                   ind_b = 0;
%                end
%                if i-2 <= len2+1
%                   ind_a = strfind(refer,a(i));
%                else
%                   ind_a = 0;
%                end
%                if ind_b == 0 && ind_a ~= 0
%                    d = abs(d-aa_index(ind_a));
%                else if ind_b ~= 0 && ind_a == 0
%                        d = abs(aa_index(ind_b)-d);
%                     else
%                        d = abs(aa_index(ind_b)-d)+abs(d-aa_index(ind_a));
%                    end
%                end
%             else
%                d = 0;   
%             end
            s_value = score_mat(ind_x,ind_y);
            if size(s_value,2) == 0
                s_value = 0;
            end
            vector = [temp(i-1,j-1)+s_value,temp(i-1,j)-d,temp(i,j-1)-d,0];
            score = max(vector);
            temp(i,j) = score;
        end
    end
    
    final_score = max(max(temp));

end

