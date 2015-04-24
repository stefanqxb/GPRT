function [ final_score ] = score_sys(a,word)
    % a: peptide
    % word: compare subsequence
    % score : local alignment
    refer = 'ARNDCQEGHILKMFPSTWYV'; 
    score_mat = blosum(50,'EXTENDED',false);  % A R N D C Q E G H I L K M F P S T W Y V 
    word = char(word);
    len = size(a,2);
    len2 = size(word,2);
    temp = zeros(len2+1,len+1);
    d = 1; % gap open penalty 
    
    for i =2:len2+1
        for j =2:len+1
            ind_x = strfind(refer,word(i-1));   
            ind_y = strfind(refer,a(j-1));   
            s_value = score_mat(ind_x,ind_y);
            vector = [temp(i-1,j-1)+s_value,temp(i-1,j)-d,temp(i,j-1)-d,0];
            score = max(vector);
            temp(i,j) = score;
        end
    end
    
    final_score = max(max(temp));

end

