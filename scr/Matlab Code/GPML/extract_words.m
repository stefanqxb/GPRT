function [ data_set ] = extract_words(orign)
    [row,~] = size(orign);
    n = 4;
    temp =[];
    for i =1:row       
        cp = char(orign{i});
        t = strfind(cp,'.');
        cp = [ cp(1:t(1)-1),cp(t(1)+1:t(2)-1),cp(t(2)+1:end)];
        t = strfind(cp,'-');
        if t~=0
           cp = [ cp(1:t(1)-1),cp(t(1)+1:end)];
        end
        low = size(cp,2);
        pep_len = low-n+1;
        sub_temp = cell(1,pep_len);
        for j =1:pep_len
            sub_temp{j} = cp(j:j+n-1);
        end         
        ind = 1:size(sub_temp,2);
        sub_temp = sub_temp(ind);
        temp = [sub_temp,temp];  
    end
    
    data_set_temp = unique(temp);
    [~,nc] = size(data_set_temp);
    account = zeros(1,nc);
    for i =1:nc
        ind = strcmp(temp,data_set_temp(i));
        g = find(ind ~= 0);
        account(i) = size(g,2);
    end
    threshold = 8;
    indx = find(account>=threshold);
    data_set = data_set_temp(indx);
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
end

