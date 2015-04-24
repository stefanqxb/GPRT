function [ e ] = detect( e,substring)
    if isempty(e)
       e = [e,substring];
    else
       m = strcmp(e,substring);
       t = find(m==1);
       if size(t,2)==0
          e = [e,substring];
       end           
    end

end

