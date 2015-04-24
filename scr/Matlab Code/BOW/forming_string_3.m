function [ string ] = forming_string_3( list )
   [x,~] = size(list);
   string = cell(1,x);
   temp = [];
   for i =1:x
       cor_z = fix(list(i)/400)+1;
       rest = mod(list(i),400);
       if rest == 0
           cor_y = 20;
           cor_x = 20;
       else
           cor_y = fix(rest/20)+1;
           cor_x = mod(rest,20);
           if cor_x ==0
               cor_x = 20;
           end
       end         
       chart = lookup_table(cor_z);
       temp = [temp,chart];
       chart = lookup_table(cor_y);
       temp = [temp,chart];
       chart = lookup_table(cor_x);
       temp = [temp,chart];
       string{1,i} = temp;
       temp = [];
   end


end

