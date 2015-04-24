function [ string ] = forming_string( list )
   [x,y] = size(list);
   string = cell(1,x);
   temp = [];
   for i =1:x
       for j=1:y
           chart = lookup_table( list(i,j) );
           temp = [temp,chart];
       end
       string{1,i} = temp;
       temp = [];
   end
   

end

