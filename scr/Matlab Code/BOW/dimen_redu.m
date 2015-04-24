function [ result ] = dimen_redu( aa_bow )
   [~,low] = size(aa_bow);
   m =1:5:low;
   result=[];
   for i =1:size(m,2)
       temp = aa_bow{1,m(i)};
       temp = cellstr(temp);
       result = [result, temp];
   end
    

end

