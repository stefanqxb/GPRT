function [ feature,counter ] = count_freq( aa_bow )
   
   [row,~] = size(aa_bow);
   e = [];
   for i=1:row
       low2 =size(aa_bow{i,:},2);
       for j =1:low2
           e = detect(e,aa_bow{i}(j));
       end
   end
   
   feature  = e;
   selected_feature = simi_measure( feature );
%    counter = cell(row,1);
%    for i =1:row
%        low2 =size(aa_bow{i,:},2);
%        for j =1:low2
%            temp = strcmp(aa_bow{i,j},e);
%            ind = find(temp==1);
%            counter{i}(ind) = counter{i}(ind)+1; 
%        end       
%    end
   
end

