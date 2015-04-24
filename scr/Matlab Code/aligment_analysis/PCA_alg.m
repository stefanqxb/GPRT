function [ pro_mat ] = PCA_alg( or_mat )
   aver = mean(or_mat,1);
   [r,~] = size(or_mat);
   sum1 = 0;
   ac_sum =0;
   for i=1:r 
      sum1 = sum1+(or_mat(i,:)-aver)'*(or_mat(i,:)-aver);
   end
   
   S = sum1/(r-1);
   eig_val = eig(S,'nobalance');
   l = size(eig_val,1);
   temp = sort(eig_val,'descend');
   contributed_rate = temp/sum(temp);
   for i=1:l
      ac_sum = ac_sum+contributed_rate(i);
      if ac_sum>=0.8
          num = i;
          break;
      end
   end
   [eig_vec,~] = eig(S,'nobalance');
   pro_mat = eig_vec(:,end-num:end);
end

