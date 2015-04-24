function [ rep_rete_mat,label,temp] = data_split_select( sort_rete,step,min_time)
   
rep_rete_mat = [];
label = [];
k =0:1:99;
time_point = round(min_time+step*k);
count = zeros(size(time_point));
t = size(time_point,2);
temp = zeros(size(sort_rete,1),2);
for i =1:size(sort_rete)
   for j =1:t
       if j == t
           if round(sort_rete(i)) >= time_point(j)
               count(j) = count(j)+ 1;
               temp(i,1) = sort_rete(i);
               temp(i,2) = j;
           end
       else
           if round(sort_rete(i)) >= time_point(j) && round(sort_rete(i))<time_point(j+1)
              count(j) = count(j)+ 1;
              temp(i,1) = sort_rete(i);
              temp(i,2) = j; %label
           end
       end
   end
end
count = round(15900*count/sum(count)); %normalization
for i=1:100
    [ind_x,~]= find(temp(:,2)==i);
    ram_value = randperm(size(ind_x));
    ram_value = ram_value(1:count(i));
    rep_rete_mat = [rep_rete_mat;sort_rete(ind_x(ram_value))];
    label =[label;i*ones(count(i),1)] ;
end   
end

