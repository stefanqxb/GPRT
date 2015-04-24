function [ ac_compostion ] = acid_compostion( peptide )
  [m,~] = size(peptide);
  ac_compostion = zeros(m,20);
  for i =1:m
  name = char(peptide{i});
  temp=cell(20);
  temp{1} = strfind(name,'A');
  temp{2} = strfind(name,'R');
  temp{3} = strfind(name,'N');
  temp{4} = strfind(name,'D');
  temp{5} = strfind(name,'C');
  temp{6} = strfind(name,'Q');
  temp{7} = strfind(name,'E');
  temp{8} = strfind(name,'G');
  temp{9} = strfind(name,'H');
  temp{10} = strfind(name,'I');
  temp{11} = strfind(name,'L');
  temp{12} = strfind(name,'K');
  temp{13} = strfind(name,'M');
  temp{14} = strfind(name,'F');
  temp{15} = strfind(name,'P');
  temp{16} = strfind(name,'S');
  temp{17} = strfind(name,'T');
  temp{18} = strfind(name,'W');
  temp{19} = strfind(name,'Y');
  temp{20} = strfind(name,'V');
     for j=1:20
         ac_compostion(i,j) = size(temp{j},2);
     end
      
  end
 

end

