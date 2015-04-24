function [ ac_compostion ] = acid_compostion( peptide )
  ac_compostion = zeros(22,1);
  name = char(peptide);
  temp=cell(22);
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
  temp{21} = strfind(name,'O');
  temp{22} = strfind(name,'U');
for i=1:22
    ac_compostion(i) = size(temp{i},2);
end

end

