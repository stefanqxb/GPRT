function [ ind ] = identify(a)
  refer = char('ARNDCQEGHILKMFPSTWYVOU ');
  ind = strfind(refer,a);
end

