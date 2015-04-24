function [ chart ] = lookup_table( num )
    refer = char('ARNDCQEGHILKMFPSTWYVOU ');
    chart = refer(num);
end

