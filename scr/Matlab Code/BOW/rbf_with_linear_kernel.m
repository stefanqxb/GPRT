function kernel_matrix=rbf_with_linear_kernel(descr,gamma,m)
[len,row] = size(descr,1);
kernel_1 = zeros(len,row);
kernel_2 = zeros(len,row);
for i=1:len
    for j=1:i-1
        aaaa=exp(-gamma * sum((descr(i,:)-descr(j,:)).^2));
        kernel_1(i,j)=aaaa;
        kernel_1(j,i)=aaaa;
    end
end

for i=1:len
    for j=1:i-1
        aaaa=(descr(i,:)-descr(j,:)+1)^3;
        kernel_2(i,j)=aaaa;
        kernel_2(j,i)=aaaa;
    end
end

kernel_matrix = (1-m)*kernel_1+m*kernel_2;

end

