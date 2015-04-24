function [ time_interval ] = time_95_diff(hist,max_t,min_t,step)
    [~,n] = size(hist);
    total = sum(hist);
    threshold = 0.95*total;
    count =0;
    for i=1:n
        count = hist(i) + count ;
        if count >= threshold
            counter =i;
            break;
        end
    end
    
    time_interval = counter*(max_t-min_t)/step;
    time_interval = time_interval / max_t ;
    
end

