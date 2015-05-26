__author__ = 'Administrator'

def mini_time_win(hist,max_t,min_t,step,max_total_t):
    total = sum(hist)
    threshold = round(0.95 * total)
    count = 0
    counter = 0
    for i in range(len(hist)):
        count = hist[i] + count
        if count >= threshold:
            counter = i
            break

    time_interval = counter * (max_t - min_t)/(step*max_total_t)
    return  time_interval