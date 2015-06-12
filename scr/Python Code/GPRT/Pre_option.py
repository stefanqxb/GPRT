__author__ = 'Administrator'

def option():
    request1  = raw_input('Choose a type of feature extraction method 1-BOW, 2-Elude:')
    if request1 == '1':
        request12  = raw_input('Choose the number of gram for BOW 1 or 2:')
    elif request1 == '2':
        request12 = 0
    request2  = raw_input('Number of subset, set 0 for the whole dataset:')
    return request1, request12, request2