1. run BOW_GP.py to build the model and make prediction (80% of the input will be used as training data to build the model and the rest 20% will be used as testing data)

2. Description:
   The RT predictor has two different models, one is based on SVR and the other is based on Gaussian Process.


3. Some parameters for extracting the features:
            (in selecting feature part)
     
     (1) 1--for BOW mehtod 2--for Elude method , else, error
     
     (2) 1--for 1-gram word 2--for 2-gram word (if 2 in (1), then it can be ramdon number) 
     
     (3) Choose a subset you want to test, the number is the amount of data, if choose the whole dataset, set 0 to this term  

4. Test data is in the file called  'data'