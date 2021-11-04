# -*- coding: utf-8 -*-
"""
Created on Sun Sep  5 13:07:47 2021

Run linear decoder 

@author: CT
"""
#%% Import modules

import os
import scipy.io
import numpy as np #general math stuff
import matplotlib.pyplot as plt #some plotting stuff
from sklearn import svm #This allows us to use non linear svm models if needed
from sklearn.linear_model import Perceptron #not used but another fitting technique we could use if needed
from sklearn.svm import LinearSVC #This is the linear svc model from scikit learn, can 
from sklearn.model_selection import cross_validate, ShuffleSplit #Get stuff for cross validations

## Import data
os.chdir('D:/Deuteron Data/Ready to analyze output/Amos_2021-07-29/')
full_data = scipy.io.loadmat('SVM_input.mat')

## Run model
my_labels = full_data['Labels'] #Put labels for each time point here, labels should be numeric I think
my_labels = np.squeeze(my_labels)
my_labels_shuffled = full_data['Labels_shuffled'] 
my_labels_shuffled = np.squeeze(my_labels_shuffled)
#idx = np.squeeze(np.logical_or(my_labels ==6, my_labels==5, my_labels==4))
#my_labels_no_zero = np.squeeze(my_labels[idx])

my_data = full_data['Input_matrix'] #Put neural data here, needs to be time x neurons
#my_data_no_zero = my_data[idx,:]


linSVC = LinearSVC(max_iter = 100000, tol = 0.001)#, class_weight = 'balanced') #Perceptron(max_iter = 10000, tol = 0.001) #from scikit (presumably) #set up Perceptron classifier

crossval = ShuffleSplit(n_splits = 10, test_size = 0.20)

cv_results = cross_validate(linSVC, my_data, my_labels, cv=crossval) #This stores results as dictionary
 
avg_score = np.mean(cv_results['test_score']) #take average of the cv_results to 
 
print(avg_score)

#Notes: 
# Groom receive, groom given, foraging, HIP, HIS, SP, SS yields 94% accuracy when considering both arrays. 
# When only considering array 1 (TEO), accuracy drops to 93%. When considering array #2, accuracy remains at 95%.
#Adding scratch affects performance (drops from 94% to 86%, although chance level drops to 1/8)
#Adding self-groom drives performance down to 91%
#When considering 10 classes of behaviorm accuracy is still really high (82% accuracy) when considering both arrays.

#%% Quick plotting



plot_splits = 1

#note this can create a lot of figures as it creates out one for each split.
    

if plot_splits: 
    for train_index, test_index in crossval.split(my_data):
    
    
        plt.figure()
        plt.hist(my_labels[train_index])
        plt.title('training data distribution for each split')
    
    
    for train_index, test_index in crossval.split(my_data):
    
    
        plt.figure()
        plt.hist(my_labels_no_zero[test_index])
        plt.title('testing data distribution for each split')

