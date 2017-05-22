# -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 10:40:46 2016
Most recently edited: 16J27
@author: Scott

This is the core file of the package. Includes functions for combining EC and MS 
"""

# make python2-compatible:
from __future__ import print_function
from __future__ import division

import numpy as np
import re
import os, sys    

def synchronize(Dataset_List, t_zero='start', append=1, cutit=0, 
                override=False, verbose=1):
    '''
    This will combine array data from multiple dictionaries into a single 
    dictionary with all time variables aligned according to absolute time.
    Data will be retained where the time spans overlap, unless cutit = 0, in 
    which case all data will be retained, but with t=0 at the start of the overlap.
    if t_zero is specified, however, t=0 will be set to t_zero seconds after midnight
    (added 16J29:) if append=1, data columns of the same name will be joined
    '''
    if verbose:
        print('\n\nfunction \'synchronize\' at your service!')
        
    if type(Dataset_List) is dict:
        print('''The first argument to synchronize should be a list of datasets! 
                You have instead input a dictionary as the first argument. 
                I will assume that the first two arguments are the datasets you
                would like to synchronize with standard settings.''')
        Dataset_List = [Dataset_List, t_zero]
        t_zero = 'start'
    
                              #prepare to collect:
    recstarts = []            #first recorded time in each file in seconds since midnight
    t_start = 0               #latest start time (start of overlap) in seconds since midnight
    t_finish = 60*60*24*7     #earliest finish time (finish of overlap) in seconds since midnight 
    t_first = 60*60*24*7      #earliest timestamp in seconds since midnight
    t_last = 0                #latest timestamp in seconds since midnight

    
    Combined_Data = {'data_type':'combined', 'data_cols':[]}
    title_combined = ''
    
    #go through once to generate the title and get the start and end times
    for nd, Dataset in enumerate(Dataset_List):
        Dataset['combining_number'] = nd
        title_combined += Dataset['title'] + '__as_' + str(nd) + '__and___'
        #Dataset = numerize(Dataset)    #16I28: the dDataset should already be numerized by importdata
        
        t_0 = timestamp_to_seconds(Dataset['timestamp'])
        
        t_f = 0
        t_s = 60*60*24*7
        
        for col in Dataset['data_cols']:
            if is_time(col):
                t_s = min(t_s, t_0 + Dataset[col][0])   #earliest start of time data in dataset
                t_f = max(t_f, t_0 + Dataset[col][-1])  #latest finish of time data in dataset
        
        recstarts += [t_s]               #first recorded time
    
        t_first = min([t_first, t_0])    #earliest timestamp  
        t_last = max([t_last, t_0])      #latest timestamp 
        t_start = max([t_start, t_s])    #latest start of time variable overall
        t_finish = min([t_finish, t_f])  #earliest finish of time variable overall
    
    title_combined = title_combined[:-6]
    Combined_Data['title'] = title_combined
    Combined_Data['timestamp'] = seconds_to_timestamp(t_start)
    Combined_Data['tspan'] =    [t_start, t_finish] #overlap start and finish times as seconds since midnight
    Combined_Data['tspan_1'] = [t_start - t_first, t_finish - t_first]    # start and finish times as seconds since earliest start
    if t_zero == 'start':
        t_zero = t_start
    elif t_zero == 'first':
        t_zero = t_first
    elif t_zero == 'last':
        t_zero = t_last
    if verbose:
        print('first: ' + str(t_first) + ', last: ' + str(t_last) + 
        ', start: ' + str(t_start) + ', finish: ' + str(t_finish))
    Combined_Data['tspan_2'] = [t_start - t_zero, t_finish - t_zero]    #start and finish times of overlap as seconds since zero point   
    
    if t_start > t_finish and not override:
        print('No overlap. Check your files.\n')
        offerquit()
    
    I_sort = np.argsort(recstarts)
    Dataset_List = [Dataset_List[I] for I in I_sort]       
        #sort by first recorded absolute time. This is so that successive datasets of same type can be joined,
        #with the time variables increasing all the way through 
        #(note: EC lab techniques started together have same timestamp but different recstart)
    
    #and again to synchronize the data and put it into the combined dictionary
    for Dataset in Dataset_List:
        nd = Dataset['combining_number']           
        #this way names in Combined_Data match the order the datasets are input with
        t_0 = timestamp_to_seconds(Dataset['timestamp'])
        offset = t_0 - t_zero
        
        #first figure out where I need to cut, by getting the indeces striclty corresponding to times lying within the overlap
            #a kind of ugly way to get the number of rows in dataset:
        I_keep = [I for (I, v) in enumerate(Dataset[Dataset['data_cols'][0]])]     
            #figure out what to keep:
        for col in Dataset['data_cols']:
            if is_time(col):
                t = Dataset[col] + t_0 #absolute time
                I_keep = [I for I in I_keep if t_start < t[I] < t_finish]
        
        #then cut, and put it in the new data set
        for col in Dataset['data_cols']:
            data = Dataset[col] 
            if cutit:           #cut data to only return where it overlaps
                data = data[I_keep] 
            if is_time(col):
                data = data + offset
            if col in Combined_Data:
                if append:
                    Combined_Data[col] = np.append(Combined_Data[col], data)
                    continue                    
                col = col + '_' + str()
            Combined_Data[col] = data
            Combined_Data['data_cols'].append(col)  
            
        #keep all of the metadata from the original datasets (added 16J27)
        for col, value in Dataset.items():
            if col not in Dataset['data_cols'] and col not in ['combining_number', 'data_cols']:     #fixed 16J29
                if col in Combined_Data.keys():
                    Combined_Data[col + '_' + str(nd)] = value
                else:
                    Combined_Data[col] = value
           
    if verbose:
        print('function \'synchronize\' finsihed!\n\n')   
    
    return Combined_Data        

    
def cut(x, y, tspan):
    if tspan is None:
        return x, y
    I_keep = [I for (I, x_I) in enumerate(x) if tspan[0]<x_I<tspan[-1]]
    x = x.copy()[I_keep]
    y = y.copy()[I_keep]
    if np.size(x) == 0:
        print('\nfunction \'cut\' received an empty input\n')
        offerquit()
    elif np.size(I_keep) == 0:
        print ('\nWarning! cutting like this leaves an empty dataset!\n' +
               'x goes from ' + str(x[0]) + ' to ' + str(x[-1]) + 
                ' and tspan = ' + str(tspan) + '\n')
        offerquit()
    return x, y

def offerquit():
    yn = input('continue? y/n\n')
    if yn == 'n':
        raise SystemExit
    
def time_cut(MS_Data_0, tspan, verbose=1):
    '''
    cuts an MS data set, retaining the portion of the data set within a specified
    time interval. Does nothing to the EC data.
    '''
    if verbose:
        print('\n\nfunction \'time_cut\' at your service! \n Time cutting ' + MS_Data_0['title'])
    MS_Data = MS_Data_0.copy() #otherwise I cut the original dataset!
    length = 10**9
    for col in MS_Data['data_cols']:
        if is_time(col):
            x = MS_Data[col].copy()
            if col[-2:] == '-x': 
                ycols = [col[:-2] +'-y',]      #then this is an MS time variable. 
                #This assumes column names haven't been fucked with by synchronization
            else:           #then MS_data is in fact MS_and_EC data
                #break       #will need to change this when EC file lasts longer than tspan
                ycols = [c for c in MS_Data['data_cols'] if c != col if c[-2:] != '-x' if c[-2:] != '-y']
                #Assuming that everything that doesn't end in -x and -y is EC data
            
            I_start_object = np.where(x>tspan[0])
            if len(I_start_object[0])>0:
                if verbose:
                    print('cutting ' + col + ' at start')
                I_start = I_start_object[0][0]
            else:
                I_start = 0
            I_finish_object = np.where(x>tspan[1])
            if len(I_finish_object[0])>0:
                if verbose:
                    print('cutting ' + col + ' at finish')
                I_finish = I_finish_object[0][0]
            else:
                I_finish = len(x)
            x = x[I_start:I_finish]
            MS_Data[col]=x
            for ycol in ycols:
                y = MS_Data[ycol].copy()
                y = y[I_start:I_finish]       
                MS_Data[ycol]=y
            if col[-2:] == '-x':
                length = min(length, len(x))
            #and now, to make sure all of the QMS columns are still the same length:
    for col in MS_Data['data_cols']:
        if col[-2:] == '-x' or col[-2:] == '-y':
            MS_Data[col] = MS_Data[col][:length]
    if verbose:
        print('function \'time_cut\' finished!\n\n')    
    return MS_Data
  
    
def is_time(col, verbose = 0):
    '''
    determines if a column header is a time variable, 1 for yes 0 for no
    '''
    if verbose:
        print('\nfunction \'is_time\' checking \'' + col + '\'!')
    if col[0:4]=='time':
        return 1
    if col[-2:]=='-x': 
        return 1         
    #in case it the time marker is just burried in a number suffix:
    ending_object = re.search(r'_[0-9][0-9]*\Z',col) 
    if ending_object:
        col = col[:ending_object.start()]
        return is_time(col)
    if verbose:
        print('...not time')
    return 0

def timestamp_to_seconds(timestamp):
    '''
    seconds since midnight derived from timestamp hh:mm:ss
    '''
    h = int(timestamp[0:2])
    m = int(timestamp[3:5])
    s = int(timestamp[6:8])
    seconds = 60**2 *h + 60 *m + s
    return seconds
    
def seconds_to_timestamp(seconds):
    '''
    timestamp hh:mm:ss derived from seconds since midnight
    '''
    h = int(seconds/60**2)
    seconds = seconds - 60**2 *h
    m = int(seconds/60)
    seconds = seconds - 60 *m
    s = int(seconds)
    timestamp = '{0:2d}:{1:2d}:{2:2d}'.format(h,m,s)
    timestamp = timestamp.replace(' ','0')
    return timestamp



if __name__ == '__main__':
    
    from Data_Importing import import_data #honestly, I would just have everything in one module if you could fold code in spyder3
    from Plotting import plot_vs_time
    import os    
    
    default_directory = os.path.abspath(os.path.join(os.getcwd(), os.pardir))     
#    
    CA_Data = import_data(default_directory, data_type='EC')
    MS_Data = import_data(default_directory, data_type='MS')
    
    CA_and_MS = synchronize([CA_Data,MS_Data], cutit = 1)
    
    plot_vs_time(CA_and_MS,
                 cols_1=[('M4-x','M4-y'),('M18-x','M18-y'),('M28-x','M28-y'),('M32-x','M32-y')],  
#                 cols_2=[('time/s', '<I>/mA'),],       #for CV
                 cols_2=[('time/s', 'I/mA'),],          #for CA
                 )
    
    
    