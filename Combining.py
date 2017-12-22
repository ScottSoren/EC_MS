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
import os #, sys    

def synchronize(Dataset_List, t_zero='start', append=None, cutit=False, 
                override=False, verbose=1):
    '''
    This will combine array data from multiple dictionaries into a single 
    dictionary with all time variables aligned according to absolute time.
    Data will be retained where the time spans overlap, unless cutit = False, in 
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
    
    if append is None: #by default, append if all datasets are same type, 17G28
        append = len({d['data_type'] for d in Dataset_List}) == 1
        override = True #added 17J28
    if verbose:
        print('append is ' + str(append))
        
                              #prepare to collect some data in the first loop:
    recstarts = []            #first recorded time in each file in seconds since midnight
    t_start = 0               #latest start time (start of overlap) in seconds since midnight
    t_finish = 60*60*24*7     #earliest finish time (finish of overlap) in seconds since midnight 
    t_first = 60*60*24*7      #earliest timestamp in seconds since midnight
    t_last = 0                #latest timestamp in seconds since midnight
    hasdata = {}              #'combining number' of dataset with 0 if its empty or 1 if it has data
    
    Combined_Data = {'data_type':'combined', 'data_cols':[]}
    title_combined = ''
    
    #go through once to generate the title and get the start and end times of the files and of the overlap
    if verbose:
        print('---------- syncrhonize entering first loop -----------')

    for nd, Dataset in enumerate(Dataset_List):
        Dataset['combining_number'] = nd
        if 'data_cols' not in Dataset or len(Dataset['data_cols']) == 0:
            print(Dataset['title'] + ' is empty')
            hasdata[nd] = 0
            recstarts += [60*60*24*7] # don't want dataset list to be shortened when sorted later according to recstarts!
            continue
        hasdata[nd] = 1    
        title_combined += Dataset['title'] + '__as_' + str(nd) + '__and___'
        if verbose:
            print('working on ' + Dataset['title'])
            print('timestamp is ' + Dataset['timestamp'])
        t_0 = timestamp_to_seconds(Dataset['timestamp'])
        
        t_f = 0
        t_s = 60*60*24*7
        
        for col in Dataset['data_cols']:
            if is_time(col):
                try:
                    t_s = min(t_s, t_0 + Dataset[col][0])   #earliest start of time data in dataset
                    t_f = max(t_f, t_0 + Dataset[col][-1])  #latest finish of time data in dataset
                except IndexError:
                    print(Dataset['title'] + ' may be an empty file.')
                    hasdata[nd] = 0
                    
        if hasdata[nd] == 0:
            continue
        recstarts += [t_s]               #first recorded time
    
        t_first = min([t_first, t_0])    #earliest timestamp  
        t_last = max([t_last, t_0])      #latest timestamp 
        t_start = max([t_start, t_s])    #latest start of time variable overall
        t_finish = min([t_finish, t_f])  #earliest finish of time variable overall
    
    title_combined = title_combined[:-6]
    Combined_Data['title'] = title_combined
    Combined_Data['tspan_0'] =    [t_start, t_finish] #overlap start and finish times as seconds since midnight
    Combined_Data['tspan_1'] = [t_start - t_first, t_finish - t_first]    # start and finish times as seconds since earliest start   

    if t_zero == 'start':
        t_zero = t_start
    elif t_zero == 'first':
        t_zero = t_first
    elif t_zero == 'last':
        t_zero = t_last
        
    Combined_Data['first'] = t_first - t_zero
    Combined_Data['last'] = t_last - t_zero
    Combined_Data['start'] = t_start - t_zero
    Combined_Data['finish'] = t_finish - t_zero
    
    Combined_Data['timestamp'] = seconds_to_timestamp(t_zero) #we want that timestamp refers to t=0, duh!

    Combined_Data['tspan'] = [t_start - t_zero, t_finish - t_zero]    #start and finish times of overlap as seconds since t=0  
    Combined_Data['tspan_2'] = Combined_Data['tspan'] #old code calls this tspan_2.
    
    if verbose:
        print('first: ' + str(t_first) + ', last: ' + str(t_last) + 
        ', start: ' + str(t_start) + ', finish: ' + str(t_finish))
        
    
    if t_start > t_finish and not override:
        print('No overlap. Check your files.\n')
        offerquit()
    
    print(hasdata)
    N_notempty = len([1 for v in hasdata.values() if v==1])
    if N_notempty == 0:
        print('First loop indicates that no files have data!!! Synchronize will return an empty dataset!!!')
    elif N_notempty == 1:
        print('First loop indicates that only one dataset has data! Synchronize will just return that dataset!')
        Combined_Data = next(Dataset_List[nd] for nd, v in hasdata.items() if v==1)
        print('\nfunction \'synchronize\' finished!\n\n')
        return Combined_Data
    
        #this is how I make sure to combine the datasets in the right order!     
    I_sort = np.argsort(recstarts)    
    Dataset_List = [Dataset_List[I] for I in I_sort]       
        #sort by first recorded absolute time. This is so that successive datasets of same type can be joined,
        #with the time variables increasing all the way through 
        #(note: EC lab techniques started together have same timestamp but different recstart)
    
    #and loop again to synchronize the data and put it into the combined dictionary.
    
    if (append and 
        'EC' in {d['data_type'] for d in Dataset_List}): #condition added 17I21
        #add a datacolumn that can be used to separate again afterwards by file
        #as of now (17G28) can't think of an obvious way to keep previous file numbers
        Combined_Data['file number'] = [] #'file_number' renamed 'file number' for consistency 17H09

    if verbose:
        print('---------- syncrhonize entering second loop -----------')
        
    for i, Dataset in enumerate(Dataset_List):
        nd = Dataset['combining_number']
        if hasdata[nd] == 0:
            if verbose:
                print('skipping this dataset, because its combining number, ' + str(nd) + ' is in the empty files list')
            continue
        elif verbose:
            print('proceding with this dataset, with combining number ' + str(nd))
            
        if verbose:
            print('cols in ' + Dataset['title'] + ':\n' + str(Dataset['data_cols']))
            print('cols in Combined_Data:\n' + str(Combined_Data['data_cols']))
        #this way names in Combined_Data match the order the datasets are input with
        t_0 = timestamp_to_seconds(Dataset['timestamp'])
        offset = t_0 - t_zero
        
        #first figure out where I need to cut, by getting the indeces striclty corresponding to times lying within the overlap
            
        I_keep = {}    
            #figure out what to keep:
        for col in Dataset['data_cols']:        
                #fixed up a bit 17C22, but this whole function should just be rewritten.
            if is_time(col):
                t = Dataset[col] + t_0 #absolute time
                I_keep[col] = [I for (I, t_I) in enumerate(t) if t_start < t_I < t_finish]
        
        #then cut, and put it in the new data set
        #first, get the present length of 'time/s' to see how much fill I need in the case of appending EC data from different methods.
        if append:
            print('append is True')
            if 'time/s' in Dataset:
                l1 = len(Dataset['time/s'])
                fn = np.array([i]*l1)
                if Dataset['data_type'] == 'EC':
                    Combined_Data['file number'] = np.append(Combined_Data['file number'], fn) 
                    print('len(Combined_Data[\'file number\']) = ' + str(len(Combined_Data['file number'])))
            else:
                print('\'time/s\' in Dataset is False')
            if 'time/s' in Combined_Data:
                l0 = len(Combined_Data['time/s'])     
            else:
                l0 = 0
        for col in Dataset['data_cols']:
            #print(col)
            data = Dataset[col] 
            if cutit:           #cut data to only return where it overlaps
                data = data[I_keep[get_time_col(col)]]  
                    #fixed up a bit 17C22, but this whole function should just be rewritten. 
            if is_time(col):
                data = data + offset
            
            if append:
                if col in Combined_Data:
                    data_0 = Combined_Data[col]
                else:
                    data_0 = np.array([])
                
                if get_time_col(col) == 'time/s': #rewritten 17G26 for fig3 of sniffer paper
                    #I got l0 before entering the present loop.
                    l2 = len(data_0)
                    if l0>l2:
                        fill = np.array([0]*(l0-l2))
                        #print('filling ' + col + ' with ' + str(len(fill)) + ' zeros')
                        data_0 = np.append(data_0, fill)
                
                #print('len(data_0) = ' + str(len(data_0)) + ', len(data) = ' + str(len(data)))
                data = np.append(data_0, data)
                #print('len(data_0) = ' + str(len(data_0)) + ', len(data) = ' + str(len(data)))
                    
            else:
                if col in Combined_Data:
                    print('conflicting versions of ' + col + '. adding subscripts.')
                    col = col + '_' + str(nd)                        
                    
            Combined_Data[col] = data
            if col not in Combined_Data['data_cols']:
                Combined_Data['data_cols'].append(col)  
        #offerquit()          
        
        #keep all of the metadata from the original datasets (added 16J27)
        for col, value in Dataset.items():
            if col not in Dataset['data_cols'] and col not in ['combining_number', 'data_cols']:     #fixed 16J29
                if col in Combined_Data.keys():
                    Combined_Data[col + '_' + str(nd)] = value
                else:
                    Combined_Data[col] = value
                    
    #17G28: There's still a problem if the last dataset is missing columns! Fixing now.  
    if 'time/s' in Combined_Data:
        Combined_Data['tspan_EC'] = [Combined_Data['time/s'][0], Combined_Data['time/s'][-1]] #this is nice to use
        l1 = len(Combined_Data['time/s'])
        for col in Combined_Data['data_cols']:   
            if get_time_col(col) == 'time/s': #rewritten 17G26 for fig3 of sniffer paper
                l2 = len(Combined_Data[col])
                if l1>l2:
                    fill = np.array([0]*(l1-l2))
                    print('filling ' + col + ' with ' + str(len(fill)) + ' zeros')
                    Combined_Data[col] = np.append(Combined_Data[col], fill)
        if append and 'file number' in Combined_Data.keys():
            if not 'file number' in Combined_Data['data_cols']:    #these somehow get added twice somehow, which sometimes gives problems
                Combined_Data['data_cols'].append('file number') #for new code            
            Combined_Data['file_number'] = Combined_Data['file number'] #for old code
            if not 'file_number' in Combined_Data['data_cols']:
                Combined_Data['data_cols'].append('file_number') #for old code

            
    if len(Combined_Data['data_cols']) == 0:
        print('The input did not have recognizeable data! Synchronize is returning an empty dataset!')    
        
    if verbose:
        print('function \'synchronize\' finsihed!\n\n')   
    
    return Combined_Data        


    
def cut(x, y, tspan=None, returnindeces=False, override=False):
    if tspan is None:
        return x, y
    
    if np.size(x) == 0:
        print('\nfunction \'cut\' received an empty input\n')
        if not override:
            offerquit()
        if returnindeces:
            return x, y, [] #I_keep #new 17H09
        return x, y        
    #I_keep = [I for (I, x_I) in enumerate(x) if tspan[0]<x_I<tspan[-1]] #list comprehensions are slow
    mask = np.logical_and(tspan[0]<x, x<tspan[-1]) #vectorized logic is fast
    
    if True not in mask and not override:
        print ('\nWarning! cutting like this leaves an empty dataset!\n' +
               'x goes from ' + str(x[0]) + ' to ' + str(x[-1]) + 
                ' and tspan = ' + str(tspan) + '\n')
        offerquit()
        
    x = x.copy()[mask]
    y = y.copy()[mask]
    
    if returnindeces:
        return x, y, mask #I_keep #new 17H09
    return x, y


def cut_dataset(dataset_0, tspan=None, verbose=True, override=False):
    '''
    Makes a time-cut of a dataset. Written 17H09.
    Unlike time_cut, does not ensure all MS data columns are the same length.
    '''
    if verbose:
        print('\n\nfunction \'cut dataset\' at your service!\n') 
    dataset = dataset_0.copy()
    if tspan is None:
        return dataset
    #print(dataset['title'])
    # I need to cunt non-time variables irst
    indeces = {} #I imagine storing indeces improves performance
    for col in dataset['data_cols']:
        if verbose:
            print(col + ', length = ' + str(len(dataset[col])))
        if is_time(col): #I actually don't need this. 
            continue #yes I do, otherwise it cuts it twice.
        timecol = get_time_col(col)
        if verbose:
            print(timecol + ', timecol length = ' + str(len(dataset[timecol])))
        if timecol in indeces.keys():
            if verbose:
                print('already got indeces, len = ' + str(len(indeces[timecol])))
            if timecol == col:
                continue
            dataset[col] = dataset[col].copy()[indeces[timecol]]
        else:
            if verbose:
                print('about to cut!')
            dataset[timecol], dataset[col], indeces[timecol] = \
                cut(dataset[timecol], dataset[col],
                    tspan=tspan, returnindeces=True, override=override) 
    tspan0 = dataset['tspan']
    dataset['tspan'] = [max([tspan[0], tspan0[0]]), min([tspan[-1], tspan0[-1]])]
    if verbose:   
        print('\nfunction \'cut dataset\' finsihed!\n\n') 
    return dataset


def offerquit():
    yn = input('continue? y/n\n')
    if yn == 'n':
        raise SystemExit
    
def time_cut(MS_Data_0, tspan, verbose=1):
    '''
    cuts an MS data set, retaining the portion of the data set within a specified
    time interval. Does nothing to the EC data.
    # antiquated. Replaced by cut_dataset
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
  
    
def is_time(col, verbose=0):
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

def is_MS_data(col, verbose=False):
    if re.search(r'^M[0-9]', col):
        return True
    return False

def is_EC_data(col, verbose=False):
    return not is_MS_data(col)
    
def get_type(col, verbose=False):
    if is_MS_data(col, verbose):
        return 'MS'
    if is_EC_data(col, verbose):
        return 'EC'
    return None

def get_time_col(col, verbose=False):
    if is_time(col):
        time_col = col
    elif is_MS_data(col):
        time_col = col.replace('-y','-x')
    elif is_EC_data(col): 
        time_col = 'time/s'
    else:
        print('don\'t know what ' + col + ' is or what it\'s time col is.')
        time_col = None
    if verbose:
        print('\'' + col + '\' should correspond to time col \'' + str(time_col) +'\'')
    return time_col

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

def dayshift(dataset, days=1):
    ''' Can work for up to 4 days. After that, hh becomes hhh... 
    '''
    dataset['timestamp'] = seconds_to_timestamp(timestamp_to_seconds(dataset['timestamp']) + days*24*60*60) 
    return dataset

def sort_time(dataset, data_type='EC', verbose=False):
    #17K11: This now operates on the original dictionary, so
    #that I don't need to read the return.
        if verbose:
            print('\nfunction \'sort_time\' at your service!\n\n')
        
        if 'NOTES' in dataset.keys():
            dataset['NOTES'] += '\nTime-Sorted\n'
        else: 
            dataset['NOTES'] = 'Time-Sorted\n'
        
        if data_type == 'all':
            data_type = ['EC','MS']
        elif type(data_type) is str:
            data_type = [data_type]

        sort_indeces = {} #will store sort indeces of the time variables
        data_cols = dataset['data_cols'].copy()
        dataset['data_cols'] = []
        for col in data_cols:
            if verbose:
                print('working on ' + col)
            data = dataset[col] #do I need the copy?
            if get_type(col) in data_type: #retuns 'EC' or 'MS', else I don't know what it is.
                time_col = get_time_col(col, verbose)
                if time_col in sort_indeces.keys():
                    indeces = sort_indeces[time_col]
                else:
                    print('getting indeces to sort ' + time_col)
                    indeces = np.argsort(dataset[time_col])
                    sort_indeces[time_col] = indeces
                if len(data) != len(indeces):
                    if verbose:
                        print(col + ' is not the same length as its time variable!\n' +
                              col + ' will not be included in the time-sorted dataset.')
                else:
                    dataset[col] = data[indeces]
                    dataset['data_cols'] += [col]
                    print('sorted ' + col + '!')
            else: #just keep it without sorting.
                dataset['data_cols'] += [col]
                dataset[col] = data
                

        if verbose:
            print('\nfunction \'sort_time\' finished!\n\n')    
        
        #return dataset#, sort_indeces  #sort indeces are useless, 17J11
    #if I need to read the return for normal use, then I don't want sort_indeces
    
    

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
    
    
    