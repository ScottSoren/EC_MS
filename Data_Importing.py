# -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 10:40:46 2016

@author: soren

This file holds methods for importing and organizing experimental data. Made as
part of EC_MS, but saved seperately just because spyder3 doesn't have code folding
and I don't want to get lost in the scroll bar.

"""

import os
import re
import codecs
from copy  import deepcopy



def import_text(full_path_name='default', verbose=1):   
    '''
    This method will import the full text of a file selected by user input as a 
    list of lines
    '''    
    if verbose:
        print('\n\nfunction \'import_text\' at your service!')
    
    if full_path_name == 'input':
        full_path_name = input('Enter full path for the file name as \'directory/file.extension\'')
    elif full_path_name == 'default':
        full_path_name = '/home/soren/Desktop/Sniffer_Experiments/O18_NiNPs/00_python/test_files/02_O16_to_O18_10_CA_C01.mpt'
    
    if full_path_name[0] == '~':
        full_path_name = '/home/soren' + full_path_name[1:]
        
    file_name_object = re.search(r'/[^/]*\Z', full_path_name)
    file_name = file_name_object.group()[1:]
    directory_name = full_path_name[:file_name_object.start()]

    original_directory = os.getcwd()
    os.chdir(directory_name)
    ls_string = str(os.listdir())    
    if verbose:
        print('\n' + full_path_name + '\n ls: \n' + ls_string + '\n')
    if len(file_name)==0:
        print('\n' + full_path_name + '\n ls: \n' + ls_string + '\n')
        file_name = input('Directory given. Enter the full name of the file to import\n')
    if verbose:   
        print('importing data from ' + file_name )


    possible_encodings = ['utf8','iso8859_15']  
    #mpt files seem to be the latter encoding, even though they refer to themselves as ascii
    for encoding_type in possible_encodings:    
        try:
            file_object = codecs.open(file_name, 'r', encoding = encoding_type)
            file_lines = file_object.readlines()
            if verbose:
                print('Was able to readlines() with encoding ' + encoding_type)
            break
        except UnicodeDecodeError:
            if verbose:
                print('Shit, some encoding problem in readlines() for ' + encoding_type)
    else:
        print('couldn\'t read ' + file_name + '\n ... may by due to an encoding issue')
        
    os.chdir(original_directory)
    if verbose:
        print('function \'import_text\' finished!\n\n')
    return file_lines




def text_to_data(file_lines, title='get_from_file',
                 data_type='EC', N_blank=10, verbose=1):
    '''
    This method will organize data in the lines of text from an EC or MS file 
    into a dictionary as follows (plus a few more keys)
    {'title':title, 'header':header, 'colheader1':[data1], 'colheader2':[data2]...}
    '''    
    if verbose:
        print('\n\nfunction \'text_to_data\' at your service!')
    
    #disect header
    N_lines = len(file_lines)        #number of header lines
    N_head = N_lines                 #this will change when I find the line that tells me how ling the header is
    header_string = ''              
    n_blank = 0                      #header ends with a number of blank lines in MS
        
    DataDict = {}
    
    for nl, line in enumerate(file_lines):
        
        if nl < N_head - 1:            #we're in the header
        
            if data_type == 'EC':
                if title == 'get_from_file':
                    if re.search('File :',line):
                        title_object = re.search(r'[\S]*\Z',line.strip())
                        title = title_object.group()
                        if verbose:
                            print('name \'' + title + '\' found in line ' + str(nl))
                if re.search(r'[Number ]*header lines',line):
                    N_head_object = re.search(r'[0-9][0-9]*',line)
                    N_head = int(N_head_object.group())
                    if verbose:
                        print('N_head \'' + str(N_head) + '\' found in line ' + str(nl))
                elif re.search('Acquisition started',line):
                    timestamp_object = re.search(r'[\S]*\Z',line.strip())
                    timestamp = timestamp_object.group()
                    if verbose:
                        print('timestamp \'' + timestamp + '\' found in line ' + str(nl))        
                header_string = header_string + line
            
            elif data_type == 'MS':
                if len(line.strip())==0:
                    n_blank += 1
                    if n_blank>N_blank and len(file_lines[nl+1].strip())>0:
                        N_head = nl+2
                else:
                    n_blank = 0    
                    if title == 'get_from_file':
                        object1 = re.search(r'"Comment"[\s]*"[^"]*',line)
                        if object1:
                            string1 = object1.group()
                            title_object = re.search(r'[\S]*\Z',string1.strip())
                            title = title_object.group()[1:]
                            if verbose:
                                print('name \'' + title + '\' found in line ' + str(nl))
                    object2 = re.search(r'"Recorded at"[\s]*"[^"]*',line)
                    if object2:
                        string2 = object2.group()
                        timestamp_object = re.search(r'[\S]*\Z',string2.strip()) 
                        timestamp = timestamp_object.group()
                        if verbose:
                            print('timestamp \'' + timestamp + '\' found in line ' + str(nl))                 
                header_string = header_string + line
                
            
        elif nl == N_head - 1:      #then it is the column-header line
               #(EC-lab includes the column-header line in header lines)
            #col_header_line = line
            col_headers = line.strip().split('\t')
            DataDict['N_col']=len(col_headers) 
            DataDict['data_cols'] = deepcopy(col_headers)   #will store names of columns containing data
            for col in col_headers:
                DataDict[col] = []              #data will go here    
            header_string = header_string+line #include this line in the header
            if verbose:
                print('Data starting on line ' + str(N_head) + '\n')
            

        else:                   # data, baby!
            line_data = line.strip().split('\t')
            if not len(line_data) == len(col_headers):
                if verbose:
                    print(list(zip(col_headers,line_data)))
                    print('Mismatch between col_headers and data on line ' + str(nl) + ' of ' + title)
                if nl == N_lines - 1:
                    print('mismatch due to an incomplete last line of ' + title + '. I will discard the last line.')
                    break
            for col, data in zip(col_headers,line_data):
                if col in DataDict['data_cols']:
                    try:
                        data = float(data)
                    except ValueError:
                        if verbose:
                            print(list(zip(col_headers,line_data)))
                            print(title + ' in text_to_data: \nRemoved \'' + str(col) +'\' from data columns because of value \'' + 
                                str(data) + '\' at line ' + str(nl) +'\n')
                        DataDict['data_cols'].remove(col)
                            
                DataDict[col].append(data)
                    
            
    DataDict['title'] = title
    DataDict['header'] = header_string
    DataDict['timestamp'] = timestamp
    DataDict['data_type'] = data_type

    if verbose:
        print('function \'text_to_data\' finished!\n\n')    
    return DataDict


def import_data(full_path_name='default', title='get_from_file',
                 data_type='EC', N_blank=10, verbose=1):
                 
    file_lines = import_text(full_path_name, verbose)
    DataDict = text_to_data(file_lines, title, data_type, N_blank, verbose)

    return DataDict



if __name__ == '__main__':
    
    default_directory = '/home/soren/Desktop/Sniffer_Experiments/O18_NiNPs/00_python/test_files/'    
    
    CA_file = default_directory + '02_O16_to_O18_10_CA_C01.mpt'
    CV_file = default_directory + '02_O16_to_O18_06_CVA_C01.mpt'
    MS_file = default_directory + 'QMS_data.txt'
    
    CA_Data = import_data(CA_file)
    CV_Data = import_data(CV_file)
    MS_Data = import_data(MS_file,data_type='MS')
    
    
    
    