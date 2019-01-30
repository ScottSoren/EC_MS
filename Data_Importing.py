#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 10 03:17:41 2017

@author: scott
"""
import os, platform, re, codecs
import time, datetime, pytz  #all three seem necessary for dealing with timezones.
import numpy as np

float_match = '[-]?\d+[\.]?\d*(e[-]?\d+)?'     #matches floats like '-3.54e4' or '7' or '245.13' or '1e-15'
#note, no white space included on the ends! Seems to work fine.
timestamp_match = '([0-9]{2}:){2}[0-9]{2}'     #matches timestamps like '14:23:01'
date_match = '([0-9]{2}[/-]){2}[0-9]{4}'          #matches dates like '01/15/2018' or '09-07-2016'
#^  mm/dd/yyyy, as EC lab does
date_match_2 = '[0-9]{4}([/-][0-9]{2}){2}'          #matches dates like '01/15/2018' or '09-07-2016'
#^  yyyy/mm/dd, as cinfdata does
# older EC-Lab seems to have dashes in date, and newer has slashes.
# Both seem to save month before day, regardless of where the data was taken or .mpt exported.


def numerize(data):
    for col in data['data_cols']: #numerize!
        data[col] = np.array(data[col])

def get_empty_set(cols, **kwargs):
        #get the colheaders and make space for data
    data = {}
    data.update(kwargs)
    for col in cols:
        data[col] = []
    data['data_cols'] = cols
    return data

def parse_timezone(tz=None):
    '''
    Gets a timezone object from a timezone string. Includes some abbreviations
    useful for Scott. If the input is not a string, it is returned as is.
    '''
    abbreviations = {'CA':'US/Pacific', 'DK':'Europe/Copenhagen',}
    if tz in abbreviations:
        return pytz.timezone(abbreviations[tz])
    elif type(tz) is str:
        return pytz.timezone(tz)
    else:
        return tz

def timestamp_to_epoch_time(timestamp, date='today', tz=None, verbose=True):
    '''
    Possibly overly idiot-proof way to convert a number of timestamps read
    from my data into epoch (unix) time.

    tz is the Timezone, which is only strictly necessary when synchronizing
    data taken at different places or accross dst at a place with different dst
    implementation than the local (dst is stupid!).
    If timezone is a number, it is interpreted as the offset from GMT of the
    data.

    The epoch time is referred to here and elsewhere as tstamp.
    '''
    if tz is not None:
        tz = parse_timezone(tz)
        if verbose:
            print('getting epoch time given a timestamp local to ' + str(tz))
        epoch = pytz.utc.localize(datetime.datetime.utcfromtimestamp(0))
    if timestamp == 'now':
        return time.time()
    elif type(timestamp) is time.struct_time:
        if verbose:
            print('\'timestamp_to_unix_time\' revieved a time.struct_time object. ' +
                  'Returning the corresponding epoch time.')
    elif type(timestamp) is not str:
        if verbose:
            print('timestamp_to_unix_time\' didn\'t receive a string. ' +
                  'Received:\n' + str(timestamp) + ' \nReturning the argument.')
        return timestamp
    if len(timestamp) > 8 and date=='today':
        if verbose:
            print('\'timestamp_to_unix_time\' is assuming' +
                  ' the date and timestamp are input in the same string.')
        try:
            if tz is None:
                struct = time.strptime(timestamp)
                tstamp = time.mktime(struct)
            else:
                dt_naive = datetime.datetime.strptime(timestamp, '%a %b %d %H:%M:%S %Y')
                dt = tz.localize(dt_naive)
                tstamp = (dt - epoch).total_seconds()
            if verbose:
                print('timestamp went straight into time.strptime()! Returning based on that.')
            return tstamp
        except ValueError:
            if verbose:
                print('bit \'' + timestamp + '\' is not formatted like ' +
                      'time.strptime() likes it. Checking another format')
            try:
                date = re.search(date_match, timestamp).group()
                if verbose:
                    print('Matched the date with \'' + date_match + '\'.')
            except AttributeError:
                if verbose:
                    print('Couldn\'t match with \'' + date_match +
                          '\'. Assuming you want today.')
            try:
                timestamp = re.search(timestamp_match, timestamp).group()
                if verbose:
                    print('Matched the time with \'' + timestamp_match + '\'.' )
            except AttributeError:
                if verbose:
                    print('I got no clue what you\'re talking about, dude, ' +
                          'when you say ' + timestamp + '. It didn\'t match \.' +
                          timestamp_match + '\'. Assuming you want 00:00:00')
                timestamp = '00:00:00'
    #h, m, s = (int(n) for n in timestamp.split(':'))
    #D, M, Y = (int(n) for n in re.split('[-/]', date))
    if date == 'today':
        date = time.strftime('%m/%d/%Y')
    if tz is None:
        if '-' in date:
            date = date.replace('-','/') #18D08
        struct = time.strptime(date + ' ' + timestamp, '%m/%d/%Y %H:%M:%S')
        tstamp = time.mktime(struct)
    else:
        dt_naive = datetime.datetime.strptime(date + ' ' + timestamp, '%m/%d/%Y %H:%M:%S')
        dt = tz.localize(dt_naive)
        tstamp = (dt - epoch).total_seconds()

    return tstamp

def epoch_time_to_timestamp(tstamp, tz=None, verbose=True):
    '''
    tz is the Timezone, which is only strictly necessary when synchronizing
    data taken at different places or accross dst at a place with different dst
    implementation than the local (dst is stupid!).
    If timezone is a number, it is interpreted as the offset from GMT of the
    data in hours (+1 for Denmark, -8 for California)
    '''
    if tz is None:
        struct = time.localtime(tstamp)
    else:
        tz = parse_timezone(tz)
        if verbose:
            print('getting the timestamp local to ' + str(tz) + ' from epoch time.')
        dt_utc = datetime.datetime.utcfromtimestamp(tstamp)
        dt_tz = tz.fromutc(dt_utc)
        struct = dt_tz.timetuple()
    hh = str(struct.tm_hour)
    if len(hh) == 1:
        hh = '0' + hh
    mm = str(struct.tm_min)
    if len(hh) == 1:
        mm = '0' + mm
    ss = str(struct.tm_sec)
    if len(hh) == 1:
        ss = '0' + ss
    timestamp = hh + ':' + mm + ':' + ss
    return timestamp

def timetag_to_timestamp(filename):
    '''
    Converts a time tag of format _<hh>h<mm>m<ss>_ to timestamp of format
    <hh>:<mm>:<ss>, which is what synchronize reads (I should maybe change
    synchronize to work with the unix epoch timestamp instead...)
    The time tag is something we write in the file names to give an approximate
    sense of when a measurement is started. It can be used as the measurement
    start time if there really is no better way.
    I can't beleive SPEC doesn't save time. I will pressure them to fix this.
    '''
    hm_match = re.search(r'_[0-9]{2}h[0-9]{2}', filename)
    hm_str = hm_match.group()
    hh = hm_str[1:3]
    mm = hm_str[-2:]
    ss_match = re.search(r'_[0-9]{2}h[0-9]{2}m[0-9]{2}', filename)
    if ss_match is None:
        ss = '00'
    else:
        ss = ss_match.group()[-2:]
    return hh + ':' + mm + ':' + ss


def get_creation_timestamp(filepath):
    '''
    Returns creation timestamp of a file in the format that
    combining.syncrhonize reads.
    The timestamp is local time, not absolute time.
    We need to move to epoch time everywhere!!!
    '''
    t = get_creation_time(filepath)
    struct = time.localtime(t)
    hh = str(struct.tm_hour)
    mm = str(struct.tm_minute)
    ss = str(struct.tm_second)
    return hh + ':' + mm + ':' + ss


def get_creation_time(filepath, verbose=True):
    """
    Try to get the date that a file was created, falling back to when it was
    last modified if that isn't possible.
    See http://stackoverflow.com/a/39501288/1709587 for explanation.
    """
    if platform.system() == 'Windows':
        tstamp = os.path.getctime(filepath)
        if verbose:
            print('In Windows. Using os.path.getctime(\'' + filepath + '\') as tstamp.')
    else:
        stat = os.stat(filepath)
        try:
            tstamp = stat.st_birthtime
            if verbose:
                print('In linux. Using os.stat(\'' + filepath + '\').st_birthtime as tstamp.')
        except AttributeError:
            # We're probably on Linux. No easy way to get creation dates here,
            # so we'll settle for when its content was last modified.
            tstamp = stat.st_mtime
            if verbose:
                print('Couldn\'t get creation time! Returing modified time.\n' +
                  'In linux. Using os.stat(\'' + filepath + '\').st_mtime as tstamp.')
    return tstamp

def timestamp_from_file(filepath, verbose=True):
    a = re.search('[0-9]{2}h[0-9]{2}', filepath)
    if a is None:
        if verbose:
            print('trying to read creation time')
        timestamp = get_creation_timestamp(filepath)
    else:
        if verbose:
            print('getting timestamp from filename ' + filepath)
        timestamp = timetag_to_timestamp(filepath)
    return timestamp

def load_from_csv(filepath, multiset=False, timestamp=None, verbose=True):
    '''
    This function is made a bit more complicated by the fact that some csvs
    seem to have multiple datasets appended, with a new col_header line as the
    only indication. If multiset=True, this will separate them and return them
    as a list.
    if timestamp = None, the timestamp will be the date created
    I hate that SPEC doesn't save absolute time in a useful way.
    '''
    if verbose:
        print('function \'load_from_csv\' at your service!')
    if timestamp is None:
        a = re.search('[0-9]{2}h[0-9]{2}', filepath)
        if a is None:
            print('trying to read creation time')
            timestamp = get_creation_timestamp(filepath)
        else:
            print('getting timestamp from filename ' + filepath)
            timestamp = timetag_to_timestamp(filepath)

    with open(filepath,'r') as f: # read the file!
        lines = f.readlines()
    colheaders = [col.strip() for col in lines[0].split(',')]
    data = get_empty_set(colheaders, title=filepath,
                         timestamp=timestamp, data_type='SPEC')
    datasets = []

    for line in lines[1:]: #put data in lists!
        vals = [val.strip() for val in line.split(',')]
        not_data = []
        newline = {}
        for col, val in zip(colheaders, vals):
            if col in data['data_cols']:
                try:
                    val = float(val)
                except ValueError:
                    print('value ' + val + ' of col ' + col + ' is not data.')
                    not_data += [col]
            newline[col] = val
        if len(not_data) == len(data['data_cols']):
            print('it looks like there is another data set appended!')
            if multiset:
                print('continuing to next set.')
                numerize(data)
                datasets += [data.copy()]
                colheaders = [val.strip() for val in vals]
                data = get_empty_set(colheaders,
                         timestamp=timestamp, data_type='SPEC')
                continue
            else:
                print('returning first set.')
                numerize(data)
                return data
        else:
            for col in not_data:
                data['data_cols'].remove(col)
                print('column ' + col + ' removed from \'data_cols \'.')

        for col, val in zip(colheaders, vals):
            data[col] += [newline[col]]

    numerize(data)
    datasets += [data]
    if verbose:
        print('function \'load_from_csv\' finished!')
    if multiset:
        return datasets
    return data


def read_macro(file):
    with open(file) as macro:
        lines = macro.readlines()
    lines = remove_comments(lines)
    settings = {'tth':[], 'alpha':[], 'savepath':[], 'newfile':[], 'measurements':[]}
    for line in lines:
        #print(line)
        tth_match = re.search('umv tth ' + float_match, line)
        if tth_match:
            #print('got tth!')
            settings['tth'] += [float(tth_match.group()[7:])]
            continue
        alpha_match = re.search('umv th ' + float_match, line)
        if alpha_match:
            settings['alpha'] += [float(alpha_match.group()[6:])]
            continue
        if 'pd savepath' in line:
            settings['savepath'] += [line[12:]]
            continue
        if 'newfile ' in line:
            settings['newfile'] += [line[8:]]
            continue
        if '_timescan ' in line or 'ascan ' in line or 'pdascan ' in line:
            settings['measurements'] += [line]
            continue
    return settings

def remove_comments(lines):
    new_lines = []
    for line in lines:
        if '#' in line:
            line = re.search('^.*\#', line).group()[:-1]
            if re.search(r'\w', line): #to drop lines that only have comments
                new_lines += [line]
        else:
            new_lines += [line] #I don't want to get rid of empty lines here
    return new_lines



def import_EC_data(full_path_name, title='get_from_file',
                 data_type='EC', N_blank=10, verbose=True,
                 header_string=None, timestamp=None, ):
    file_lines = import_text(full_path_name, verbose=verbose)
    dataset = text_to_data(file_lines, title='get_from_file',
                           data_type='EC', N_blank=10, verbose=True,
                           header_string=None, timestamp=None)
    numerize(dataset)
    return dataset

'''
The following couple functions are adapted from EC_MS on 17L09
as last commited to EC_MS with code c1c6efa
They might benifit from a full rewrite, but not now.
'''


def import_text(full_path_name='current', verbose=True):
    '''
    This method will import the full text of a file selected by user input as a
    list of lines.
    When I first wrote it for EC_MS, way back in the day, I made it so you can
    call it without any arguments, and then input. probably unecessary.
    '''
    if verbose:
        print('\n\nfunction \'import_text\' at your service!\n')

    if full_path_name == 'input':
        full_path_name = input('Enter full path for the file name as \'directory' + os.sep + 'file.extension\'')
    if full_path_name == 'current':
        full_path_name = os.getcwd()

    [directory_name, file_name] = os.path.split(full_path_name)
    original_directory = os.getcwd()
    os.chdir(directory_name)

    if os.path.isdir(full_path_name) and not os.path.isfile(file_name):
        directory_name = full_path_name
        os.chdir(directory_name)
        ls_string = str(os.listdir())
        print('\n' + full_path_name + '\n ls: \n' + ls_string + '\n')
        file_name = input('Directory given. Enter the full name of the file to import\n')

    if verbose:
        print('directory: ' + directory_name)
        print('importing data from ' + file_name )

    possible_encodings = ['utf8','iso8859_15']
    #mpt files seem to be the latter encoding, even though they refer to themselves as ascii
    for encoding_type in possible_encodings:
        try:
            with codecs.open(file_name, 'r', encoding = encoding_type) as file_object:
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
        print('\nfunction \'import_text\' finished!\n\n')
    return file_lines


def text_to_data(file_lines, title='get_from_file',
                 timestamp=None, date='today', tstamp=None, tz=None,
                 data_type='EC', N_blank=10, sep=None,
                 header_string=None, verbose=True):
    '''
    This method will organize data in the lines of text from a file useful for
    electropy into a dictionary as follows (plus a few more keys)
    {'title':title, 'header':header, 'timestamp':timestamp,
     'data_cols':[colheader1, colheader2, ...],
     colheader1:[data1], colheader2:[data2]...}
     So far made to work with SPEC files (.csv), XAS files (.dat),
     EC_Lab files (.mpt).
     If you need to import cinfdata files (.txt), you are in the wrong place.
     This is EC_Xray's import_data, not EC_MS!!!
    '''
    if verbose:
        print('\n\nfunction \'text_to_data\' at your service!\n')

    from .Combining import get_timecol

    #disect header
    N_lines = len(file_lines)        #number of header lines
    N_head = N_lines                 #this will change when I find the line that tells me how ling the header is
    header_string = ''

    dataset = {}
    commacols = []                   #will catch if data is recorded with commas as decimals.

    loop = False

    if data_type == 'SPEC':
        N_head = 1 #the column headers are the first line
        if sep is None:
            sep = ','
    elif data_type == 'XAS':
        datacollines = False # column headers are on multiple lines, this will be True for those lines
        col_headers = []   #this is the most natural place to initiate the vector

    if sep is None:  #EC, XAS, and MS data all work with '\t'
        sep = '\t'

    n_blank = 0
    for nl, line in enumerate(file_lines):
        l = line.strip()
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
                elif timestamp is None and re.search('Acquisition started',line):
                    timestamp_object = re.search(timestamp_match,l)
                    timestamp = timestamp_object.group()
                    date_object = re.search(date_match, l)
                    date = date_object.group()
                    if verbose:
                        print('timestamp \'' + timestamp + '\' found in line ' + str(nl))
                elif re.search('Number of loops', line): #Then I want to add a loop number variable to data_cols
                    loop = True
                    dataset['loop number'] = []
                elif re.search('Loop', line):
                    n = int(re.search(r'^Loop \d+', line).group()[5:])
                    start = int(re.search(r'number \d+', line).group()[7:])
                    finish = int(re.search(r'to \d+', line).group()[3:])
                    N = finish - start + 1
                    dataset['loop number'] += N * [n]

                header_string = header_string + line


            elif data_type == 'XAS':
                if datacollines:
                    if l == '': #then we're done geting the column headers
                        datacollines = False
                        #header = False  #and ready for data.
                        N_head = nl + 1 # so that the next line is data!
                        dataset['data_cols'] = col_headers.copy()
                        if verbose:
                            print('data col lines finish on line' + str(nl) +
                                  '. the next line should be data')
                    else:
                        col_headers += [l]
                        dataset[l] = []
                elif l == 'Data:':  #then we're ready to get the column headers
                    datacollines = True
                    if verbose:
                        print('data col lines start on line ' + str(nl+1))
                a = re.search(timestamp_match, l)
                if timestamp is None and a is not None:
                    timestamp = a.group()
                    d = re.search(date_match, l)
                    if d is not None:
                        date = d.group()
                    tstamp = timestamp_to_epoch_time(l, tz=tz, verbose=verbose) #the XAS data is saved with time.ctime()
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
                        timestamp_object = re.search(timestamp_match, string2.strip())
                        timestamp = timestamp_object.group()
                        date_object = re.search(date_match_2, string2.strip())
                        date = date_object.group()
                        date = date[5:7] + '/' + date[-2:] + '/' + date[:4]
                        # ^convert yyyy-mm-dd to dd-mm-yyyy
                        if verbose:
                            print('timestamp \'' + timestamp + '\' found in line ' + str(nl))
                header_string = header_string + line

        elif nl == N_head - 1:      #then it is the column-header line
               #(EC-lab includes the column-header line in header lines)
            #col_header_line = line
            col_headers = [col.strip() for col in l.split(sep=sep)]
            dataset['N_col']=len(col_headers)
            dataset['data_cols'] = col_headers.copy()  #will store names of columns containing data
            #DataDict['data_cols'] = col_headers.copy()
            for col in col_headers:
                dataset[col] = []              #data will go here
            header_string = header_string + line #include this line in the header
            if verbose:
                print('Data starting on line ' + str(N_head) + '\n')


        else:                   # data, baby!
            line_data = [dat.strip() for dat in l.split(sep=sep)]
            if not len(line_data) == len(col_headers):
                #print('Mismatch between col_headers and data on line ' + str(nl) + ' of ' + title) #debugging
                pass
            for col, data in zip(col_headers, line_data):
                if col in dataset['data_cols']:  #why would it not?
                    if get_timecol(col) not in col_headers:
                        print('Missing time column ' + get_timecol(col) +
                              ' on line ' + str(nl) + '. Dropping ' + col)
                    try:
                        data = float(data)
                    except ValueError:
                        if data == '':
                            continue        #added 17C22 to deal with data acquisition crashes.
                        try:
                            if verbose and not col in commacols:
                                print('ValueError on value ' + data + ' in column ' + col + ' line ' + str(nl) +
                                      '\n Checking if you''re using commas as decimals in that column... ')
                            data = data.replace('.','')
                            # ^ in case there's also '.' as thousands separator, just get rid of it.
                            data = data.replace(',','.')       #put '.' as decimals
                            data = float(data)
                            if not col in commacols:
                                if verbose:
                                    print('... and you were, dumbass. I''ll fix it.')
                                commacols += [col]
                        except ValueError:
                            if verbose :
                                print(list(zip(col_headers,line_data)))
                                print(title + ' in text_to_data: \nRemoved \'' + str(col) +
                                      '\' from data columns because of value \'' +
                                    str(data) + '\' at line ' + str(nl) +'\n')
                            dataset['data_cols'].remove(col)

                dataset[col].append(data)

    if loop:
        dataset['data_cols'] += ['loop number']
    dataset['title'] = title
    dataset['header'] = header_string
    dataset['timestamp'] = timestamp
    dataset['date'] = date
    if tstamp is None:
        tstamp = timestamp_to_epoch_time(timestamp, date, tz=tz, verbose=verbose)
    dataset['timezone'] = tz
    dataset['tstamp'] = tstamp
    #UNIX epoch time, for proper synchronization!
    dataset['data_type'] = data_type

    if data_type == 'EC':           #so that synchronize can combine current data from different EC-lab techniques
        if '<I>/mA' in dataset['data_cols'] and 'I/mA' not in dataset['data_cols']:
            dataset['data_cols'].append('I/mA')
            dataset['I/mA'] = dataset['<I>/mA'] #so that synchronize can combine current data from different EC-lab techniques
        if '<Ewe>/V' in dataset['data_cols'] and 'Ewe/V' not in dataset['data_cols']:
            dataset['data_cols'].append('Ewe/V')
            dataset['Ewe/V'] = dataset['<Ewe>/V']

    if verbose:
        print('\nfunction \'text_to_data\' finished!\n\n')
    return dataset


def import_data(*args, **kwargs):
    print('\'import_data\' is now called \'load_from_file\'!\n' +
          'Remember that next time, goof.')
    return load_from_file(*args, **kwargs)

def load_from_file(full_path_name='current', title='file', tstamp=None, timestamp=None,
                 data_type='EC', N_blank=10, tz=None, verbose=True):
    '''
    This method will organize the data in a file useful for
    electropy into a dictionary as follows (plus a few more keys)
    {'title':title, 'header':header, 'timestamp':timestamp,
    'data_cols':[colheader1, colheader2, ...],
    colheader1:[data1], colheader2:[data2]...}
    So far made to work with SPEC files (.csv), XAS files (.dat),
    EC_Lab files (.mpt).
    If you need to import cinfdata files (.txt), you are in the wrong place.
    Use EC_MS's import_data instead!!!
    '''
    if verbose:
        print('\n\nfunction \'load_from_file\' at your service!\n')
    if title == 'file':
        folder, title = os.path.split(full_path_name)
    file_lines = import_text(full_path_name, verbose)
    dataset = text_to_data(file_lines=file_lines, title=title, data_type=data_type,
                           timestamp=timestamp, N_blank=N_blank, tz=tz, tstamp=tstamp,
                           verbose=verbose)
    if tstamp is not None: #then it overrides whatever test_to_data came up with.
        dataset['tstamp'] = tstamp
    elif dataset['tstamp'] is None:
        dataset['tstamp'] = get_creation_time(full_path_name, verbose=verbose)
    numerize(dataset)

    if verbose:
        print('\nfunction \'load_from_file\' finished!\n\n')
    return dataset


def import_EC_set(directory, EC_file=None, tag='01',
                  verbose=True, tz=None):
    if verbose:
        print('\n\nfunction \'load_EC_set\' at your service!\n')
    from .Combining import synchronize, sort_time

    lslist = os.listdir(directory)

    if EC_file is None:
        EC_file = [f for f in lslist if f[:2] == tag and f[-4:] == '.mpt']
    elif type(EC_file) is str:
        EC_file = [EC_file]
    EC_datas = []
    for f in EC_file:
        try:
            EC_datas += [load_from_file(directory + os.sep + f, data_type='EC', tz=tz, verbose=verbose)]
        except OSError:
            print('problem with ' + f + '. Continuing.')
    EC_data = synchronize(EC_datas, verbose=verbose, append=True, t_zero='first', tz=tz)
    if 'loop number' in EC_data['data_cols']:
        sort_time(EC_data, verbose=verbose) #note, sort_time no longer returns!

    if verbose:
         print('\nfunction \'load_EC_set\' finished!\n\n')
    return EC_data

def load_EC_set(*args, **kwargs):
    print('load_EC_set has been renamed import_EC_set')
    return import_EC_set(*args, **kwargs)


def download_cinfdata_set(group_id=None, grouping_column=None, **kwargs):

    if grouping_column is None:
        grouping_column, group_id = kwargs.popitem()

    from .Combining import synchronize

    try:
        from cinfdata import Cinfdata
    except ImportError:
        print('the cinfdata module must be on your python path. It\'s here: \n' +
              'https://github.com/CINF/cinf_database/blob/master/cinfdata.py')

    try:
        cinfd = Cinfdata('sniffer', grouping_column=grouping_column,
                         allow_wildcards=True,
                         label_column='mass_label')
    except:
        raise  #untill I know exactly which error I'm trying to catch.
        print('couldn\'t connect. You should run gstm')
        #os.system('gstm')
        raise RuntimeError('Couldn\'t connect to cinfdata!')

    #obj = cinfd.get_metadata_group('2018-03-30 14:13:17')

    obj = cinfd.get_metadata_group(group_id)
    idlists = {} # keys will be time as string. values will be corresponding id's

    for key, value in obj.items():
        #label = value['mass_label']
        #print(label)
        timestamp = str(value['time'])
        if timestamp not in idlists:
            idlists[timestamp] = []
        idlists[timestamp] += [value['id']]

    datasets = {}
    for timestamp, idlist in idlists.items():

        if len(idlist) == 0:
            print('No data associated with timestamp \'' + timestamp + '\'.')
            continue

        dataset = {'title':timestamp, 'data_type':'MS'}

        metadatas = dict([(i, cinfd.get_metadata(i)) for i in idlist])

        unixtimes = [metadatas[i]['unixtime'] for i in idlist]
        if len(set(unixtimes)) > 1:
            msg = 'unix times don\'t match for timestamp \'' + timestamp + '\'!'
            raise ValueError(msg)

        dataset['tstamp'] = unixtimes[0]
        dataset['timestamp'] = metadatas[idlist[0]]['time'].strftime('%H:%M:%S')

        labels = [metadatas[i]['mass_label'] for i in idlist]
        if 'Mass Scan' in labels:
            dataset['scan_type'] = 'mass'
        else:
            dataset['scan_type'] = 'time'

        dataset['data_cols'] = []
        dataset['timecols'] = []
        for i in idlist: #avoiding id since it's got a builtin meaning
            data = cinfd.get_data(i)
            label = metadatas[i]['mass_label']
            if len(data.shape) == 1:
                dataset[label] = data
                dataset['data_cols'] += [label]
            elif data.shape[1] == 2:
                x = data[:, 0]
                y = data[:, 1]
                x_label = label + '-x'
                y_label = label + '-y'
                dataset['timecols'] += [(x_label, y_label)]
                dataset[x_label] = x * 1e-3 # cinfdata saves time in ms!!!
                dataset[y_label] = y

                dataset['data_cols'] += [x_label, y_label]

        datasets[timestamp] = dataset


    timescans = [dataset for dataset in datasets.values() if dataset['scan_type'] == 'time']


    combined = synchronize(timescans, t_zero='first')

    return combined


def get_xy(data, xcol=None, ycol=None, label=None, ):
    '''
    '''
    if xcol is None:
        xcol = label + '-x'
    if ycol is None:
        ycol = label + '-y'
    return data[xcol], data[ycol]

def set_xy(data, x, y, xcol=None, ycol=None, label=None, ):
    '''
    '''
    if xcol is None:
        xcol = label + '-x'
    if ycol is None:
        ycol = label + '-y'
    data[xcol] = x
    data[ycol] = y

def remove_repeats(data, xcol=None, ycol=None, label=None):
    '''
    '''
    x_0, y_0 = get_xy(data, xcol, ycol, label)
    x, y = only_while_increasing(x_0, y_0)
    set_xy(data, x=x, y=y, xcol=xcol, ycol=ycol, label=label)

def only_while_increasing(x=None, y=None):
    '''
    removes the repeats if a dataset goes back and repeats, as happens in
    the Analog In anomoly first observed 18D05.
    Does so in a vectorized way (only loops over "cliff points" where x falls)
    x is monotonically increasing in the returned data
    '''

    x_up = np.append(x[1:], x[-1]+1) # x shifted up one, so that x_up[i] = x[i+1]
    cliff_points = np.where(x_up < x)[0]      # points right before a drop in x
    mask = np.tile(True, np.size(x))

    for point in cliff_points:
        x_cliff = x[point]
        mask[point:] = np.logical_and(mask[point:], x[point:] > x_cliff)

    return x[mask], y[mask]


def import_set(directory, MS_file='QMS.txt', MS_data=None, t_zero='start',
               EC_file=None, tag='01',
               cutit=False, cut_buffer=60,
               verbose=True, override=False):
    from .Combining import synchronize, sort_time

    if verbose:
        print('\n\nfunction import_set at your service!\n')


    lslist = os.listdir(directory)

    if MS_data is None:
        if type(MS_file) is str:
            MS_file = [MS_file]
        MS_datas = [load_from_file(directory + os.sep + f,
                                data_type='MS', verbose=verbose)
                    for f in MS_file]
        MS_data = synchronize(MS_datas, verbose=verbose)
        if len(MS_datas) > 1:
            sort_time(MS_data)

    if EC_file is None:
        EC_file = [f for f in lslist if f[:2] == tag and f[-4:] == '.mpt']
    elif type(EC_file) is str:
        EC_file = [EC_file]
    EC_datas = [load_from_file(directory + os.sep + f, verbose=verbose, data_type='EC')
                for f in EC_file]
    EC_data = synchronize(EC_datas, verbose=verbose)
    if 'loop number' in EC_data['data_cols']:
        sort_time(EC_data, verbose=verbose) #note, sort_time no longer returns!

    data = synchronize([MS_data, EC_data], t_zero=t_zero, verbose=verbose,
                       override=override, cutit=cutit, cut_buffer=cut_buffer)
    if verbose:
         print('\nfunction import_set finished!\n\n')
    return data

def save_as_text(filename, dataset, cols=[], mols=[], tspan='all', header=None,
                 N_chars=None, timecols={}, **kwargs):
    '''
    kwargs is fed directly to Molecule.get_flux()
    '''
    from .Combining import get_timecol, cut

    lines = []
    if type(header) is list:
        lines += header
    elif type(header) is str:
        lines += [header]

    if N_chars is None:
        N_chars = max([len(col) for col in cols])

    col_header = ''
    i_col = 0
    columns = []
    datas = {}
    for col in cols:
        if col in timecols:
            tcol = timecols[col]
        else:
            tcol = get_timecol(col)
        if tcol in dataset and tcol not in columns: # don't want same tcol twice
            col_header += ('{0:>' + str(N_chars) + 's},\t').format(tcol)
            columns += [tcol]
            i_col += 1
        if col in dataset and col not in columns: # don't want same tcol twice
            col_header += ('{0:>' + str(N_chars) + 's},\t').format(col)
            columns += [col]
            i_col += 1
        else:
            print(col + ' not in dataset. ignoring it.')
            continue
        if tcol in columns:
            x, y = dataset[tcol].copy(), dataset[col].copy()
            if tspan is not False and not tspan=='all':
                x, y = cut(x, y, tspan=tspan)
            datas[tcol], datas[col] = x, y
        else:
            print('timecol \'' + tcol + '\' for col \'' + col +
                  '\' is not in dataset, so can\'t cut it.')
            datas[col] = dataset[col].copy()

    for mol in mols:
        tcol = mol.name + '_' + mol.primary + '-x'
        col = mol.name + '_' + mol.primary + '-y'
        x, y = mol.get_flux(dataset, tspan=tspan, **kwargs)
        datas[tcol] = x
        datas[col] = y
        col_header += ('{0:>' + str(N_chars) + 's},\t').format(tcol)
        columns += [tcol]
        col_header += ('{0:>' + str(N_chars) + 's},\t').format(col)
        columns += [col]

    lines += [col_header + '\n']

    i_data = 0
    finished = False
    while not finished:
        N_unfinished = 0
        line = ''
        for col in columns:
            try:
                d = datas[col][i_data]
                line += ('{0:>' + str(N_chars) + '.6g},\t').format(d)
                N_unfinished += 1
            except IndexError:
                line += ' '*N_chars + ',\t'
        if N_unfinished == 0:
            finished = True
        else:
            lines += [line + '\n']
            i_data += 1

    with open(filename, 'w') as f:
        f.writelines(lines)


