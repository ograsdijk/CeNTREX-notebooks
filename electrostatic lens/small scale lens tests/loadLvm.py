import pandas as pd
import linecache
import datetime as dt


def loadlvm_electrostaticlenstests(fname):
    data = pd.read_table(fname, header=22)
    data = data.iloc[:, [0, 1, 2, 3, 4]]
    data.columns = ['timestamp', 'V1', 'V2', 'I1', 'I2']
    units = ['Seconds'] + linecache.getline(fname, 18).replace('\n', '').split('\t')[1:-1]
    date = linecache.getline(fname, 10).replace('\n', '').split('\t')[1]
    time = linecache.getline(fname, 11).replace('\n', '').split('\t')[1]
    starttime = dt.datetime.strptime(date + ' ' + time[:8], '%Y/%m/%d %H:%M:%S')
    data['time'] = dt.timedelta(seconds=1) * data.timestamp.values + starttime
    return data
