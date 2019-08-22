import pandas as pd
import os
import datetime as dt


def load_radeyeg20(fname):
    data = pd.read_table(fname, header=3, encoding='cp1250')
    data = data.rename(columns={'Unit': 'Unit Dose Rate', 'Unit.1': 'Unit Dose'})
    d = data['mm/dd/yyyy'] + ' ' + data['hh:mm:ss']
    data['time'] = pd.to_datetime(d, infer_datetime_format=True)
    if os.name != 'nt':
        timecreated = dt.datetime.fromtimestamp(os.stat(fname).st_birthtime)
        if timecreated < dt.datetime(2018, 2, 21):
            data.time -= (data.time.iloc[-1] - data.time.iloc[0])
            data.time -= data.time.iloc[0]
            data.time += timecreated - (data.time.iloc[-1] - data.time.iloc[0])
    return data
