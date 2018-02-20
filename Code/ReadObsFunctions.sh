def read_daily_precip(file_name, file_colnames=None, header='infer', delimiter='\s+'):
    # read in a daily precipitation data set
    
    # if file_colnames are supplied, use header=None
    if ps.notnull(file_colnames):
        header=None
    
    # read in the data
    daily_data=pd.read_table(file_name, delimiter=delimiter, header=header) 
    
    # set columns, if header=None
    if pd.notnull(file_colnames):
        daily_data.columns=file_colnames
    else:
        file_colnames=list(daily_data.columns)
    
    # calculate cfs to cms conversion, or vice versa
    if 'precip_m' in daily_data.columns:
        precip_m=daily_data['precip_m']
        precip_mm=precip_m*1000
    
    # determine the datetime
    date_index=[file_colnames.index(each) for each in ['year','month','day']]
    row_dates=pd.to_datetime(daily_data[date_index])
    
    # generate the daily_flow and set the datetime as row indices
    daily_precip=pd.concat([precip_m, precip_mm],axis=1)
    daily_precip.set_index(row_dates, inplace=True)
    daily_precip.columns=['precip_m', 'precip_mm']
    return(daily_precip)


def read_daily_snotel(file_name, file_colnames=None, usecols=None, delimiter=',', header='infer'):
    # read in a daily SNOTEL observation data set
    
    # if file_colnames are supplied, use header=None
    if file_colnames is not None:
        header=None
    
    # read in the data
    daily_data=pd.read_table(file_name, usecols=usecols, header=header, delimiter=delimiter)
    
    # reset the colnames
    daily_data.columns=['Date', 'Tmax_C', 'Tmin_C', 'Tavg_C', 'Precip_mm']
    
    # transform the data
    daily_data['Tmax_C']=(daily_data['Tmax_C'] -32)/1.8
    daily_data['Tmin_C']=(daily_data['Tmin_C'] -32)/1.8
    daily_data['Tavg_C']=(daily_data['Tavg_C'] -32)/1.8
    daily_data['Precip_mm']=daily_data['Precip_mm'] *25.4
    
    # determine the datetime
    row_dates=pd.to_datetime(daily_data.Date)
    
    # generate the daily_flow and set the datetime as row indices
    daily_snotel=daily_data[['Tmax_C', 'Tmin_C', 'Tavg_C', 'Precip_mm']]
    daily_snotel.set_index(row_dates, inplace=True)
    return(daily_snotel)


def read_daily_coop(file_name, file_colnames=None, usecols=None, delimiter=',', header='infer'):
    # read in a daily COOP observation data set
    
    # if file_colnames are supplied, use header=None
    if file_colnames is not None:
        header=None
    
    # read in the data
    daily_data=pd.read_table(file_name, usecols=usecols, header=header, delimiter=delimiter, 
                             date_parser=lambda x: pd.datetime.strptime(x, '%Y%m%d'), 
                             parse_dates=[0], 
                             na_values=-9999)
    
    # reset the colnames
    daily_data.columns=['Date', 'Precip_mm','Tmax_C', 'Tmin_C', 'Tavg_C']
    
    # transform the data
    daily_data['Tmax_C']=(daily_data['Tmax_C'] -32)/1.8
    daily_data['Tmin_C']=(daily_data['Tmin_C'] -32)/1.8
    daily_data['Tavg_C']=(daily_data['Tavg_C'] -32)/1.8
    daily_data['Precip_mm']=daily_data['Precip_mm'] *25.4
    
    # determine the datetime
    row_dates=pd.to_datetime(daily_data.Date)
    
    # generate the daily_flow and set the datetime as row indices
    daily_coop=daily_data[['Precip_mm','Tmax_C', 'Tmin_C', 'Tavg_C']]
    daily_coop.set_index(row_dates, inplace=True)
    return(daily_coop)
