from csv import writer
import os.path


def wavelet_data_saver(SAVE_PATH, csvname,errorbar_packed,t):
    solarmax=errorbar_packed[0]
    solarmin=errorbar_packed[1]
    max_power_period=errorbar_packed[2]
    max_power_min=errorbar_packed[3]
    max_power_max=errorbar_packed[4]
    index_max_power=errorbar_packed[5]
    solarmax_period=errorbar_packed[6]
    solarmax_period_min=errorbar_packed[7]
    solarmax_period_max=errorbar_packed[8]
    index_solarmax=errorbar_packed[9]

    print("******************************************************************************")  
    print('Solarmax: ',solarmax, ', Solarmin: ',solarmin)
    print("Period with max power: %.2f d"%max_power_period)
    print('Min Period at CI : %.2f d'%max_power_min)
    print('Max Period at CI : %.2f d'%max_power_max)
    print("******************************************************************************")

    Header=[
        '#Solarmax',
        '#Solarmin',
        '#Period at Max Power',
        '#Min Period at Max Power',
        '#Max Period at Max Power',
        '#Time at Max Power',
        '#Period at Solar Max',
        '#Min Period at Solarmax',
        '#Max Period at Solarmax'
    ]

    data_list = [solarmax, solarmin, max_power_period, max_power_min, max_power_max, t[index_max_power], solarmax_period, solarmax_period_min, solarmax_period_max, index_solarmax]
    file_exists = os.path.isfile(SAVE_PATH+'/'+csvname+'.csv')
    """
    with open(SAVE_PATH+'/'+csvname+'.csv', 'a') as f_object:
        if not file_exists:
            writer_object = writer(f_object)
            writer_object.writerow([Header])

        writer_object = writer(f_object)
        writer_object.writerow(data_list)
        f_object.close()
    print('CSV output saved to: ',SAVE_PATH+'/'+csvname+'.csv')
    """
    print('Analysis Complete')
    print("******************************************************************************")
    return solarmax, solarmin, max_power_period, max_power_min, max_power_max, t[index_max_power], solarmax_period, solarmax_period_min, solarmax_period_max, index_solarmax
