import requests

with open('apikey.txt', 'r') as file_key:
    key = file_key.readline()

stationIDs = ['06180']
parameterIDs = ['temp_mean_past1h', 'wind_speed_past1h'] 
# Available Parameter Values: https://confluence.govcloud.dk/display/FDAPI/Meteorological+Observation+Data#MeteorologicalObservationData-Parameters

MAX_LIMIT = 300000 # Maximum number of datapoints allowed per API-call

for stationID in stationIDs:

    for parameterID in parameterIDs:
    
        # This is NOT the proper way to do get-requests in python with requests, but quicker to implement perhaps
        url = 'https://dmigw.govcloud.dk/v2/metObs/collections/observation/items?'
        url += 'stationId=' + stationID
        url += '&datetime=2018-01-01T00:00:00Z/2023-11-28T00:00:00Z'
        url += '&parameterId=' + parameterID
        url += '&limit=' + str(MAX_LIMIT)
        url += '&api-key=' + key

        r = requests.get(url=url)

        data = r.json()
    
        file_name = 'station=' + stationID + '_param=' + parameterID + '.csv'

        with open(file_name, 'w') as file_data:
            # Write Header
            file_data.write('datetime, ' + parameterID + '\n')
            
            # Write Data Values
            for idx in range(len(data['features'])-1, -1, -1): 
                file_data.write( data['features'][idx]['properties']['observed'] + ' , ' + str(data['features'][idx]['properties']['value']) + '\n' )