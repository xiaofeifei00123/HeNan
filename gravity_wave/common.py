class Common():
    def __init__(self, ):

        self.areaA = {
            'lat1':34.4,
            'lat2':34.9,
            'lon1':113.0,
            'lon2':113.7,
        }        
        self.areaB = {
            'lat1':35.3,
            'lat2':35.8,
            'lon1':113.5,
            'lon2':114.2,
        }        

        self.station = {
                'A': {
                    'abbreviation':'A',
                    'lat': 33.6,
                    'lon': 112.1,
                },
                'B': {
                    'abbreviation':'B',
                    'lat': 34.5,
                    'lon': 112.9,
                },
                'C': {
                    'abbreviation':'C',
                    'lat': 35.15,
                    'lon': 113.5,
                },
                'D': {
                    'abbreviation':'D',
                    'lat': 35.75,
                    'lon': 114.0,
                },
            }
        self.cross_start = [111, 34.5]
        self.cross_end = [114.5, 32.5]