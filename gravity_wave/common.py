class Common():
    def __init__(self, ):

        self.areaA = {
            'lat1':33.5,
            'lat2':33.8,
            'lon1':112.4,
            'lon2':112.8,
        }        
        self.areaB = {
            'lat1':34.4,
            'lat2':34.9,
            'lon1':113.0,
            'lon2':113.7,
        }        
        self.areaC = {
            'lat1':35.3,
            'lat2':35.8,
            'lon1':113.5,
            'lon2':114.2,
        }        
        ## 模式降水核心区
        self.areaD = {
            'lat1':34.4,
            'lat2':34.8,
            'lon1':113.0,
            'lon2':113.4,
        }        
        ## 观测降水核心区
        self.areaE = {
            'lat1':34.4,
            'lat2':34.8,
            'lon1':113.2,
            'lon2':113.6,
        }        

        self.station_line = {
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

        self.station_sta = {
                'ZhengZhou': {
                    'abbreviation':'郑州',
                    'lat': 34.76,
                    'lon': 113.65
                },
                'NanYang': {
                    'abbreviation':'南阳',
                    'lat': 33.1,
                    'lon': 112.49,
                },
                'LuShi': {
                    'abbreviation':'卢氏',
                    'lat': 34.08,
                    'lon': 111.07,
                },
            }

        self.station_zz = {
                'ZhengZhou': {
                    'abbreviation':'Zhengzhou',
                    'lat': 34.76,
                    'lon': 113.65
                },
            }
        self.cross_start = [111, 34.5]
        self.cross_end = [114.5, 32.5]