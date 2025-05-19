import ee
import math

DEG2RAD = math.pi / 180.0

def meteo_era5land(time_start, meteorology_source_daily):
    """
    Parameters
    ----------
    time_start : str
        Image property: time start of the image.
    meteorology_source_inst: ee.ImageCollection, str
        Instantaneous meteorological data.
    meteorology_source_daily :  ee.ImageCollection, str
        Daily meteorological data.

    Returns
    -------
    ee.Image

    Notes
    -----
    Accepted collections:
    Inst : ECMWF/ERA5_LAND/HOURLY
    Daily : projects/openet/assets/meteorology/era5land/na/daily
            projects/openet/assets/meteorology/era5land/sa/daily

    References
    ----------

    """

    # Get date information
    time_start = ee.Number(time_start)

    # Filtering Daily data
    meteorology_daily = (
        ee.ImageCollection(meteorology_source_daily)
        .filterDate(ee.Date(time_start).advance(-1, "day"), ee.Date(time_start).advance(1, "day"))
        .first()
    )

    # Incoming shorwave down [W m-2]
    swdown24h = meteorology_daily.select("surface_solar_radiation_downwards_sum") #.divide(1 * 60 * 60 * 24)

    swdown24h = (
        ee.ImageCollection('ECMWF/ERA5_LAND/HOURLY')
        .filterDate(ee.Date(time_start), ee.Date(time_start).advance(1, "day"))
        .select("surface_solar_radiation_downwards_hourly")
        .sum()
        .divide(1 * 60 * 60*24)
        .rename("rso_inst")
    )

    # Minimum air temperature [K]
    tmin = meteorology_daily.select("temperature_2m_min").rename("tmin")
    
    # Maximum air temperature [K]
    tmax = meteorology_daily.select("temperature_2m_max").rename("tmax")

    # Air temperature [C]
    tair_k = meteorology_daily.select('temperature_2m').rename('tair')

    # Wind speed [ m/s]
    wind_u = meteorology_daily.select('v_component_of_wind_10m')

    wind_v = meteorology_daily.select('u_component_of_wind_10m')

    wind_med = wind_u.expression(
        "sqrt(ux_u ** 2 + ux_v ** 2)",
        {"ux_u": wind_u, "ux_v": wind_v},
    ).rename("ux")

    wind_med = wind_med.expression("ux * (4.87) / log(67.8 * z - 5.42)", {"ux": wind_med, "z": 10.0}).rename("ux")

    # Dew point temperature [°K]
    tdp = meteorology_daily.select('dewpoint_temperature_2m').rename("tdp")

    # Actual vapour pressure [kPa]
    ea = tdp.expression("0.6108 * (exp((17.27 * T_air) / (T_air + 237.3)))", {"T_air": tdp.subtract(273.15)})

    # Saturated vapour pressure [kPa]
    esat = tair_k.expression("0.6108 * (exp((17.27 * T_air) / (T_air + 237.3)))", {"T_air": tair_k.subtract(273.15)})

    # Relative humidity(%)
    rh = ea.divide(esat).multiply(100).rename("RH")

    tmin = tmin.subtract(273.15).resample("bilinear")
    tmax = tmax.subtract(273.15).resample("bilinear")
    tair_k = tair_k.resample("bilinear")
    wind_med = wind_med.resample("bilinear")
    rh = rh.resample("bilinear")
    swdown24h = swdown24h.resample("bilinear")

    return [tmin, tmax, tair_k, wind_med, rh, swdown24h]



def radiation_24h(time_start, tmax, tmin, dem, rso24h):
    """
    Daily Net radiation [W m-2] - FAO56

    Parameters
    ----------
    time_start : ee.Date
        Date information of the image.
    tmax : ee.Image
        Maximum air temperature [Celsius].
    tmin : ee.Image
        Minimum air temperature [Celsius].
    elev : ee.Image
        Digital Elevation information [m].
    sun_elevation : ee.Number, int
        Sun elevation information.
    cos_terrain : ee.Image
        Solar zenith angle cos (aspect/slope).
    rso24h : ee.Image
        Daily Short wave radiation [W m-2]

    Returns
    -------
    ee.Image

    References
    ----------
    .. [FAO56] Allen, R., Pereira, L., Raes, D., & Smith, M. (1998).
       Crop evapotranspiration: Guidelines for computing crop water
       requirements. FAO Irrigation and Drainage Paper (Vol. 56).

    """

    # Convert to MJ m-2
    rs = rso24h.multiply(0.0864).rename("Rs")

    # Solar constant [MJ m-2]
    gsc = 0.0820

    # Day of the year
    doy = ee.Date(time_start).getRelative("day", "year").add(1)

    # Inverse relative distance earth-sun (FAO56 Eqn 23)
    dr = tmax.expression("1 + (0.033 * cos((2 * pi / 365) * doy))", {"doy": doy, "pi": math.pi})

    # Solar declination [rad] (FAO56 Eqn 24)
    sd = tmax.expression("0.40928 * sin(((2 * pi / 365) * doy) - 1.39)", {"doy": doy, "pi": math.pi})

    # Latitude of the image
    lat = tmax.pixelLonLat().select(["latitude"]).multiply(DEG2RAD).rename("latitude")

    #  Sunset hour angle [rad] (FAO56 Eqn 25)
    ws = tmax.expression("acos(-tan(Lat) * tan(Sd))", {"Lat": lat, "Sd": sd})

    # Extraterrestrial radiation [MJ m-2 d-1] (FAO56 Eqn 21)
    rad_a = tmax.expression("Ws * sin(Lat) * sin(Sd) + cos(Lat) * cos(Sd) * sin(Ws)", {"Ws": ws, "Lat": lat, "Sd": sd})

    ra = tmax.expression("((24 * 60) / pi) * Gsc * Dr * rad_a", {"pi": math.pi, "Gsc": gsc, "Dr": dr, "rad_a": rad_a})
    
    # Simplified clear sky solar formulation [MJ m-2 d-1] (FAO56 Eqn 37)
    rso = tmax.expression("(0.75 + 2E-5 * z) * Ra", {"z": dem, "Ra": ra})

    # Net shortwave radiation [MJ m-2 d-1] (FAO56 Eqn 38)
    rns = tmax.expression("(1 - albedo) * Rs", {"Rs": rs, "albedo": 0.23})

    # Actual vapor pressure [MJ m-2 d-1] (FAO56 Eqn 11)
    ea_max = tmax.expression(
     ' 0.6108 *(exp( (17.27 * tmax) / (tmax + 237.3)))', {
       'tmax': tmax,  
       }).rename('ea_max')
    
    ea_min = tmax.expression(
     ' 0.6108 *(exp( (17.27 * tmin) / (tmin + 237.3)))', {
       'tmin': tmin,  
       }).rename('ea_min')
    

    es_mean = tmax.expression(
     '(ea_min +ea_max)/2', {
       'ea_min': ea_min,  
       'ea_max':ea_max
       }).rename('es_mean') 

    # Cloudness function (dimensionless)
    fcd = tmax.expression(
        "(1.35 * (rs / rso) - 0.35)",{
        'rs':rs,
        'rso':rso
    })

    # Net longwave radiation [MJ m-2 d-1] (FAO56 Eqn 39)
    rnl = tmax.expression(
        "4.901E-9 * ((tmax ** 4 + tmin ** 4) / 2) * (0.34 - 0.14 * sqrt(ea)) * fcd",
        {"tmax": tmax.add(273.15), 
         "tmin": tmin.add(273.15), "ea": es_mean, "fcd":fcd.clamp(0.05,1)},
    )

    # Net radiation [MJ m-2 d-1] (FAO56 Eqn 40)
    rn = tmax.expression("rnd - rnl", {"rnd": rns, "rnl": rnl})

    return rn.rename("rad_24h")


def etr(tmin, tmax, tair, ws, dem, rad_24h):
    
    """
    Reference Evapotranspiration - ASCE 2005

    Parameters
    ----------
    tmin : ee.Image
        Minimum air temperature [Celsius]
    tmax : ee.Image
        Maximum air temperature [Celsius].
    tair : ee.Image
        Air temperature [K].
    ws : ee.Image
        Wind speed [m s-1]
    dem : ee.Image
        Digital Elevation information [m].
    rad_24h : ee.Image
        Net daily radiation [MJ m-2 d-1]

    Returns
    -------
    ee.Image

    References
    ----------
    .. [ASCE 2005] Task Committee on Standardization of Reference Evapotranspiration
        Environmental and Water Resources Institute of the American Society of Civil Engineers. 
        January 2005.

    """

    # Atm pressure (kPa)
    press = dem.expression(
    '101.3 * pow(((293 - (0.0065 * dem))/ 293),5.26) ', {
       'dem' : dem,
       }).rename('p_atm')
    
    # psychrometric constant (kPa °C-1)
    const_psci=dem.expression(
            '0.000665* press ',{
                    'press':press}).rename('cte_psi')
    
    # max actual vapor pressure (kPa)  
    ea_max = dem.expression(
     ' 0.6108 *(exp( (17.27 * tmax) / (tmax + 237.3)))', {
       'tmax': tmax,  
       }).rename('ea_max')
    
    # min actual vapor pressure (kPa)
    ea_min = dem.expression(
     ' 0.6108 *(exp( (17.27 * tmin) / (tmin + 237.3)))', {
       'tmin': tmin,  
       }).rename('ea_min')
    
    # saturation vapor pressure (kPa)
    es_mean = dem.expression(
     '(ea_min +ea_max)/2', {
       'ea_min': ea_min,  
       'ea_max':ea_max
       }).rename('es_mean') 
    
    # slope of saturation (kpa °C-1)
    delta = dem.expression(
            '2503*exp((17.27*tair)/(tair+237.3))/(tair+237.3)**2',{
            'tair':tair.subtract(273.15)}).rename('delta')
    
    # Reference evapotranspiration for grass (mm day-1)
    etr_grass_24h =dem.expression(
            '(0.408*delta*(rn)+ cte_psi*(cn/(tair + 273))*ws*(es -ea))/'+
                                          '(delta + cte_psi*(1 + cd*ws))',{
            'delta':delta,
            'rn':rad_24h,
            'cte_psi':const_psci,
            'cn':ee.Number(900),
            'cd':ee.Number(0.34),
            'tair':tair.subtract(273.15),
            'ws':ws,
            'es':es_mean,
            'ea':ea_min}).rename('eto24h')
   
    # Reference evapotranspiration for alfafa (mm day-1)
    etr_alfafa_24h =dem.expression(
            '(0.408*delta*(rn)+ cte_psi*(cn/(tair + 273))*ws*(es -ea))/(delta + cte_psi*(1 + cd*ws))',{
                    'delta':delta,
                    'rn':rad_24h,
                    'cte_psi':const_psci,
                    'cn':ee.Number(1600),
                    'cd':ee.Number(0.38),
                    'tair':tair.subtract(273.15),
                    'ws':ws,
                    'es':es_mean,
                    'ea':ea_min}).rename('etr24h')
  
    return ee.Image.cat([etr_grass_24h,etr_alfafa_24h])

    