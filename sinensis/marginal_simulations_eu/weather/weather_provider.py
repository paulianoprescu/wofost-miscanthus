import os
import re
import logging
import requests
import yaml
import json
import pandas as pd

from pcse.base import WeatherDataContainer, WeatherDataProvider
from pcse.util import reference_ET, wind10to2
from pcse.exceptions import PCSEError
from pcse.input import WOFOST72SiteDataProvider
from pcse.base import ParameterProvider
from pcse.engine import Engine
from pcse.models import Wofost72_WLP_CWB, Wofost72_PP
import config

no_conv = lambda x: x
mm_to_cm = lambda x: x/10.
kJ_to_J = lambda x: x*1000.
to_date = lambda x: x.date()

class WOFOSTWebWeatherDataProvider(WeatherDataProvider):
    variable_renaming = [("temperature_max", "TMAX", no_conv),
                         ("temperature_min", "TMIN", no_conv),
                         ("temperature_avg",  "TEMP", no_conv),
                         ("vapourpressure", "VAP", no_conv),
                         ("windspeed", "WIND", wind10to2),
                         ("precipitation", "RAIN",mm_to_cm),
                         ("radiation", "IRRAD", no_conv),
                         ("snowdepth", "SNOWDEPTH", no_conv),
                         ("day", "DAY", to_date)]
    angstA = 0.25
    angstB = 0.45
    ETmodel = "PM"

    def __init__(self, inputs):
        WeatherDataProvider.__init__(self)
        url = f'{config.weather_host}/api/v1/get_agera5'
        r = requests.get(url, params=inputs)
        data = json.loads(r.text)
        self.elevation = data["location_info"]["grid_agera5_elevation"]
        self.longitude = data["location_info"]["grid_agera5_longitude"]
        self.latitude = data["location_info"]["grid_agera5_latitude"]
        df = pd.DataFrame(data["weather_variables"])
        df["day"] = pd.to_datetime(df.day)
        df_new = self._rename_columns(df)
        for daily_weather in df_new.itertuples(index=False):
            self._make_WeatherDataContainer(daily_weather)

    def _rename_columns(self, df):
        """Renames columns and applies conversions for use in PCSE models
        """
        df_new = pd.DataFrame()
        for old_name, new_name, conversion in self.variable_renaming:
            df_new[new_name] = df[old_name].apply(conversion)

        return df_new

    def _make_WeatherDataContainer(self, daily_weather):
        """Converts a record of daily weather data into a WeatherDataContainer and stores it.
        """
        t = daily_weather._asdict()
        t.update({"LAT": self.latitude, "LON": self.longitude, "ELEV": self.elevation})

        # Reference evapotranspiration in mm/day
        try:
            E0, ES0, ET0 = reference_ET(ANGSTA=self.angstA, ANGSTB=self.angstB, ETMODEL=self.ETmodel, **t)
        except ValueError as e:
            msg = (f"Failed to calculate reference ET values on {daily_weather.DAY}. " +
                   f"With input values:\n {t}\n Due to error: {e}")
            raise PCSEError(msg)

        # update record with ET values value convert to cm/day
        t.update({"E0": E0 / 10., "ES0": ES0 / 10., "ET0": ET0 / 10.})

        # Build and store weather data container from dict 't' for this day
        wdc = WeatherDataContainer(**t)
        self._store_WeatherDataContainer(wdc, daily_weather.DAY)