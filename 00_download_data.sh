#!/bin/bash
mkdir -p data/time_series/CSSE
wget https://github.com/CSSEGISandData/COVID-19/raw/master/csse_covid_19_data/csse_covid_19_time_series/time_series_19-covid-Confirmed.csv -O data/time_series/CSSE/time_series_19-covid-Confirmed.csv
wget https://github.com/CSSEGISandData/COVID-19/blob/master/csse_covid_19_data/csse_covid_19_time_series/time_series_19-covid-Deaths.csv -O data/time_series/CSSE/time_series_19-covid-Deaths.csv
wget https://github.com/CSSEGISandData/COVID-19/blob/master/csse_covid_19_data/csse_covid_19_time_series/time_series_19-covid-Recovered.csv -O data/time_series/CSSE/time_series_19-covid-Recovered.csv

