#!/bin/bash
mkdir -p data/time_series/CSSE
wget https://github.com/CSSEGISandData/COVID-19/blob/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv -O data/time_series/CSSE/time_series_covid19_confirmed_global.csv
wget https://github.com/CSSEGISandData/COVID-19/blob/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_global.csv -O data/time_series/CSSE/time_series_covid19_deaths_global.csv
wget https://github.com/CSSEGISandData/COVID-19/blob/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_recovered_global.csv -O data/time_series/CSSE/time_series_covid19_recovered_global.csv
wget https://github.com/CSSEGISandData/COVID-19/blob/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_US.csv -O data/time_series/time_series_covid19_confirmed_US.csv
wget https://github.com/CSSEGISandData/COVID-19/blob/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_US.csv -O data/time_series/time_series_covid19_deaths_US.csv
