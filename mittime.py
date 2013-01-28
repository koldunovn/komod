#!/bin/python
# -*- coding: utf-8 -*-


"""Komod time module
Convert timesteps to time and vice versa

steptotime(start_date_time, time_step, present_step)  - convert step to time
timetostep(start_date_time,time_step, present_time)   - convert time to step


Nikolay Koldunov 17 February 2010
"""

# -------------------------------------------------

import time

def steptotime(start_date_time,time_step, present_step):
	"""  Converts timestep of MITgcm model to date   
   	Usage:
     		time_for_model(start_date_time, time_step, present_step )
   	Input:
     		start_date_time  - should be string in form of YYYYDDMMhhmmss
     		time_step        - model timestep (deltaT) in seconds
     		present_step -   - in steps :)
	Output:
		string that represents time in format 
			

        """
	start_date_time_in_python_format = time.strptime(start_date_time,"%Y%m%d%H%M%S")
	start_date_time_in_seconds = time.mktime(start_date_time_in_python_format)
	seconds_from_start_date = time_step*present_step
	seconds_to_time = time.localtime(start_date_time_in_seconds+seconds_from_start_date)
	time_to_string = time.strftime("%Y-%m-%d %H:%M:%S", seconds_to_time)

	return time_to_string

def timetostep(start_date_time,time_step, present_time):
	"""  Converts date to MITgcm model timestep   
   	Usage:
     		timetostep(start_date_time, time_step, present_time )
   	Input:
     		start_date_time - should be string in form of YYYYDDMMhhmmss
     		time_step - model timestep (deltaT) in seconds
     		present_time - should be string in form of YYYYDDMMhhmmss
	Output: time step 
		
        """
	start_date_time_in_python_format = time.strptime(start_date_time,"%Y%m%d%H%M%S")
	start_date_time_in_seconds = time.mktime(start_date_time_in_python_format)
	present_time_in_python_format   = time.strptime(present_time,"%Y%m%d%H%M%S")
	present_time_in_seconds = time.mktime(present_time_in_python_format)

	seconds_from_start_date = present_time_in_seconds  - start_date_time_in_seconds
	present_time_step = seconds_from_start_date/time_step

	return present_time_step

