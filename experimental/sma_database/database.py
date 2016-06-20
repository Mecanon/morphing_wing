# -*- coding: utf-8 -*-
"""
Code using pytables to store SMA data

Created on Thu Feb 11 14:05:32 2016

Must haves:
- Need to adjustable for different sources
- Able to store diverve properties
- Store experimental analysis
- Learn how to not modify if exists
- Easy to create an input
- Embebbed fitting
- Multiple SMA models?
- Must have units
@author: Pedro
"""

import numpy as np
import tables as tb

experimental = False
#
f = tb.openFile('sma_data.h5', 'w')

wire_group = f.createGroup('/', 'wires', 'Data and porperties for wires')

#==============================================================================
# Data from source: variable input according to source
#==============================================================================
# In case of data from experiment
# Create source/experiment table (date, equipment, operator, vendor)
def source(object = None, origin = "paper"):
    if object == None:
      raise Exception("No object defined")
      
    elif origin == "paper":
        object.year = 
        object.author = 
        object.institution = 
    elif origin == "experiment":
        object.laboratory = 
        object.year = 
        object.month = 
        object.day = 
        object.equipment = 
        object.operator = 
        object.fabricant = 
        
class source_exp(tb.IsDescription):
    laboratory = tb.StringCol(20)   #20-character String
    year = tb.Int64Col()            #Signed 64-bit integer
    month = tb.Int64Col()            #Signed 64-bit integer
    day = tb.Int64Col()            #Signed 64-bit integer    
    equipment = tb.StringCol(20)   #20-character String
    operator = tb.StringCol(20)   #20-character String
    fabricant = tb.StringCol(20)   #20-character String

# In case of data from paper or etc
class source_paper(tb.IsDescription):
    laboratory = tb.StringCol(20)   #20-character String
    year = tb.Int64Col()            #Signed 64-bit integer
    month = tb.Int64Col()            #Signed 64-bit integer
    day = tb.Int64Col()            #Signed 64-bit integer    
    equipment = tb.StringCol(20)   #20-character String
    operator = tb.StringCol(20)   #20-character String
    fabricant = tb.StringCol(20)   #20-character String
    
f.createTable(wire_group, 'source', source_exp, 'Data about the experiment')

##==============================================================================
## Properties
##==============================================================================
class properties(tb.IsDescription):
    idnumber = tb.Int64Col()
    E1 = tb.Float64Col()
    density = tb.Float64Col()

f. createTable(wire_group, 'properties', properties, 'All material properties')   
BREAK
#==============================================================================
# Experiments
#==============================================================================
experiments_group = f.createGroup(wire_group, 'experiments', 'Contains all' + 
                                  ' experimental data')

electrical_group = f.createGroup('/wires/experiments/', 'electrical_characterization', 
                                 'All data and calculated plots for ' + 
                                 'electrical characterization')
f.createArray(electrical_group, 'data', properties, 'All material properties')
tensile_group = f.createGroup('/wires/experiments/', 'tensile_test',
                              'All data from tensile tests')
dsc_group = f.createGroup('/wires/experiments/', 'DSC', 
                          'Data from the differential scanning calorimetry.')


#source_dtype = np.dtype([('laboratory', 'S20'),
#               ('year', int),
#               ('month', int),
#               ('day', int),
#               ('equipment', 'S20'),
#               ('operator', 'S20'),
#               ('fabricant', 'S20')])
#
##In case of data from a paper or etc
#source_dtype = np.dtype([('laboratory', 'S20'),
#               ('year', int),
#               ('month', int),
#               ('day', int),
#               ('equipment', 'S20'),
#               ('operator', 'S20'),
#               ('fabricant', 'S20')] )
#
#source = np.array([('Mecanon', 2016, 02, 12, 'Instron', 'Eduardo', 
#                    'Sandinox')], dtype=source_dtype)  
#

#f.createTable('/wires', 'source', source)
## Create material properties table
#prop_dtype = np.dtype([('id', int),
#                      ('E1', float),
#                      ('density', float)])
#
#properties = np.array([(1,10.,100.)], dtype=prop_dtype)
#f.createTable('/wires', 'properties', properties)
#f.flush()