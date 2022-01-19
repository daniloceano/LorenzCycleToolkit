#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script for computating vertical integrations and partial derivatives

Created by Danilo Couto de Souza
Universidade de São Paulo (USP)
Instituto de Astornomia, Ciências Atmosféricas e Geociências
São Paulo - Brazil

danilo.oceano@gmail.com

"""

import numpy as np
from metpy.units import units


def Differentiate(Data,Axis,AxisName):
    """
    Computates the partial derivative of some data.
    
    It used a second-order finite difference scheme for intermediate levels
    and first-order forward/backward schemes for bottom/top levels
    
    Parameters
    ----------
    Data: xarray.Dataset
        Data containing the variable to be integrated   
    AxisName: string
        the indexer used for the xarray coordante which the differentiation will
        be performed. For example, if AxisName='latitude', the function will 
        differentiate along the latitude coordinate, performing a spatial
        differentiation.
    """
    try:
        DataUnits = Data.metpy.units
        AxisUnits = Axis.metpy.units
        Data = Data.metpy.dequantify()
    except:
        pass
    
    steps = Data[AxisName]
    # StepsUnits = Data[AxisName].metpy.units
    DifArray = Data*np.nan
    # DifArray = Data
    for step in range(len(steps)):
        # Forward difference for first step
        if step == 0:
            diff = Data.sel(**{AxisName:steps[step+1]})-Data.sel(**{AxisName:steps[step]})
            h = steps[step+1].metpy.convert_units(str(AxisUnits))-steps[step].metpy.convert_units(str(AxisUnits))
            DifArray.loc[dict({AxisName:DifArray[AxisName][step]})] = diff/h
        # Backward difference for last step
        elif step == len(steps)-1:
            diff = Data.sel(**{AxisName:steps[step-1]})-Data.sel(**{AxisName:steps[step]})
            h = steps[step-1].metpy.convert_units(str(AxisUnits))-steps[step].metpy.convert_units(str(AxisUnits))
            DifArray.loc[dict({AxisName:DifArray[AxisName][step]})] = diff/h
        # Centred difference for all othe steps
        else:
            diff = Data.sel(**{AxisName:steps[step+1]})-Data.sel(**{AxisName:steps[step-1]})
            h = steps[step+1].metpy.convert_units(str(AxisUnits))-steps[step-1].metpy.convert_units(str(AxisUnits))
            DifArray.loc[dict({AxisName:DifArray[AxisName][step]})] = diff/(2*h)
    return (DataUnits/AxisUnits)*DifArray