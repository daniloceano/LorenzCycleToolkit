#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 20 10:05:52 2022

@author: danilocoutodsouza
"""
import aux
from EnergyContents import EnergyContents
from ConversionTerms import ConversionTerms
from BoxData import BoxData
from metpy.units import units

import dataclasses
import sys
from typing import List, Any

USAGE = f"Usage: python {sys.argv[0]} [--help] | file fvar min_lon, max_lon, min_lat, max_lat]"

@dataclasses.dataclass
class Arguments:
    file: str
    fvar: str
    min_lon: float
    max_lon: float
    min_lat: float
    max_lat: float

def check_type(obj):
    for field in dataclasses.fields(obj):
        value = getattr(obj, field.name)
        print(
            f"Value: {value}, "
            f"Expected type {field.type} for {field.name}, "
            f"got {type(value)}"
        )
        if type(value) != field.type:
            print("Type Error")
        else:
            print("Type Ok")

def validate(args: List[str]):
    if len(args) > 2:
        for i in range(2,len(args)):
                args[i] = float(args[i])
    try:
        arguments = Arguments(*args)
    except TypeError:
        raise SystemExit(USAGE)
    check_type(arguments)

def tests_args() -> None:
    args = sys.argv[1:]
    if not args:
        raise SystemExit(USAGE)

    if args[0] == "--help":
        print(USAGE)
    else:
        validate(args)
        
        
def main():
    
    tests_args()
    
    data = aux.get_data(*sys.argv[1:])    
    print(data)

if __name__ == "__main__":
    main()