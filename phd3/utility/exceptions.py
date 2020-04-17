#!/usr/bin/env python3
"""
Author  ==>> Matthew R. Hennefarth
Date    ==>> April 16, 2020
"""

__all__ = [
    'Alarm',
    'ParameterError',
    'DefineError'
]

class Alarm(Exception):
    pass

class ParameterError(Exception):
    pass

class DefineError(Exception):
    pass
