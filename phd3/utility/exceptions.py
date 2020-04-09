#!/usr/bin/env python3

"""Custom exceptions for the turbopy python package"""

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