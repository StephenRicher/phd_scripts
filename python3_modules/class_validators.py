#!/usr/bin/env python3

import re

class IntRange():
    ''' Descriptor to validate attribute format with custom
    integer range. '''
    
    def __init__(self, minimum, maximum):
        self.minimum = minimum
        self.maximum = maximum
    
    def __get__(self, instance, owner):
        return instance.__dict__[self.name]

    def __set__(self, instance, value):
        if not isinstance(value, int):
            raise TypeError(f'Error: {self.name} must be integer.')
        elif not self.minimum <= value <= self.maximum:
            raise ValueError(f'Error: {self.name} out of range.')
        instance.__dict__[self.name] = value

    def __set_name__(self, owner, name):
        self.name = name

class RegexMatch():

    def __init__(self, regex):
        self.regex = regex
    
    def __get__(self, instance, owner):
        return instance.__dict__[self.name]

    def __set__(self, instance, value):
        if not isinstance(value, str):
            raise TypeError(f'Error: {self.name} must be string.')
        elif not re.match(self.regex, value):
            raise ValueError(f'Invalid format in {self.name}.')
        instance.__dict__[self.name] = value

    def __set_name__(self, owner, name):
        self.name = name
