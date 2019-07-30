#!/usr/bin/env python3

import sys, logging

def handle_exception(exc_type, exc_value, exc_traceback):
    
    ''' Redirect uncaught exceptions (excluding KeyboardInterrupt) 
        through the logging module. 
    '''
    
    fun_name = sys._getframe().f_code.co_name
    log = logging.getLogger(f'{__name__}.{fun_name}')
    
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return
    (log.critical(
        "Uncaught exception",
        exc_info = (exc_type, exc_value, exc_traceback)))
