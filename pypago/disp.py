# -*- coding: utf-8 -*-

"""
Module that handles output of errors and log messages
"""

from __future__ import print_function
import sys
import linecache


class PypagoErrors(Exception):

    """
    This class is made to handle errors/warnings in |pypago|

    :param str message: Message to print in the error (value of the exception)

    :param str type: A string to indicate what kind of errors is raised
        (should be 'Warning' or 'Error').

    """

    def __init__(self, message):

        """
        Initialisation of the Exception class

        :param str message: Message to print in the error (value of the
           exception)
        :param str type: A string to indicate what kind of errors is raised
           (should be 'Warning' or 'Error').
        """

        # initialisation as suggested by PyLint
        super(PypagoErrors, self).__init__()

        self.message = '\n@@@@@@@@@@@@@@@@@@ PypagoErrors @@@@@@@@@@@@@@@@@@\n'
        self.message += message
        self.message += '\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n'

        self.value = self._make_str()

    def __str__(self):

        """ Function that redefines the __str__ function """

        return self.value

    def _make_str(self):

        """
        Function which handles the creation of an error/warning message
        containing the name of the exception, of the :file:`.py` file,
        the line number and the line content that raised the exception/warning.

        Inspired from `Apogentus answer
        <http://stackoverflow.com/questions/14519177/python-exception-handling-line-number>`
        """

        # sys.exc_info() returns info on the exception that is handled
        # tb is the traceback, in which we are interested
        # if no exception are provided, then arg returns a tuple of None
        arg = sys.exc_info()

        if arg[2] is None:

            filename = __name__ + '.py'
            strout = '{}'.format(self.message)
        else:

            # f is the frame associated with the traceback
            tb_frame = arg[2].tb_frame

            # this is the line number where the exception has been raised
            lineno = arg[2].tb_lineno

            # this is the filename
            filename = tb_frame.f_code.co_filename

            # allows to recover the line responsible for the warning
            linecache.checkcache(filename)
            line = linecache.getline(filename, lineno, tb_frame.f_globals)
            strout = '{}  (in {}, line {}: "{}")'.format(self.message,
                                                         filename, lineno,
                                                         line.strip())
        return strout
