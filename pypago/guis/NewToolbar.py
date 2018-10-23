# -*- coding: utf-8 -*-

""" Module that handles the Toolbar used in the Guis"""

from matplotlib.backends.backend_tkagg import NavigationToolbar2TkAgg


class NewToolbar(NavigationToolbar2TkAgg):

    """
    Toolbar of the Matplotlib window (inherited from the
    :py:class:`NavigationToolbar2TkAgg`)

    Modified in order to remove the display of the mouse position.
    barrier.n
    """

    def set_message(self, msg):

        """ Function that overwrite the default set_message function """
        pass
