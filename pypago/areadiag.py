
""" Module that handles domain calculation into domains """


import numpy as np
from sample_param import *
try:
    from param import *
except:
    pass
from pypago.disp import PypagoErrors
import pypago.misc
import pypago.secdiag as secdiag


def compute_tracer_conv(area, sections, componame, velname='vecv'):

    r""" Calculation of tracer convergence within the domain:

        .. math::

            \sum_{i=1}^N\left[\iint_{S_{i}} [U T]\ dl\  dz\right]

        with *i* the section index.

    :param pypago.area.Area area: Closed domain
    
    :param list sections: List of pypago.sections.GridSection objects,
     containing the sections that define the domain

    :param str componame: Name of the tracer variable
    
    :param str velname: Name of the velocity field

    :return: A numpy array containing the tracer convergence.

    """

    cpt = 1

    # loop over the section domains 
    for name, sign in zip(area.secnames, area.signs):

        # recovering the sectiion
        idsec = pypago.misc.findsecnum(sections, name)

        trans = secdiag.net_tracer_trans(sections[idsec], componame, velname)

        if cpt:
            # initialisation of the ocean heat convergence
            heatconv = np.zeros(len(trans), dtype=np.float)
            cpt = 0

        heatconv += sign * trans

    return heatconv


def compute_tracer_conv_trans(area, sections, velname):

    r""" Calculation of tracer convergence within the domain:

        .. math::

            \sum_{i=1}^N\left[\iint_{S_{i}} [U T]\ dl\  dz\right]

        with *i* the section index.
       
        In this function, the :samp:`velname` argument is already a transport.

        .. note::
         
            Use this function if the heat transport has been stored by the model,
            since it provides a better precision than offline calculation

    :param pypago.area.Area area: Closed domain
    
    :param list sections: List of pypago.sections.GridSection objects,
     containing the sections that define the domain

    :param str velname: Name of the transport field.

    :return: A numpy array containing the tracer convergence.

    """

    cpt = 1

    # loop over the section domains 
    for name, sign in zip(area.secnames, area.signs):

        # recovering the sectiion
        idsec = pypago.misc.findsecnum(sections, name)

        trans = secdiag.net_volume_trans(sections[idsec], velname)

        if cpt:
            # initialisation of the ocean heat convergence
            heatconv = np.zeros(len(trans), dtype=np.float)
            cpt = 0

        heatconv += sign * trans

    return heatconv



def volume_content(area, componame):

    """ Computation of volume integrated tracer content:

        .. math::
            HC = \iiint_V T dV

    """

    try:
        vect = getattr(area, componame)
        vect = np.ma.masked_invalid(vect)
    except:
        message = 'The %s attribute does not exist' %velname
        raise PypagoErrors(message)

    if vect.ndim != 3:
        message = 'The variable must be 3D. This program '
        message += 'will be stopped'
        raise PypagoErrors(message)

    ntime, nz, nspace = vect.shape

    volume = area.volume
    volume = np.ma.masked_invalid(volume)

    hc = np.sum(volume * vect, axis=(1, 2))

    return hc
    

def surface_content(area, componame):
    
    """ Computation of surface integrated tracer content. 
        
        - If the tracer field is 2D (time, space)

        .. math::
            HC = \iint_S T dS
        
        - If the tracer field is 3D (time, space)

        .. math::
            HC = \iint_S T(k=0) dS

    """

    try:
        vect = getattr(area, componame)
        vect = np.ma.masked_invalid(vect)
    except:
        message = 'The %s attribute does not exist' %velname
        raise PypagoErrors(message)
    
    if vect.ndim not in [2, 3]:
        message = 'The variable must be 2D or 3D. This program '
        message += 'will be stopped'
        raise PypagoErrors(message)
    
    surf = area.surface
    surf = np.ma.masked_invalid(surf)

    if vect.ndim == 3:
        vect = vect[:, 0, :]
        

    ntime, nspace =  vect.shape
    surf = np.tile(surf, (ntime, 1))
    output = np.sum(surf * vect, axis=1)
    

    return output


if __name__ == '__main__':

    import pypago.pyio

    data = pypago.pyio.load('/home/nbarrier/Python/pago/trunk/doc_pypago/data/natl_datadom.pygo')
    area = data[0]
    print area

    hc = volume_content(area, 'temp')

    hf = surface_content(area, 'hf')
    
    sst = surface_content(area, 'temp')

    
