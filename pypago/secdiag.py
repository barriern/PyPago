
""" Module that handles transport calculation across gridded sections """

# note: if heat transport: multiply by rho0_C
# if fwater: (salt - S0)/S0

import numpy as np
from pypago.disp import PypagoErrors


def net_volume_trans(section, velname='vecv'):

    r"""
    Computation of net volume transport.

    .. math::
        \iint_{S_{o}} U\ dl\ dz

    .. note::

        This function can also be used to compute the net tracer transport
        across a section if the `velname` argument is a transport
        (see also :py:func:`pypago.areadiag.compute_tracer_conv_trans`)

    :param pypago.sections.GridSection section: Gridded section
    :param str name: Name of the velocity attribute

    :return: A numpy array containing the net volume transport
    """

    try:
        vecv = getattr(section, velname)
    except:
        message = 'The %s attribute does not exist' % velname
        raise PypagoErrors(message)

    areavect = np.tile(section.areavect, (vecv.shape[0], 1, 1))

    output = vecv * areavect

    # mask the data
    output = np.ma.masked_invalid(output)

    output = np.ma.sum(output, axis=(1, 2))

    return output


def net_tracer_trans(section, componame, velname='vecv'):

    r"""
    Computation of net tracer transport:

    .. math::
        \iint_{S_{o}} [UT]\ dl\ dz

    :param pypago.sections.GridSection section: Gridded section
    :param str componame: Name of the variable whose transport to compute
    :param str velname: Name of the velocity field
    :return: A numpy array containing the net tracer transport
    """

    try:
        vecv = getattr(section, velname)
        vecv = np.ma.masked_invalid(vecv)
    except:
        message = 'Velocity field %s does not exist' % velname
        raise PypagoErrors(message)

    try:
        vect = getattr(section, componame)
        vect = np.ma.masked_invalid(vect)
    except:
        message = 'Tracer field %s does not exist' % componame
        raise PypagoErrors(message)

    areavect = np.tile(section.areavect, (vecv.shape[0], 1, 1))

    # mask the data
    vect = np.ma.masked_invalid(vect)
    vecv = np.ma.masked_invalid(vecv)
    areavect = np.ma.masked_invalid(areavect)

    output = np.sum(vecv*vect*areavect, axis=(1, 2))

    return output


def total_tracer_trans(section, componame, velname='vecv'):

    r"""
    Computation of total tracer transport, i.e without the net mass transport

    .. math::
        \iint_{S_{o}} [U_{nonet}T]\ dl\ dz

    :param pypago.sections.GridSection section: Gridded section
    :param str componame: Name of the variable whose transport to compute
    :param str velname: Name of the velocity field
    :return: A numpy array containing the net tracer transport
    """

    vecv = remove_spatial_mean(section, velname)

    try:
        vect = getattr(section, componame)
        vect = np.ma.masked_invalid(vect)
    except:
        message = 'Tracer field %s does not exist' % componame
        raise PypagoErrors(message)

    areavect = np.tile(section.areavect, (vecv.shape[0], 1, 1))
    areavect = np.ma.masked_invalid(areavect)
    output = np.ma.sum(vecv*vect*areavect, axis=(1, 2))

    return output


def remove_spatial_mean(section, componame):

    r"""
    Correction of a field by removing the mean:

    .. math::
        T_{nonet} = T - \dfrac{\iint_{S_{o}} T\ dl\ dz}{\iint_{S_{o}}dl\ dz}

    :param pypago.sections.GridSection section: Gridded section
    :param str componame: Name of the tracer field

    """

    try:
        vect = getattr(section, componame)
        vect = np.ma.masked_invalid(vect)
    except:
        message = 'Tracer field %s does not exist' % componame
        raise PypagoErrors(message)

    areavect = np.ma.masked_invalid(section.areavect)
    areavect = np.tile(areavect, (vect.shape[0], 1, 1))   # time, z, l

    nx = areavect.shape[1]
    ny = areavect.shape[2]

    temp = np.sum(areavect * vect, axis=(1, 2)) / np.sum(areavect, axis=(1, 2))  # time
    temp = np.tile(temp, (nx, ny, 1))    # z, l, time
    temp = np.transpose(temp, (2, 0, 1))    # time, z, l

    output = vect - temp

    return output


def overturning_calc(section, velname='vecv'):

    r""" Computation of overturning component:

        .. math::
            U^{ovt}(z) = \dfrac{\int_{l=0}^L U_{nonet}(z, l) dS(z, l)}{\int_{l=0}^L dS(z, l)}

    :param pypago.sections.GridSection section: Gridded section
    :param str velname: Name of the velocity field

    """

    vect = remove_spatial_mean(section, velname)  # time, z, l

    area = section.areavect
    area = np.ma.masked_invalid(area)
    area = np.tile(area, (vect.shape[0], 1, 1))   # time, z, l

    output = np.sum(area * vect, axis=2) / np.sum(area, axis=2)  # time, z

    return output


def overturning_volume_transport(section, velname='vecv'):

    r""" Computation of overturning volume transport:

        .. math::
            U^{ovt}(z) & = \frac{\int_{l=0}^L U_{nonet}(z, l) dS(z, l)}{\int_{l=0}^L dS(z, l)}

            Area(z) & = \int_{l=0}^L dS(z, l)

            V_{ovt} & = max_k\left[\sum_{i=0}^{k}  U_T^{ovt}(z_k) Area(z_k)\right]

    :param pypago.sections.GridSection section: Gridded section
    :param str velname: Name of the velocity field

    """

    # extraction of section area
    area = section.areavect
    area = np.ma.masked_invalid(area)  # nz, nl

    # calculation of zonally averaged field
    vecv_ov = overturning_calc(section, velname)   # time, nz
    area = np.sum(area, axis=1)   # nz
    area = np.tile(area, (vecv_ov.shape[0], 1))   # time, nz

    temp = vecv_ov*area  # time, nz
    temp = np.cumsum(temp, axis=1)   # time, nz

    output = np.max(temp, axis=1)  # ntime

    return output


def overturning_tracer_transport(section, componame, velname='vecv'):

    r""" Computation of overturning volume transport:

        .. math::
            U^{ovt}(z) & = \frac{\int_{l=0}^L U_{nonet}(z, l) dS(z, l)}{\int_{l=0}^L dS(z, l)}

            T^{ovt}(z) & = \frac{\int_{l=0}^L T_{nonet}(z, l) dS(z, l)}{\int_{l=0}^L dS(z, l)}

            Area(z) & = \int_{l=0}^L dS(z, l)

            OVT & = \sum_{z=0}^H Area(z) U^{ovt}(z) T^{ovt}(z)

    :param pypago.sections.GridSection section: Gridded section
    :param str componame: Name of the tracer field
    :param str veloname: Name of the velocity field

    """

    vecv_ov = overturning_calc(section, velname)   # time, z
    vect_ov = overturning_calc(section, componame)   # time, z

    ntime = vecv_ov.shape[0]

    area = section.areavect
    area = np.ma.masked_invalid(area)
    area = np.tile(area, (ntime, 1, 1))   # time, z, l

    temp = vecv_ov * vect_ov * np.sum(area, axis=2)   # time, z

    output = np.sum(temp, axis=1)

    return output


def horizontal_calc(section, velname='vecv'):

    r""" Calculation of horizontal component

        .. math::
            U_{hor} = U_{nonet} - U_{ovt}

        :param pypago.sections.GridSection section: Gridded section
        :param str velname: Name of the velocity field
        :return: A numpy array containing the net tracer transport

    """

    vecv = remove_spatial_mean(section, velname)  # time, z, l
    nl = vecv.shape[2]

    ovmt = overturning_calc(section, velname)   # time, z
    ovmt = np.tile(ovmt, (nl, 1, 1))  # nl, time, z
    ovmt = np.transpose(ovmt, (1, 2, 0))

    output = vecv - ovmt

    return output


def barotropic_calc(section, velname='vecv'):

    r""" Calculation of barotropic component

        .. math::
            U_{bar}(l) = \frac{\int_{z=0}^H U_{hor}(z,l) dS(z, l)}{\int_{z=0}^HdS(z, l)}

        :param pypago.sections.GridSection section: Gridded section
        :param str velname: Name of the velocity field
        :return: A numpy array containing the net tracer transport

    """

    area = section.areavect
    area = np.ma.masked_invalid(area)   # z, l

    vecv = horizontal_calc(section, velname)  # time, z, l
    ntime = vecv.shape[0]

    area = np.tile(area, (ntime, 1, 1))

    output = np.sum(vecv * area, axis=1)/np.sum(area, axis=1)

    return output


def barotropic_volume_transport(section, velname='vecv'):

    r""" Computation of barotropic volume transport

        .. math::
            Area(l) & = \int_{z=0}^{H} dS(z,l)

            V_{bar}(l) & = U_{bar}(l) Area(l)

        :param pypago.sections.GridSection section: Gridded section
        :param str velname: Name of the velocity field
        :return: A numpy array containing the net tracer transport

    """

    vecv_bp = barotropic_calc(section, velname)

    ntime = vecv_bp.shape[0]

    area = section.areavect
    area = np.ma.masked_invalid(area)
    area = np.tile(area, (ntime, 1, 1))

    output = vecv_bp * np.sum(area, axis=1)

    return output


def barotropic_tracer_transport(section, componame, velname='vecv'):

    r""" Computation of baroptropic tracer transport

        .. math::
            U_{bar}(l) & = \frac{\int_{z=0}^H U_{hor}(z,l) dS(z, l)}{\int_{z=0}^HdS(z, l)}

            T_{bar}(l) & = \frac{\int_{z=0}^H T_{hor}(z,l) dS(z, l)}{\int_{z=0}^HdS(z, l)}

            Area(l) & = \int_{z=0}^H dS(z,l)

            BP & = \sum_{l=0}^L U_{bar}(l) T_{bar}(l) Area(l)

        :param pypago.sections.GridSection section: Gridded section
        :param str componame: Name of the variable whose transport to compute
        :param str velname: Name of the velocity field
        :return: A numpy array containing the net tracer transport

    """

    vecv_bp = barotropic_calc(section, velname)
    vect_bp = barotropic_calc(section, componame)

    ntime = vecv_bp.shape[0]

    area = section.areavect
    area = np.ma.masked_invalid(area)  # nz, nl
    area = np.sum(area, axis=0)  # nl

    output = vecv_bp * vect_bp * np.tile(area, (ntime, 1))  # ntime, nl
    output = np.sum(output, axis=1)  # ntime

    return output


def baroclinic_calc(section, componame):

    r""" Computation of baroclinic component

        .. math::
            U_{bc}  = U_{hor} - U_{bp}

        :param pypago.sections.GridSection section: Gridded section
        :param str componame: Name of the variable whose transport to compute
        :return: A numpy array containing the net tracer transport

    """

    vecv_ho = horizontal_calc(section, componame)  # time, z, l
    vecv_bp = barotropic_calc(section, componame)  # time, l

    nz = vecv_ho.shape[1]

    vecv_bp = np.tile(vecv_bp, (nz, 1, 1))  # z, time, l
    vecv_bp = np.transpose(vecv_bp, (1, 0, 2))  # time, z, l

    vecv_bc = vecv_ho - vecv_bp

    return vecv_bc


def baroclinic_tracer_transport(section, componame, velname='vecv'):

    r""" Computation of tracer transport

        .. math::
            U_{bc} & = U_{hor} - U_{bp}

            T_{bc} & = T_{hor} - T_{bp}

            BC & = \iint_{S}  U_{bc} T_{bc} dS

        :param pypago.sections.GridSection section: Gridded section
        :param str componame: Name of the variable whose transport to compute
        :param str velname: Name of the velocity field
        :return: A numpy array containing the net tracer transport

    """

    vecv_bc = baroclinic_calc(section, velname)  # time, nz, nl
    vect_bc = baroclinic_calc(section, componame)   # time, nz, nl

    ntime = vecv_bc.shape[0]

    area = section.areavect   # nz, nl
    area = np.ma.masked_invalid(area)
    area = np.tile(area, (ntime, 1, 1))

    output = np.sum(vecv_bc * vect_bc * area, axis=(1, 2))

    return output
