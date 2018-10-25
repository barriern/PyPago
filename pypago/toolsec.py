''' Toolbox for the definition of gridded sections '''

import numpy as np
try:
    from param import ee
except ImportError:
    from pypago.sample_param import ee
from pypago.misc import findsecnum, PypagoErrors

def secingrid(fromi, fromj, toi, toj):

    """
    Calculates the sequence of grid points to go from
    (`fromj`, `fromi`) to (`toj`, `toi`), by following the straight line
    that joints the two end points. Adjacent grid points have one side in common.

    :param int fromi: start i index
    :param int fromj: start j index
    :param int toi: end i index
    :param int toj: end i index

    """

    # this offset is only for consistency with the |matlab| veci indices
    off = 1

    fromi = np.atleast_1d(fromi) + off
    fromj = np.atleast_1d(fromj) + off
    toi = np.atleast_1d(toi) + off
    toj = np.atleast_1d(toj) + off

    if toi - fromi == 0:

        vecj = np.linspace(fromj, toj, np.abs(toj - fromj) + 1).astype(np.int)
        veci = fromi * np.ones(len(vecj))

    else:

        if toj - fromj == 0:
            veci = np.linspace(fromi, toi, np.abs(toi - fromi) + 1).astype(np.int)
            vecj = fromj * np.ones(len(veci))

        else:
            veci = np.copy(fromi)
            vecj = np.copy(fromj)

            a = (toj - fromj) / ((toi - fromi).astype(np.float))
            b = toj - toi*a

            while not (veci[-1] == toi) & (vecj[-1] == toj):
                newi = veci[-1] + np.sign(toi - fromi)
                newj = vecj[-1] + np.sign(toj - fromj)
                y = a * newi + b
                x = (newj - b) / (a.astype(np.float))

                if np.abs(x - veci[-1]) >= np.abs(y - vecj[-1]):
                    veci = np.append(veci, newi)
                    vecj = np.append(vecj, vecj[-1])
                else:
                    veci = np.append(veci, veci[-1])
                    vecj = np.append(vecj, newj)

    # correction back to python (removing off=1)
    veci = veci.astype(np.int) - off
    vecj = vecj.astype(np.int) - off

    return vecj, veci

def lengthinsec(veci, vecj, faces, dw, dn):

    """
    Determines the length of section edges
    depending on the cell face.

    :param numpy.array veci: i index of the section grid cells
    :param numpy.array vecj: j index of the section grid cells
    :param numpy.array faces: face of the section grid cells
    :param numpy.array dw: scale factors on the western faces
       of the model cells
    :param numpy.array dn: scale factors on the northern faces
       of the model cells
    :return: length along the section
    """

    lengthvect = np.empty(len(veci))

    for l in xrange(0, len(veci)):
        if faces[l] == 'N':
            lengthvect[l] = np.squeeze(dn[vecj[l], veci[l]])
        elif faces[l] == 'W':
            lengthvect[l] = np.squeeze(dw[vecj[l], veci[l]])

    return lengthvect

def locingrid(lon, lat, matlon, matlat):

    """
    Locates the nearest grid point (lon,lat)
    within the 2D grid defined by the matlon
    and matlat arrays.

    :param numpy.array lon: Array that contains the longitude
       of the sections' points
    :param numpy.array lat: Array that contains the latitude
       of the sections' points
    :param numpy.array matlon: Matrix that contains the longitude
       of the model
    :param numpy.array matlon: Matrix that contains the latitude
       of the model
    :return: the j and i indexes that correspond to the section
       endpoints on the grid
    :rtype: list
    """

    a = np.max(np.abs(matlon[1:, :] - matlon[:-1, :]))
    b = np.max(np.abs(matlon[:, 1:] - matlon[:, :-1]))
    step_x = np.max([a, b])

    c = np.max(np.abs(matlat[1:, :] - matlat[:-1, :]))
    d = np.max(np.abs(matlat[:, 1:] - matlat[:, :-1]))
    step_y = np.max([c, d])

    vecj = np.empty(len(lon))
    veci = np.empty(len(lon))

    for l in xrange(0, len(lon)):

        if lon[l] < np.min(matlon):
            lon[l] = lon[l] + 360

        if lon[l] > np.max(matlon):
            lon[l] = lon[l] - 360

        [j, i] = np.nonzero((matlon <= lon[l] + step_x) &
                            (matlon >= lon[l] - step_x) &
                            (matlat <= lat[l] + step_y) &
                            (matlat >= lat[l] - step_y))

        if len(j) > 1:
            dis = np.empty(len(j))
            for ii in xrange(0, len(j)):
                dis[ii] = np.abs(distance(lat[l], lon[l],
                                          matlat[j[ii], i[ii]],
                                          matlon[j[ii], i[ii]], ee))
            ind = np.argmin(dis)
        elif len(j) == 1:
            ind = 0
        else:
            message = 'Problem in locingrid. Please check that '
            message += '%s N and %s E are within grid boundaries.' %(lat[l], lon[l])
            raise IOError(message)

        vecj[l] = j[ind]
        veci[l] = i[ind]

    veci = veci.astype(np.int)
    vecj = vecj.astype(np.int)

    return vecj, veci

def distance(lat1, lon1, lat2, lon2, geoid):

    """
    Computes the distance between the points
    (lon1, lat1) and (lon2, lat2). This function
    used the :py:func:`mpl_toolkits.basemap.pyproj.Geod.inv`
    function.

    :param float lat1: latitude of the first point
    :param float lon1: longitude of the first point
    :param float lat2: latitude of the second point
    :param float lon2: longitude of the second point
    :param mpl_toolkits.basemap.pyproj.Geod geoid: :py:class:`mpl_toolkits.basemap.pyproj.Geod` object used to compute the distance
    :return: distance between the two points
    :rtype: float
    """

    output = geoid.inv(lon1, lat1, lon2, lat2)

    return output[-1]

def areainsec(veci, vecj, faces, areaw, arean):

    """
    Determines the area of a |pago| section given its indexes and
    faces

    :param numpy.array veci: i-indexes of the section
    :param numpy.array vecj: j-indexes of the section
    :param numpy.array faces: array that contain the faces of the section
       ('W' or 'N')
    :param numpy.array areaW: 3D-array (z, lat, lon) that contains the area of
       the western faces of the grid cells
    :param numpy.array areaN: 3D-array (z, lat, lon) that contains the area of
       the northern faces of the grid cells
    :return: 2D (z, l) area of the |pago| section.
       If faces[p] == 'N', area[:, p] = areaN[:, vecj[l], veci[l]]
    :rtype: numpy.array
    """

    a = arean.shape[0]
    areavect = np.zeros((a, len(veci)))

    for l in xrange(0, len(veci)):
        if faces[l] == 'N':
            areavect[:, l] = arean[:, vecj[l], veci[l]]
        elif faces[l] == 'W':
            areavect[:, l] = areaw[:, vecj[l], veci[l]]

    return areavect


def consec(veci1, vecj1, faces1, orient1, veci2, vecj2, faces2, orient2):

    """
    Function that joins multiple segments
    belonging to one section. At first does a
    preprocessing of the file and then calls the
    :py:func:`pypago.sec._consec._consec`
    function
    """

    # barrier.n --- modified 2015-08-18
    # better "syntax"
    # correction of a huge bug: any replaced by all
    # correction nbarrier: should be < instead of <=
    if np.all(veci2 < veci1[-2]):
        [faces, veci, vecj, orient] = conseclr2(veci2[::-1], vecj2[::-1],
                                                faces2[::-1], orient2[::-1],
                                                veci1[::-1], vecj1[::-1],
                                                faces1[::-1], orient1[::-1])
        finalfaces = faces[::-1]
        finalveci = veci[::-1]
        finalvecj = vecj[::-1]
        finalorient = orient[::-1]

    else:
        [finalfaces, finalveci, finalvecj, finalorient] = \
            conseclr2(veci1, vecj1, faces1, orient1,
                      veci2, vecj2, faces2, orient2)

    return finalfaces, finalveci, finalvecj, finalorient


def conseclr2(veci1, vecj1, faces1, orient1,
              veci2, vecj2, faces2, orient2):

    """
    Function that joins multiple segments
    belonging to one section.
    """

    # JD PAGO WHOI

    if faces1[-1] == 'N':  # vec1 ends with a north face
        if faces2[0] == 'N':  # vec2 starts with a north face
            veci2 = veci2[1:]
            vecj2 = vecj2[1:]
            faces2 = faces2[1:]
            orient2 = orient2[1:]

            if veci1[-2] <= veci1[-1]:
                if veci2[0] <= veci1[-1]:
                    veci1 = veci1[:-1]
                    vecj1 = vecj1[:-1]
                    faces1 = faces1[:-1]
                    orient1 = orient1[:-1]

                while (veci2[0] == veci1[-1]) & (vecj2[0] == vecj1[-1]):
                    veci2 = veci2[1:]
                    vecj2 = vecj2[1:]
                    faces2 = faces2[1:]
                    orient2 = orient2[1:]
                    veci1 = veci1[:-1]
                    vecj1 = vecj1[:-1]
                    faces1 = faces1[:-1]
                    orient1 = orient1[:-1]

            else:

                if veci2[0] >= veci1[-1]:
                    veci1 = veci1[:-1]
                    vecj1 = vecj1[:-1]
                    faces1 = faces1[:-1]
                    orient1 = orient1[:-1]

                while (veci2[1-1] == veci1[-1]) & (vecj2[0] == vecj1[-1]):
                    veci2 = veci2[1:]
                    vecj2 = vecj2[1:]
                    faces2 = faces2[1:]
                    orient2 = orient2[1:]
                    veci1 = veci1[:-1]
                    vecj1 = vecj1[:-1]
                    faces1 = faces1[:-1]
                    orient1 = orient1[:-1]

        else:  # vec2 starts with a west face
            # vec2 starts with a west face then a north face at the same point
            if (faces2[1] == 'N') & (veci2[1] == veci2[0]) & (vecj2[1] == vecj2[0]):

                veci2 = veci2[2:]
                vecj2 = vecj2[2:]
                faces2 = faces2[2:]
                orient2 = orient2[2:]

            # vec2 starts with a west face then west face again
            # (to the north or to the south)
            else:
                veci1 = veci1[:-1]
                vecj1 = vecj1[:-1]
                faces1 = faces1[:-1]
                orient1 = orient1[:-1]

                if vecj2[1] > vecj2[0]:
                    veci2 = veci2[1:]
                    vecj2 = vecj2[1:]
                    faces2 = faces2[1:]
                    orient2 = orient2[1:]

                while (veci1[-1] == veci2[0]) & (vecj1[-1] == vecj2[0]):
                    veci1 = veci1[:-1]
                    vecj1 = vecj1[:-1]
                    faces1 = faces1[:-1]
                    orient1 = orient1[:-1]
                    veci2 = veci2[1:]
                    vecj2 = vecj2[1:]
                    faces2 = faces2[1:]
                    orient2 = orient2[1:]

    else:
        if (veci1[-1] == veci2[0]) & (vecj1[-1] == vecj2[0]):
            veci1 = veci1[:-1]
            vecj1 = vecj1[:-1]
            faces1 = faces1[:-1]
            orient1 = orient1[:-1]

    finalveci = np.append(veci1, veci2)
    finalvecj = np.append(vecj1, vecj2)
    finalfaces = np.append(faces1, faces2)
    finalorient = np.append(orient1, orient2)

    return finalfaces, finalveci, finalvecj, finalorient

def nodouble(veci, vecj):

    """
    Function that removes the double in the veci and vecj arrays
    of a section. The output list are first initialised by the first
    elements of the veci and vecj arrays. Then if the n-1th elements
    is different from the nth element, it is added to the output list

    :param numpy.array veci: i indexes of the section
    :param numpy.array vecj: j indexes of the section
    :return: the veci and vecj indexes but with the double elements
    :rtype: list
    """

    newveci = veci[0]
    newvecj = vecj[0]

    for l in xrange(1, len(veci)):
        if not (veci[l] == veci[l-1]) & (vecj[l] == vecj[l-1]):
            newvecj = np.append(newvecj, vecj[l])
            newveci = np.append(newveci, veci[l])
    return newvecj, newveci



def correct_sections(sec, secname, offset, position):

    """
    Function that allows to correct section staircases.
    Especially useful in order to correct sections' junctions. It extracts
    the section variables (`veci`, `vecj`, etc) along
    the length coordinates, from `offset` to `end` or from
    `start` to `offset`, depending on the value of `position` argument.

    :param list sec: section list that contains the |pypago|
       objects, containing the veci, vecj, vect, ... obtained
       from the :py:func:`pypago.sec.finalisesections`

    :param type secname: section points to correct

    :param int offset: number of points to remove

    :param str position: {'start','end'}
       whether we remove the first (position='start')
       or last (position='end') points.

    """

    nw = findsecnum(sec, secname)

    if position == 'end':

        sec[nw].veci = sec[nw].veci[:-offset]
        sec[nw].vecj = sec[nw].vecj[:-offset]
        sec[nw].vect = sec[nw].vect[:, :, :-offset]
        sec[nw].vecs = sec[nw].vecs[:, :, :-offset]
        sec[nw].vecv = sec[nw].vecv[:, :, :-offset]
        sec[nw].faces = sec[nw].faces[:-offset]
        sec[nw].orient = sec[nw].orient[:-offset]
        sec[nw].areavect = sec[nw].areavect[:, :-offset]
        sec[nw].depthvect = sec[nw].depthvect[:, :-offset]
        sec[nw].lvect = sec[nw].lvect[:-offset]
        sec[nw].lengthvect = sec[nw].lengthvect[:-offset]

    else:
        sec[nw].veci = sec[nw].veci[offset:]
        sec[nw].vecj = sec[nw].vecj[offset:]
        sec[nw].vect = sec[nw].vect[:, :, offset:]
        sec[nw].vecs = sec[nw].vecs[:, :, offset:]
        sec[nw].vecv = sec[nw].vecv[:, :, offset:]
        sec[nw].faces = sec[nw].faces[offset:]
        sec[nw].orient = sec[nw].orient[offset:]
        sec[nw].lengthvect = sec[nw].lengthvect[offset:]
        sec[nw].lvect = sec[nw].lvect[offset:]
        sec[nw].areavect = sec[nw].areavect[:, offset:]
        sec[nw].depthvect = sec[nw].depthvect[:, offset:]

    return sec


def facesinsec(veci, vecj, dire):

    """
    Function which defines the sequence of faces that allow to join
    two section endpoints (i.e for a segment)

    By convention, we tend to favor the north and west side of each grid cell
    as a 'face' of a sequence of gridpoints.

    `veci` and `vecj` are reconstructed to i) drop the points which faces are
    not used, and ii) double the grid points which faces are used twice.

    Faces is a vector of letters:

    - N means that the north face is to be used
    - W means that the west face is to be used

    Orientation is a vector of +1 or -1 that shall be applied to velocities
    in order to make sure that the transport is integrated in the direction
    dir, which is a two letter vector: NE, NW, SE or SW.

    Calculation assumes that vectors are originally pointing toward N or E.

    :param numpy.array veci: i index determined from
       :py:func:`pypago.sec._secingrid`
    :param numpy.array vecj: j index determined from
       :py:func:`pypago.sec._secingrid`
    :param str dir: orientation of the segment
    :return: a tuple that contains the sequence of faces, i and j indexes
       and the orientation for each face
    :rtype: tuple
    :raises ValueError: if either veci or vecj are not monotonous
    :raises IndexError: if length of veci is 1
    """

    if len(veci) == 1:
        error = PypagoErrors('only one point in sequence; no face should be used')
        raise error

    if not monotonous(veci):
        error = PypagoErrors('sequence veci must be monotonous')
        raise error

    if not monotonous(vecj):
        error = PypagoErrors('sequence vecj must be monotonous')
        raise error

    if veci[-1] == veci[0]:
        newveci = np.copy(veci)
        newvecj = np.copy(vecj)
        faces = ['W'] * len(newveci)
    else:

        changedir = 0
        if veci[-1] < veci[0]:
            veci = veci[::-1]
            vecj = vecj[::-1]
            changedir = 1

        indvec = 0
        newveci = veci[indvec]
        newvecj = vecj[indvec]
        if vecj[indvec+1] == vecj[indvec]:
            faces = np.array(['N'])
        else:
            faces = np.array(['W'])
        indvec += 1

        while indvec+1 < len(veci):

            if veci[indvec] == veci[indvec-1]:

                if vecj[indvec] > vecj[indvec-1]:
                    newveci = np.append(newveci, veci[indvec])
                    newvecj = np.append(newvecj, vecj[indvec])
                    faces = np.append(faces, np.array(['W']))

                    if veci[indvec+1] > veci[indvec]:
                        newveci = np.append(newveci, veci[indvec])
                        newvecj = np.append(newvecj, vecj[indvec])
                        faces = np.append(faces, np.array(['N']))

                else:
                    newveci = np.append(newveci, veci[indvec])
                    newvecj = np.append(newvecj, vecj[indvec])

                    if veci[indvec+1] > veci[indvec]:
                        faces = np.append(faces, np.array(['N']))
                    else:
                        faces = np.append(faces, np.array(['W']))

            else:
                if veci[indvec+1] == veci[indvec]:

                    if vecj[indvec+1] < vecj[indvec]:
                        newveci = np.append(newveci, veci[indvec])
                        newvecj = np.append(newvecj, vecj[indvec])
                        faces = np.append(faces, np.array('W'))

                else:
                    newveci = np.append(newveci, veci[indvec])
                    newvecj = np.append(newvecj, vecj[indvec])
                    faces = np.append(faces, np.array('N'))

            indvec = indvec+1

        # when indvec==length(veci)
        # note that new sequence ends with N
        if veci[indvec] == veci[indvec-1]:
            newveci = np.append(newveci, veci[indvec])
            newvecj = np.append(newvecj, vecj[indvec])

            if vecj[indvec] > vecj[indvec-1]:
                newveci = np.append(newveci, veci[indvec])
                newvecj = np.append(newvecj, vecj[indvec])
                faces = np.append(faces, np.array(['W', 'N']))

            else:
                faces = np.append(faces, 'N')

        else:
            newveci = np.append(newveci, veci[indvec])
            newvecj = np.append(newvecj, vecj[indvec])
            faces = np.append(faces, 'N')

        if changedir:
            newveci = newveci[::-1]
            newvecj = newvecj[::-1]
            faces = faces[::-1]

    orientation = make_orientation(newveci, faces, dire)

    return faces, newveci, newvecj, orientation


def make_orientation(newveci, faces, dire):

    """
    Function that allows to extract the orientation of the
    faces, depending on which face and on the orientation of the
    segment

    :param numpy.array newveci: array that contains the
       section indices (used only for recovering the size of output)
    :param numpy.array faces: array that contains the faces of the
       section cells
    :param str dire: a string that contains the orientation of the segment
    """

    orientation = np.zeros(len(newveci))
    for l in xrange(0, len(newveci)):

        if faces[l] == 'W':
            if dire[1] == 'W':
                orientation[l] = -1
            else:
                orientation[l] = 1

        if faces[l] == 'N':
            if dire[0] == 'N':
                orientation[l] = 1
            else:
                orientation[l] = -1

    return orientation


def monotonous(vect):

    """ Determines whether the input array is monotonous

    :param numpy.array vect: Input array
    :return: A boolean (True or False)
    :rtype: bool
    :raises IndexError: if length of vect is 1

    >>> # monotonous
    >>> vect = np.array([1, 2, 3])
    >>> monotonous(vect)
    True

    >>> # constant
    >>> vect = np.array([0, 0, 0])
    >>> monotonous(vect)
    True

    .. todo:: following fail

    >>> # bad input - only 1 item
    >>> vect = np.array([1])
    >>> monotonous(vect)
    Traceback (most recent call last):
       ...
    raise IndexError

    """

    if len(vect) == 1:
        error = PypagoErrors('only one point in sequence; no face should be used')
        raise error

    diff = vect[1:] - vect[0:-1]
    res = (np.all(diff >= 0)) or (np.all(diff <= 0))
    return res
