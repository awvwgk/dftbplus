#! /usr/bin/env python3
#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2006 - 2025  DFTB+ developers group                           #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#
#
'''Potential/charge data conversion to cube format.'''

import argparse
import dptools.gridsio as gridsio
import dptools.grids as grids
import dptools.common as common
import numpy as np


USAGE = """Convert the potential or charge data file from dftb+ transport
calculation to a Gaussian Cube format.
"""


def main(cmdlineargs=None):
    '''Main driver for makecube.

    Args:
        cmdlineargs: List of command line arguments. When None, arguments in
            sys.argv are parsed (Default: None).
    '''
    args = parse_arguments(cmdlineargs)
    makecube(args)


def parse_arguments(cmdlineargs=None):
    '''Parses command line arguments.

    Args:
        cmdlineargs: List of command line arguments. When None, arguments in
            sys.argv are parsed (Default: None).
    '''
    parser = argparse.ArgumentParser(description=USAGE)

    helpstring = 'file containing X vector coordinates (default Xvector.dat)'
    parser.add_argument('--xvec', default='Xvector.dat',
                        help=helpstring, dest='xvec')
    helpstring = 'file containing Y vector coordinates (default Yvector.dat)'
    parser.add_argument('--yvec', default='Yvector.dat',
                        help=helpstring, dest='yvec')
    helpstring = 'file containing Z vector coordinates(default Zvector.dat)'
    parser.add_argument('--zvec', default='Zvector.dat',
                        help=helpstring, dest='zvec')
    helpstring = "Reference, If a present, the difference between the input" \
                 "file and the reference is printed"
    parser.add_argument('--reference', help=helpstring, dest='reference')
    helpstring = 'Output format (cube or vtk). If not specified, guessed' \
                 'from file extension'
    parser.add_argument('--format', help=helpstring, dest='format')
    helpstring = 'file containing grid data'
    parser.add_argument('source', help=helpstring)
    helpstring = 'output cube/vtk file name'
    parser.add_argument('dest', help=helpstring)

    args = parser.parse_args(cmdlineargs)

    return args


def makecube(args):
    '''Converts the potential or charge data file from DFTB+ transport
    calculations to a Gaussian Cube format.

    Args:
        args: Containing the obtained parsed arguments.
    '''
    # Build the grid
    with open(args.xvec, 'r') as xvecfile:
        xvec = np.array(xvecfile.read().split(), dtype=float)
    with open(args.yvec, 'r') as yvecfile:
        yvec = np.array(yvecfile.read().split(), dtype=float)
    with open(args.zvec, 'r') as zvecfile:
        zvec = np.array(zvecfile.read().split(), dtype=float)

    origin = (xvec[0], yvec[0], zvec[0])
    xres = xvec[1] - xvec[0]
    yres = yvec[1] - yvec[0]
    zres = zvec[1] - zvec[0]
    basis = ((xres, 0.0, 0.0), (0.0, yres, 0.0), (0.0, 0.0, zres))
    ranges = ((0, len(xvec)), (0, len(yvec)), (0, len(zvec)))
    grid = grids.Grid(origin, basis, ranges)

    # Build the data object
    sourcefile = common.openfile(args.source)
    vec = np.array(sourcefile.read().split(), dtype=float)
    if args.reference is not None:
        referencefile = common.openfile(args.reference)
        ref = np.array(referencefile.read().split(), dtype=float)
        vec -= ref
    # Note: potential.dat is stored in the order x0y0z0...x0y0zn, x0y1z0...
    # which correspond to a row major order (right index changing faster)
    # in a A(x,y,z) notation
    griddata = grids.GridData(grid, vec.reshape(grid.shape, order='C'))

    if args.format is None:
        if args.dest.endswith('.cub') or args.dest.endswith('.cube'):
            outformat = 'cube'
        elif args.dest.endswith('.vtk'):
            outformat = 'vtk'

    if outformat == 'vtk':
        gridsio.scalarvtk(args.dest, griddata,
                          varname=args.dest.split('.')[0])
    elif outformat == 'cube':
        gridsio.cube(args.dest, griddata)
    else:
        raise ValueError('Cannot determine output format in '
                         'makecube:unknown file extension')


if __name__ == "__main__":
    main()
