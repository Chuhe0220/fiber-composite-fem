# -*- coding: utf-8 -*-
"""
Created on Tue Nov 30 09:19:04 2021

@author: Hannes Grimm-Strele
"""

from DatabaseHDF5 import *

DB = Database('Database.h5',write=False)
ModelAttrLinearElastic = {DATABASE_MODELTYPE:'LinearElasticThermalExpansion'}
numberOfPoints = DB.getHighestNumber(DB.DB,DATABASE_POINT)


for index in range(1, numberOfPoints+1):
    fiberOrientation = DB.getPoint(index).attrs[DATABASE_FIBERORIENTATION]
    fiberOrientation_zz = 1.0 - fiberOrientation[0] - fiberOrientation[1]
    try:
        fiberVolumeFraction = DB.getPoint(index).attrs[DATABASE_FIBERVOLUMEFRACTION]
        print("Data point {0}: fiber orientation = ({1:6.4f}, {2:6.4f}, {3:6.4f}), fiber volume fraction = {4:4.1f} %".format(index, fiberOrientation[0], fiberOrientation[1], fiberOrientation_zz, fiberVolumeFraction))
    except:
        print("Data point {0}: fiber orientation = ({1:6.4f}, {2:6.4f}, {3:6.4f})".format(index, fiberOrientation[0], fiberOrientation[1], fiberOrientation_zz))

    data = DB.getData(DATABASE_MODEL_MATRIX, index, ModelAttrLinearElastic)
    stiffnessTensor = data[b'C']
    print("Stiffness tensor = {0}".format(stiffnessTensor))
    print("")

