#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import os
import numpy as np
from scipy.spatial import Delaunay
import matplotlib.tri as tri
import copy
#import mpi4py
#from mpi4py import MPI
import h5py
import xml.etree.ElementTree as ET

encoding = 'utf-8'

DATABASE_FIBERORIENTATION =    'FiberOrientation'
DATABASE_FIBERVOLUMEFRACTION = 'FiberVolumeFraction'
DATABASE_FIBERLENGTH =         'FiberLength'
DATABASE_TEMPERATURE =         'Temperature'
DATABASE_MODELTYPE =           'Type'
DATABASE_MODELGEOMETRY =       'Geometry'

DATABASE_MESH =                 'Mesh'
DATABASE_CONNECTIVITY =         'Connectivity'
DATABASE_MODEL =                'Model'
DATABASE_POINT =                'Point'

DATABASE_NAME =                 'Name'
DATABASE_NUMBER =               'Number'


DATABASE_MODEL_MATRIX =         'SystemMatrices'
DATABASE_MODEL_PODMODES =       'Modes'
DATABASE_MODEL_PODMODE =        'Mode'
DATABASE_MODEL_VOLUMEFIELD =    'VolumeField'
DATABASE_MODEL_VOLUMEFIELDS =   'VolumeFields'


def print_attrs(name, obj):
    print(name)
    for key, val in obj.attrs.items():
        print("    {0}: {1}".format(key, val))


def indent(elem, level=0):
    i = "\n" + level*"  "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for elem in elem:
            indent(elem, level+1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i

def get_Address(name, name1):
    if name ==name1:
        return name

class PointIndicator:
    def __init__(self,AttributeName,AttributeDimension,Tolerance):
        self.AttributeDimension=AttributeDimension
        self.AttributeName=np.string_(AttributeName)
        self.Tolerance=Tolerance
        self.TypeOfValue =type(np.zeros((AttributeDimension)))
        self.shape = np.zeros((AttributeDimension)).shape

class Mesh:
    def __init__(self,PI,meshgeneration):
        self.indicator=PI
        self.generator=meshgeneration
        self.Connectivity=None


class Point(object):
    def __init__(self,PI,Attributes):
        self.indicator=PI
        self.Attributes=Attributes


class Model(object):
    def __init__(self,Attributes,Data,Number=0):
        self.Attributes=Attributes
        self.Data=Data
        self.Number=Number

class Database(object):
    def __init__(self,dbfile,new=False,PointIndicator=None,write=False,compression=True):
        self.file=dbfile # filename str
        self.PIList=PointIndicator #dict PointIndicator
        self.write=write #bool
        self.new=new
        self.MeshIsAvailable=True
        self.compression=compression
        self.open()

    def close(self):
        self.DB.close()

    def compress(self):
        self.close()
        print("compress data")
        os.system("h5repack -f GZIP=4 "+self.file+" "+self.file.replace(".h5","_r.h5"))
        os.remove(self.file)
        os.rename(self.file.replace(".h5","_r.h5"),self.file)
        self.open()

    def open(self):
        if self.new:
            self.write=True

        if self.new and self.PIList==None:
            raise ValueError("can not initialize new database without point indicator")

        if not self.write:
            try:
                self.DB = h5py.File(self.file, 'r')
                self.readPointIndicator()
            except:
                raise ValueError("can not open database: "+self.file)
        else:
            if self.new:
                if(os.path.isfile(self.file)):
                    try:
                        os.remove(self.file)
                    except:
                        raise ValueError("can not delete database: "+self.file)
                print("newFile")
                self.DB = h5py.File(self.file, 'w')
                print("newFile2")
                self.writePointIndicator()
                print("newFile3")
            else:
                try:
                  self.DB = h5py.File(self.file, 'r+')
                  self.readPointIndicator()
                except:
                  raise ValueError("can not initialize new database: "+self.file)


    def readPointIndicator(self):

        self.PIList=[]
        IDList = self.DB.attrs["IDList"]
        DimensionList = self.DB.attrs["IDDimensionList"]
        ToleranceList = self.DB.attrs["IDTolerance"]

        for i in range(0,len(IDList)):
            self.PIList.append(PointIndicator(IDList[i],DimensionList[i],ToleranceList[i]))


    def writePointIndicator(self):
        IDList=[]
        DimensionList=[]
        ToleranceList=[]

        for p_i in self.PIList:
            IDList.append(p_i.AttributeName)
            DimensionList.append(p_i.AttributeDimension)
            ToleranceList.append(p_i.Tolerance)

        self.DB.attrs["IDList"] = IDList
        self.DB.attrs["IDDimensionList"] = DimensionList
        self.DB.attrs["IDTolerance"] = ToleranceList


    def getPointIndicator(self):
        return self.PIList

    def getPoint(self,someVal):
        if type(someVal) is Point:
            point=self.getPointByPoint(someVal)
        elif type(someVal) is h5py._hl.group.Group:
            point=someVal
        elif type(someVal) is dict:
            point=self.getPointByIndicatorValue(someVal)
        elif type(someVal) is int or type(someVal) is np.int32 or type(someVal) is np.int64:
            point=self.getPointByNumber(someVal)
        else:
            print("unknown value type: {0}".format(type(someVal)))
            point=-1
        return point

    def getPointByPoint(self,point):
        for key in self.DB.keys():
            if DATABASE_POINT in key:
                if np.linalg.norm(self.DB[key].attrs[self.PI[DATABASE_NAME]]-point.Attributes[self.PI[DATABASE_NAME]])<self.PI["Tolerance"]:
                    return self.DB[key]
        return -1

    def getPointByIndicatorValue(self,value):
        for key in self.DB.keys():
            if DATABASE_POINT in key:
                ispoint=True
                for i in range(0,len(self.PIList)):
                    try:
                        AttributeName=self.PIList[i].AttributeName
                        val=value[self.PIList[i].AttributeName.decode(encoding)]
                        if type(val)==list:
                            val=np.asarray(val)
                        if type(val)==float:
                            val=np.asarray([val])
                        if type(val)==int:
                            val=np.asarray([val])
                        dbval=self.DB[key].attrs[AttributeName]
                        tol=self.PIList[i].Tolerance
                    except:
                        ispoint=False
                        break

                    if self.PIList[i].shape!=val.shape:
                        ispoint=False
                        break

                    if np.linalg.norm(dbval-val)>tol:
                        ispoint=False
                        break
                if ispoint:
                    return self.DB[key]

        print("can't find point for : {0}".format(value))
        return -1

    def getPointByNumber(self,Number):
        for key in self.DB.keys():
            if DATABASE_POINT in key:
                if self.DB[key].attrs[DATABASE_NUMBER]==Number:
                    return self.DB[key]
        print("can't find point with Number {0}".format(Number))
        return -1

    def getModel(self,DatabasePoint,someVal):
        if type(someVal) is Model:
            model=self.getModelByModel(DatabasePoint,someVal)
        elif type(someVal) is int:
            model=self.getModelByNumber(DatabasePoint,someVal)
        elif type(someVal) is dict:
            model=self.getModelByDict(DatabasePoint,someVal)
        elif type(someVal) is h5py._hl.group.Group:
            return someVal
        else:
            print("unknown value type: {0}".format(type(someVal)))
            model=-1
        return model


    def getModelByModel(self,DatabasePoint,Model):
        for key in DatabasePoint.keys():
            if DATABASE_MODEL in key:
                DatabaseModel=DatabasePoint[key]
                if self.compareModel(DatabaseModel,Model):
                    return DatabaseModel
        return -1

    def getModelByDict(self,DatabasePoint,ModelDict):
        for key in DatabasePoint.keys():
            if DATABASE_MODEL in key:
                DatabaseModel=DatabasePoint[key]
                ismodel=True
                for key,val in DatabaseModel.attrs.items():
                    if key!=DATABASE_NUMBER:
                        try:
                            if DatabaseModel.attrs[key].decode(encoding)!=ModelDict[key]:
                                ismodel=False
                                break
                        except:
                            ismodel=False
                            break
                if ismodel:
                    return DatabaseModel
        print("can't find Model in point: {0} with Attrs {1}".format(DatabasePoint,ModelDict))

        return -1


    def getModelByNumber(self,DatabasePoint,Number):
        for key in DatabasePoint.keys():
            if DATABASE_MODEL in key:
                DatabaseModel=DatabasePoint[key]
                if DatabaseModel.attrs[DATABASE_NUMBER]==Number:
                    return DatabaseModel
        print("can't find Model in point: {0} with Number {1}".format(DatabasePoint,Number))
        return -1

    def getHighestNumber(self,obj,name):
        Number=0
        for key in obj.keys():
            if name in key:
                if obj[key].attrs[DATABASE_NUMBER]>Number:
                    Number=obj[key].attrs[DATABASE_NUMBER]
        return Number


    def addPoint(self,point):
        if not self.write:
            raise ValueError("(addPoint): No write intent")
        DatabasePoint = self.getPoint(point)
        if DatabasePoint == -1:
            print('Add Point: {0} to database'.format(point))
            PointNumber = self.getHighestNumber(self.DB,DATABASE_POINT)
            DatabasePoint=self.DB.create_group('Point_'+str(PointNumber+1))
            for i in range(0,len(self.PIList)):
                try:
                    val=point[self.PIList[i].AttributeName.decode(encoding)]
                except:
                    raise ValueError("point indicator mismatch")
                DatabasePoint.attrs[self.PIList[i].AttributeName]=val
            DatabasePoint.attrs[DATABASE_NUMBER]=PointNumber+1
            self.meshing()

        return DatabasePoint

    def compareModel(self,DatabaseModel,Model):

        if len(DatabaseModel.attrs) == len(Model.Attributes)+1:
            for key,val in DatabaseModel.attrs.items():
                if key!=DATABASE_NUMBER:
                    try:
                        if DatabaseModel.attrs[key]!=Model.Attributes[key]:
                            return False
                    except:
                        return False
        else:
            return False

        return True


    def addModel(self,DatabasePoint,model):
        if not self.write:
            raise ValueError("No write intent")

        DatabaseModel=self.getModel(DatabasePoint,model)

        if DatabaseModel == -1:
            print("Add model to Point : {0}".format(DatabasePoint))
            if not self.write:
                raise ValueError("No write intent")
            DatabasePoint=self.getPoint(DatabasePoint)
            ModelNumber = self.getHighestNumber(DatabasePoint,DATABASE_MODEL)
            DatabaseModel=DatabasePoint.create_group('Model_'+str(ModelNumber+1))
            if type(model) is Model:
                for key,val in model.Attributes.items():
                    DatabaseModel.attrs[key]=np.string_(val)
            if type(model) is dict:
                for key,val in model.items():
                    DatabaseModel.attrs[key]=np.string_(val)

            DatabaseModel.attrs[DATABASE_NUMBER]=ModelNumber+1
            matrices = DatabaseModel.create_group(DATABASE_MODEL_MATRIX)
            modes    = DatabaseModel.create_group(DATABASE_MODEL_PODMODES)
            volfield = DatabaseModel.create_group(DATABASE_MODEL_VOLUMEFIELDS)
        else:
            del DatabaseModel[DATABASE_MODEL_MATRIX]
            del DatabaseModel[DATABASE_MODEL_PODMODES]
            del DatabaseModel[DATABASE_MODEL_VOLUMEFIELDS]

            if type(model) is Model:
                for key,val in model.Attributes.items():
                    DatabaseModel.attrs[key]=np.string_(val)
            if type(model) is dict:
                for key,val in model.items():
                    DatabaseModel.attrs[key]=np.string_(val)

            matrices = DatabaseModel.create_group(DATABASE_MODEL_MATRIX)
            modes    = DatabaseModel.create_group(DATABASE_MODEL_PODMODES)
            volfield = DatabaseModel.create_group(DATABASE_MODEL_VOLUMEFIELDS)
        return DatabaseModel


    def updateModelData(self,pointVal,ModelVal,MatrixDict,PODMODE=False,VOLUMEFIELD=False):
        point=self.getPoint(pointVal)
        if point==-1:
            return -1
        model=self.getModel(point,ModelVal)
        if model==-1:
            if type(ModelVal) is Model:
                self.addModel(point,ModelVal)
            else:
                print("can not add new model. Type {0} != Model {1}".format(type(ModelVal), Model))
                raise ValueError("can not add new model")
                return -1
        if not self.write:
            raise ValueError("No write intent")

        if not PODMODE and not VOLUMEFIELD:
            matrices = model[DATABASE_MODEL_MATRIX]
            ListOfNames= self.ListOfMatrixNames(model)
            it = len(matrices)+1
            for key,val in MatrixDict.items():
                ListOfNames= self.ListOfMatrixNames(model)
                if key in ListOfNames:
                    dataset=self.getMatrixByName(model,key)
                    dataset.resize(np.array(val).shape)
                    dataset[...]=val
                else:
                    array = np.array(val)
                    currentShape = array.shape
                    dimension = len(currentShape)
                    noneList = []                      
                    for i in range(0,dimension):
                        noneList.append(None)
                    if self.compression:                        
                        dataset=matrices.create_dataset('SystemMatrix_'+str(it), data=val, compression='gzip',maxshape=tuple(noneList))
                    else: 
                        dataset=matrices.create_dataset('SystemMatrix_'+str(it), data=val,maxshape=tuple(noneList))
                    dataset.attrs[DATABASE_NUMBER] = it
                    dataset.attrs[DATABASE_NAME]   = np.string_(key)
                    it+=1
            return model
        elif PODMODE and not VOLUMEFIELD:
            try:
                modes = model[DATABASE_MODEL_PODMODES]
            except:
                modes = model.create_group(DATABASE_MODEL_PODMODES)
            ListOfNames= self.ListOfModeNames(model)
            it = len(modes)+1
            for key,val in MatrixDict.items():
                ListOfNames= self.ListOfModeNames(model)
                if key in ListOfNames:
                    dataset=self.getModeByName(model,key)
                    dataset[...]=val.transpose()
                else:
                    if self.compression:
                      dataset=modes.create_dataset('Mode_'+str(it), data=val.transpose(), compression='gzip')
                    else:
                      dataset=modes.create_dataset('Mode_'+str(it), data=val.transpose(), chunks=True)
                    dataset.attrs[DATABASE_NUMBER] = it
                    dataset.attrs[DATABASE_NAME]   = np.string_(key)
                    it+=1
            return model
        elif VOLUMEFIELD and not PODMODE:
            try:
                modes = model[DATABASE_MODEL_VOLUMEFIELDS]
            except:
                modes = model.create_group(DATABASE_MODEL_VOLUMEFIELDS)
            ListOfNames= self.ListOfVolumeFieldNames(model)
            it = len(modes)+1
            for key,val in MatrixDict.items():
                ListOfNames= self.ListOfVolumeFieldNames(model)
                if key in ListOfNames:
                    dataset=self.getVolumeFieldByName(model,key)
                    dataset[...]=val.transpose()
                else:
                    if self.compression:
                      dataset=modes.create_dataset('VolumeField_'+str(it),data=val.transpose(), compression='gzip')
                    else:
                      dataset=modes.create_dataset('VolumeField_'+str(it),data=val.transpose())
                    dataset.attrs[DATABASE_NUMBER] = it
                    dataset.attrs[DATABASE_NAME]   = np.string_(key)
                    it+=1
                it+=1
            return model
        else:
            raise ValueError('(updateModelData): unclassified Type of Data')



    def extendModelData(self,pointVal,ModelVal,MatrixDict,PODMODE=False,VOLUMEFIELD=False):
        point=self.getPoint(pointVal)
        if point==-1:
            return -1
        model=self.getModel(point,ModelVal)
        if model==-1:
            if type(ModelVal) is Model:
                self.addModel(point,ModelVal)
            else:
                print("can not add new model. Type {0} != Model".format(type(ModelVal)))
                raise ValueError("can not add new model")
                return -1
        if not self.write:
            raise ValueError("No write intent")

        if not PODMODE and not VOLUMEFIELD:
            matrices = model[DATABASE_MODEL_MATRIX]
            ListOfNames= self.ListOfMatrixNames(model)
            for key in MatrixDict.keys():
                if key in ListOfNames:
                    print("Unable to create link (Name already exists)")
                    raise ValueError("Matrix name: "+ key +" already exists" )
                    return -1
            it = len(matrices)+1
            for key,val in MatrixDict.items():
                ListOfNames= self.ListOfMatrixNames(model)
                if key in ListOfNames:
                    print("Unable to create link (Name already exists)")
                    raise ValueError("Matrix name: "+ key +" already exists" )
                    return -1
                if self.compression:
                  dataset=matrices.create_dataset('SystemMatrix_'+str(it), data=val, compression='gzip')
                else:
                  dataset=matrices.create_dataset('SystemMatrix_'+str(it), data=val, chunks=True)
                dataset.attrs[DATABASE_NUMBER] = it
                dataset.attrs[DATABASE_NAME]   = np.string_(key)
                it+=1
            return model

        elif PODMODE and not VOLUMEFIELD:
            try:
                modes = model[DATABASE_MODEL_PODMODES]
            except:
                modes = model.create_group(DATABASE_MODEL_PODMODES)
            ListOfNames= self.ListOfModeNames(model)
            for key in MatrixDict.keys():
                if key in ListOfNames:
                    print("Unable to create link (Name already exists)")
                    raise ValueError("Mode name: "+ key +" already exists" )
                    return -1
            it = len(modes)+1
            for key,val in MatrixDict.items():
                ListOfNames= self.ListOfModeNames(model)
                if key in ListOfNames:
                    print("Unable to create link (Name already exists)")
                    raise ValueError("Matrix name: "+ key +" already exists" )
                    return -1
                if self.compression:
                  dataset=modes.create_dataset('Mode_'+str(it), data=val.transpose(), compression='gzip')
                else:
                  dataset=modes.create_dataset('Mode_'+str(it), data=val.transpose())
                dataset.attrs[DATABASE_NUMBER] = it
                dataset.attrs[DATABASE_NAME]   = np.string_(key)
                it+=1
            return model
        elif VOLUMEFIELD and not PODMODE:
            try:
                modes = model[DATABASE_MODEL_VOLUMEFIELDS]
            except:
                modes = model.create_group(DATABASE_MODEL_VOLUMEFIELDS)
            ListOfNames= self.ListOfVolumeFieldNames(model)

            for key in MatrixDict.keys():
                if key in ListOfNames:
                    print("Unable to create link (Name already exists)")
                    raise ValueError("Mode name: "+ key +" already exists" )
                    return -1
            it = len(modes)+1
            for key,val in MatrixDict.items():
                ListOfNames= self.ListOfVolumeFieldNames(model)
                if key in ListOfNames:
                    print("Unable to create link (Name already exists)")
                    raise ValueError("Matrix name: "+ key +" already exists" )
                    return -1
                if self.compression:
                  dataset=modes.create_dataset('VolumeField_'+str(it), data=val.transpose(), compression='gzip')
                else:
                  dataset=modes.create_dataset('VolumeField_'+str(it), data=val.transpose())
                dataset.attrs[DATABASE_NUMBER] = it
                dataset.attrs[DATABASE_NAME]   = np.string_(key)
                it+=1
            return model
        else:
            raise ValueError('(extendModelData): unclassified Type of Data')



    def ListOfMatrixNames(self,model):
        matrices = model[DATABASE_MODEL_MATRIX]
        NameList=[]
        for key in matrices.keys():
            NameList.append(model[DATABASE_MODEL_MATRIX][key].attrs[DATABASE_NAME])
        return NameList

    def getMatrixByName(self,model,name):
        matrices = model[DATABASE_MODEL_MATRIX]
        for key in matrices.keys():
            if model[DATABASE_MODEL_MATRIX][key].attrs[DATABASE_NAME]== name:
                return model[DATABASE_MODEL_MATRIX][key]
        return -1

    def getModeByName(self,model,name):
        matrices = model[DATABASE_MODEL_PODMODES]
        for key in matrices.keys():
            if model[DATABASE_MODEL_PODMODES][key].attrs[DATABASE_NAME]== name:
                return model[DATABASE_MODEL_PODMODES][key]
        return -1

    def ListOfModeNames(self,model):
        modes = model[DATABASE_MODEL_PODMODES]
        NameList=[]
        for key in modes.keys():
            NameList.append(model[DATABASE_MODEL_PODMODES][key].attrs[DATABASE_NAME])
        return NameList

    def getVolumeFieldByName(self,model,name):
        matrices = model[DATABASE_MODEL_VOLUMEFIELDS]
        for key in matrices.keys():
            if model[DATABASE_MODEL_VOLUMEFIELDS][key].attrs[DATABASE_NAME]== name:
                return model[DATABASE_MODEL_VOLUMEFIELDS][key]
        return -1

    def ListOfVolumeFieldNames(self,model):
        modes = model[DATABASE_MODEL_VOLUMEFIELDS]
        NameList=[]
        for key in modes.keys():
            NameList.append(model[DATABASE_MODEL_VOLUMEFIELDS][key].attrs[DATABASE_NAME])
        return NameList


    def getData(self,TYPE,pointVal,ModelVal,name=None):
        point=self.getPoint(pointVal)
        if point==-1:
            return -1
        model=self.getModel(point,ModelVal)
        if model==-1:
            return -1
        if type(name)!=str:
            if TYPE==DATABASE_MODEL_MATRIX:
                data=self.getSystemMatrices(model)
            elif TYPE==DATABASE_MODEL_PODMODES:
                data=self.getModes(model)
            elif TYPE==DATABASE_MODEL_VOLUMEFIELDS:
                data=self.getVolumeFields(model)
        else:
            if TYPE==DATABASE_MODEL_MATRIX:
                data=self.getSystemMatrix(model,name)
            elif TYPE==DATABASE_MODEL_PODMODES:
                data=self.getMode(model,name)
            elif TYPE==DATABASE_MODEL_VOLUMEFIELDS:
                data=self.getVolumeField(model,name)
        return data

    def getInternName(self,TYPE,pointVal,ModelVal,name):
        point=self.getPoint(pointVal)
        if point==-1:
            return -1
        model=self.getModel(point,ModelVal)
        if model==-1:
            return -1

        if TYPE==DATABASE_MODEL_MATRIX:
            for key in model[DATABASE_MODEL_MATRIX].keys():
                if model[DATABASE_MODEL_MATRIX][key].attrs[DATABASE_NAME]==name:
                    return key
        elif TYPE==DATABASE_MODEL_PODMODES:
            for key in model[DATABASE_MODEL_PODMODES].keys():
                if model[DATABASE_MODEL_PODMODES][key].attrs[DATABASE_NAME]==name:
                    return key
        elif TYPE==DATABASE_MODEL_VOLUMEFIELDS:
            for key in model[DATABASE_MODEL_VOLUMEFIELDS].keys():
                if model[DATABASE_MODEL_VOLUMEFIELDS][key].attrs[DATABASE_NAME]==name:
                    return key
        else:
            return -1


    def getSystemMatrices(self,DatabaseModel):
        ret={}
        for key in DatabaseModel[DATABASE_MODEL_MATRIX].keys():
            ret[DatabaseModel[DATABASE_MODEL_MATRIX][key].attrs[DATABASE_NAME]]=DatabaseModel[DATABASE_MODEL_MATRIX][key][...]
        return ret

    def getSystemMatrix(self,DatabaseModel,name):
        print(DatabaseModel[DATABASE_MODEL_MATRIX])
        for key in DatabaseModel[DATABASE_MODEL_MATRIX].keys():
            if DatabaseModel[DATABASE_MODEL_MATRIX][key].attrs[DATABASE_NAME]==name:
                return DatabaseModel[DATABASE_MODEL_MATRIX][key][...]
        print("Dataset: {0} can't be found in Model".format(name))
        return -1

    def getModes(self,DatabaseModel):
        ret={}
        for key in DatabaseModel[DATABASE_MODEL_PODMODES].keys():
            ret[DatabaseModel[DATABASE_MODEL_PODMODES][key].attrs[DATABASE_NAME]]=DatabaseModel[DATABASE_MODEL_PODMODES][key][...].transpose()
        return ret

    def getNumberOfModes(self,DatabaseModel):
        return len(DatabaseModel[DATABASE_MODEL_PODMODES].keys())

    def getMode(self,DatabaseModel,name):
        for key in DatabaseModel[DATABASE_MODEL_PODMODES].keys():

            if DatabaseModel[DATABASE_MODEL_PODMODES][key].attrs[DATABASE_NAME]==name:
                return DatabaseModel[DATABASE_MODEL_PODMODES][key][...].transpose()
        print("Dataset: {0} can't be found in Model".format(name))
        return -1

    def getVolumeFields(self,DatabaseModel):
        ret={}
        for key in DatabaseModel[DATABASE_MODEL_VOLUMEFIELDS].keys():
            ret[DatabaseModel[DATABASE_MODEL_VOLUMEFIELDS][key].attrs[DATABASE_NAME]]=DatabaseModel[DATABASE_MODEL_VOLUMEFIELDS][key][...].transpose()
        return ret

    def getVolumeField(self,DatabaseModel,name):
        for key in DatabaseModel[DATABASE_MODEL_VOLUMEFIELDS].keys():

            if DatabaseModel[DATABASE_MODEL_VOLUMEFIELDS][key].attrs[DATABASE_NAME]==name:
                return DatabaseModel[DATABASE_MODEL_VOLUMEFIELDS][key][...].transpose()
        print("Dataset: {0} can't be found in Model".format(name))
        return -1


    def meshing(self):
        try:
            FOdim = self.PIList[0].AttributeDimension
            if FOdim==2:
                self.meshingFiberOrientation()
            else:
                print('no mesh generation implemented for indicator {0}.'.format(self.PIList[0].AttributeName))
        except:
            print('no mesh generation implemented for indicator {0}.'.format(self.PIList[0].AttributeName))



    def meshingFiberOrientation(self):

        def checkArea(mtri):
            nTS=[]
            for T in mtri.triangles:
                A=np.ones((3,3))
                cnt=0
                for i in T:
                    A[cnt,1]=mtri.x[T[cnt]]
                    A[cnt,2]=mtri.y[T[cnt]]
                    cnt=cnt+1
                    F=0.5*np.linalg.det(A)
                    if F>1e-10:
                        nTS.append(T)
                        break
            mtri.triangles=np.array(nTS)

        mySimplice =[]
        if(self.getHighestNumber(self.DB,DATABASE_POINT)>=3):
            points=[]
            pointNumbers=[]
            pointTemperatures=[]
            for p in range(1,self.getHighestNumber(self.DB,DATABASE_POINT)+1):
                point=self.getPoint(p)
                pointNumbers.append(point.attrs[DATABASE_NUMBER])
                points.append(point.attrs[self.PIList[0].AttributeName])
                pointTemperatures.append(point.attrs[self.PIList[1].AttributeName])
                                    
            for d in set(pointTemperatures):
                points_d = []
                pointNumbers_d = []
                for p in range(len(points)):
                    if d == pointTemperatures[p]:
                        points_d.append(points[p])
                        pointNumbers_d.append(pointNumbers[p])
                
                if len(points_d) < 3:
                    continue
                
                nts=np.array(points_d)
                checkONEline = True  
                cnt = 2                
                px  = points_d[0][0]
                py  = points_d[0][1]
                qx  = points_d[1][0]
                qy  = points_d[1][1]

                dpqx = px - qx
                dpqy = py - qy
                while(checkONEline and cnt<len(points_d)):
                    mpx = points_d[cnt][0]
                    mpy = points_d[cnt][1]
                    dpmx = px - mpx
                    dpmy = py - mpy
                    cross = dpqx * dpmy - dpqy * dpmx
                    if(np.abs(cross)>1e-4):
                        checkONEline= False
                    cnt=cnt+1
                if(not checkONEline):
                    mtri=tri.Triangulation(nts[:,0],nts[:,1])
                    checkArea(mtri)
                    for s in mtri.triangles:
                        myS =[]
                        for pN in s:
                            myS.append(pointNumbers_d[pN])
                        mySimplice.append(copy.deepcopy(myS))
        try:
            self.DB[DATABASE_MESH]
        except:
            self.MeshIsAvailable =False

        if(self.MeshIsAvailable):
            del self.DB[DATABASE_MESH]
        if len(mySimplice)==0:
            mySimplice=[[-1,-1,-1]]
        meshgroup = self.DB.create_group(DATABASE_MESH)
        dataset= meshgroup.create_dataset(DATABASE_CONNECTIVITY, data=np.asarray(mySimplice))

    def getMesh(self):
        try:
            self.DB[DATABASE_MESH]
            self.MeshIsAvailable =True
        except:
            self.MeshIsAvailable =False
        if(self.MeshIsAvailable):
            return self.DB[DATABASE_MESH][DATABASE_CONNECTIVITY][...]
        else:
            print("Mesh can not be found")
            return -1

    def writeXDMF(self):
        itdim=0
        root = ET.Element('Xdmf')
        root.set('Version','2.2')

        domain=ET.SubElement(root,'Domain')

        grid=ET.SubElement(domain,'Grid')
        grid.set('GridType','Uniform')
        grid.set('Collection',str(itdim))

        Topology=ET.SubElement(grid,'Topology')
        Topology.set('TopologyType','3DCORECTMesh')

        geometry=ET.SubElement(grid,"Geometry")
        geometry.set('GeometryType','ORIGIN_DXDYDZ')

        DataItemOrigin=ET.SubElement(geometry,'DataItem')
        DataItemOrigin.set('Name','Origin')
        DataItemOrigin.set('Dimensions','3')
        DataItemOrigin.set('NumberType','Float')
        DataItemOrigin.set('Precision','4')
        DataItemOrigin.set('Format','XML')
        DataItemOrigin.text='0 0 0'

        DataItemSpacing=ET.SubElement(geometry,'DataItem')
        DataItemSpacing.set('Name','Spacing')
        DataItemSpacing.set('Dimensions','3')
        DataItemSpacing.set('NumberType','Float')
        DataItemSpacing.set('Precision','4')
        DataItemSpacing.set('Format','XML')
        DataItemSpacing.text='1 1 1'


        NumberOfPoints = self.getHighestNumber(self.DB,DATABASE_POINT)
        dims=[]

        for i in range(1,NumberOfPoints+1):
            point=self.getPoint(i)
            NumberOfModels = self.getHighestNumber(point,DATABASE_MODEL)
            for j in range(1,NumberOfModels+1):
                model=self.getModel(point,j)
                modes=self.getData(DATABASE_MODEL_VOLUMEFIELDS,i,j)
                it=0
                for key,val in modes.items():

                    it+=1
                    if len(val.shape)==4:
                        dimsNew=[val.shape[3],val.shape[2],val.shape[1]]
                    if len(val.shape)==3:
                        dimsNew=[val.shape[2],val.shape[1],val.shape[0]]

                    if len(dims)>0:
                        if dims !=dimsNew:
                            itdim+=1
                            grid=ET.SubElement(domain,'Grid')
                            grid.set('GridType','Uniform')
                            grid.set('Collection',str(itdim))

                            Topology=ET.SubElement(grid,'Topology')
                            Topology.set('TopologyType','3DCORECTMesh')
                            Topology.set('Dimensions',str(dimsNew[0]+1)+' '+str(dimsNew[1]+1)+' '+str(dimsNew[2]+1))

                            geometry=ET.SubElement(grid,"Geometry")
                            geometry.set('GeometryType','ORIGIN_DXDYDZ')

                            DataItemOrigin=ET.SubElement(geometry,'DataItem')
                            DataItemOrigin.set('Name','Origin')
                            DataItemOrigin.set('Dimensions','3')
                            DataItemOrigin.set('NumberType','Float')
                            DataItemOrigin.set('Precision','4')
                            DataItemOrigin.set('Format','XML')
                            DataItemOrigin.text='0 0 0'

                            DataItemSpacing=ET.SubElement(geometry,'DataItem')
                            DataItemSpacing.set('Name','Spacing')
                            DataItemSpacing.set('Dimensions','3')
                            DataItemSpacing.set('NumberType','Float')
                            DataItemSpacing.set('Precision','4')
                            DataItemSpacing.set('Format','XML')
                            DataItemSpacing.text='1 1 1'
                    else:
                        Topology.set('Dimensions',str(dimsNew[0]+1)+' '+str(dimsNew[1]+1)+' '+str(dimsNew[2]+1))




                    #Namestring='Point_'+str(i)+'/'+model.attrs[DATABASE_MODELTYPE]+'/'+key
                    Namestring='Point_'+str(i)+'/'
                    for attrKeys in model.attrs.keys():
                      if attrKeys != DATABASE_NUMBER:
                        Namestring += attrKeys + '_' + str(model.attrs[attrKeys])+'/'
                    Namestring += key

                    Attribute=ET.SubElement(grid,'Attribute')
                    Attribute.set('Name',Namestring)
                    Attribute.set('Center','Cell')
                    DataItem=ET.SubElement(Attribute,'DataItem')
                    if len(val.shape)==3:
                        Attribute.set('AttributeType','Scalar')
                        DataItem.set('Dimensions',str(val.shape[2])+' '+str(val.shape[1])+' '+str(val.shape[0]))
                    elif val.shape[0]==1:
                        Attribute.set('AttributeType','Scalar')
                        DataItem.set('Dimensions',str(val.shape[3])+' '+str(val.shape[2])+' '+str(val.shape[1]))
                    elif val.shape[0]==3:
                        Attribute.set('AttributeType','Vector')
                        DataItem.set('Dimensions',str(val.shape[3])+' '+str(val.shape[2])+' '+str(val.shape[1])+' '+str(val.shape[0]))
                    elif val.shape[0]==6:
                        Attribute.set('AttributeType','Tensor6')
                        DataItem.set('Dimensions',str(val.shape[3])+' '+str(val.shape[2])+' '+str(val.shape[1])+' '+str(val.shape[0]))
                    elif val.shape[0]==9:
                        Attribute.set('AttributeType','Tensor')
                        DataItem.set('Dimensions',str(val.shape[3])+' '+str(val.shape[2])+' '+str(val.shape[1])+' '+str(val.shape[0]))
                    else:
                        raise ValueError('unexpected mode shape')
                    DataItem.set('NumberType','Float')
                    DataItem.set('Format','HDF')
                    DataItem.text=os.path.basename(self.file)+':/Point_'+str(i)+'/Model_'+str(j)+'/'+DATABASE_MODEL_VOLUMEFIELDS+'/'+self.getInternName(DATABASE_MODEL_VOLUMEFIELDS,point,model,key)
                    dims =dimsNew

                modes=self.getData(DATABASE_MODEL_PODMODES,i,j)
                it=0
                for key,val in modes.items():

                    it+=1
                    if len(val.shape)==4:
                        dimsNew=[val.shape[3],val.shape[2],val.shape[1]]
                    if len(val.shape)==3:
                        dimsNew=[val.shape[2],val.shape[1],val.shape[0]]

                    if len(dims)>0:
                        if dims !=dimsNew:
                            itdim+=1
                            grid=ET.SubElement(domain,'Grid')
                            grid.set('GridType','Uniform')
                            grid.set('Collection',str(itdim))

                            Topology=ET.SubElement(grid,'Topology')
                            Topology.set('TopologyType','3DCORECTMesh')
                            Topology.set('Dimensions',str(dimsNew[0]+1)+' '+str(dimsNew[1]+1)+' '+str(dimsNew[2]+1))

                            geometry=ET.SubElement(grid,"Geometry")
                            geometry.set('GeometryType','ORIGIN_DXDYDZ')

                            DataItemOrigin=ET.SubElement(geometry,'DataItem')
                            DataItemOrigin.set('Name','Origin')
                            DataItemOrigin.set('Dimensions','3')
                            DataItemOrigin.set('NumberType','Float')
                            DataItemOrigin.set('Precision','4')
                            DataItemOrigin.set('Format','XML')
                            DataItemOrigin.text='0 0 0'

                            DataItemSpacing=ET.SubElement(geometry,'DataItem')
                            DataItemSpacing.set('Name','Spacing')
                            DataItemSpacing.set('Dimensions','3')
                            DataItemSpacing.set('NumberType','Float')
                            DataItemSpacing.set('Precision','4')
                            DataItemSpacing.set('Format','XML')
                            DataItemSpacing.text='1 1 1'
                    else:
                        Topology.set('Dimensions',str(dimsNew[0]+1)+' '+str(dimsNew[1]+1)+' '+str(dimsNew[2]+1))




                    #Namestring='Point_'+str(i)+'/'+model.attrs[DATABASE_MODELTYPE]+'/'+key
                    Namestring='Point_'+str(i)+'/'
                    for attrKeys in model.attrs.keys():
                      if attrKeys != DATABASE_NUMBER:
                        Namestring += attrKeys + '_' + str(model.attrs[attrKeys])+'/'
                    Namestring += key

                    Attribute=ET.SubElement(grid,'Attribute')
                    Attribute.set('Name',Namestring)
                    Attribute.set('Center','Cell')
                    DataItem=ET.SubElement(Attribute,'DataItem')
                    if len(val.shape)==3:
                        Attribute.set('AttributeType','Scalar')
                        DataItem.set('Dimensions',str(val.shape[2])+' '+str(val.shape[1])+' '+str(val.shape[0]))
                    elif val.shape[0]==1:
                        Attribute.set('AttributeType','Scalar')
                        DataItem.set('Dimensions',str(val.shape[3])+' '+str(val.shape[2])+' '+str(val.shape[1]))
                    elif val.shape[0]==3:
                        Attribute.set('AttributeType','Vector')
                        DataItem.set('Dimensions',str(val.shape[3])+' '+str(val.shape[2])+' '+str(val.shape[1])+' '+str(val.shape[0]))
                    elif val.shape[0]==6:
                        Attribute.set('AttributeType','Tensor6')
                        DataItem.set('Dimensions',str(val.shape[3])+' '+str(val.shape[2])+' '+str(val.shape[1])+' '+str(val.shape[0]))
                    elif val.shape[0]==9:
                        Attribute.set('AttributeType','Tensor')
                        DataItem.set('Dimensions',str(val.shape[3])+' '+str(val.shape[2])+' '+str(val.shape[1])+' '+str(val.shape[0]))
                    else:
                        print(val.shape)
                        raise ValueError('unexpected mode shape')
                    DataItem.set('NumberType','Float')
                    DataItem.set('Format','HDF')
                    DataItem.text=os.path.basename(self.file)+':/Point_'+str(i)+'/Model_'+str(j)+'/'+DATABASE_MODEL_PODMODES+'/'+self.getInternName(DATABASE_MODEL_PODMODES,point,model,key)
                    dims =dimsNew


        indent(root)
        tree = ET.ElementTree(root)
        tree.write(self.file + ".xdmf", encoding="utf-8")


#PIFO=PointIndicator(DATABASE_FIBERORIENTATION,2,1e-4)

#DB=Database("/m/scratch/itwm/koebler/XDMF_J/Database.h5")


##-------------------------------------------------------------
##write
#
##Add Point
#
##define database indicator (id)
#PIFO=PointIndicator(DATABASE_FIBERORIENTATION,2,1e-4) #attribute name, attr. dimension, tolerance
#PIVF=PointIndicator(DATABASE_FIBERVOLUMEFRACTION,1,1e-2)
#
##open database
#DB=Database("/m/bulk/sms/opt2/users/koebler/t.h5",new=True,PointIndicator=[PIFO,PIVF]) #filename, indicator, write access
#
##add point to database
#DatabasePoint=DB.addPoint({DATABASE_FIBERORIENTATION:[0.8,0.1]})
#print DatabasePoint
#
##define model attr
#ModelAttr={'Type':'LinearElastic'}
#
##define model data
#ModelData={'Matrix1': np.random.rand(3,3,3,7),'Matrix1': np.random.rand(6,2,2,2)}
#
##add model to point in the database
#MyModel=DB.addModel(DatabasePoint,ModelAttr)
#
##define Mode data
#ModeData={'mode1':np.random.rand(6,4,4,4),'mode2':np.random.rand(6,4,4,4)}
#
##define VolumeFieldData data
#VolumeFieldData={'myfirstvolumefield':np.ones((3,1,2,3)),'mysecondvolumeField':np.ones((6,16,16,16))}
#VolumeFieldData['myfirstvolumefield'][0,:,:,:]*=1.0
#VolumeFieldData['myfirstvolumefield'][1,:,:,:]*=2.0
#VolumeFieldData['myfirstvolumefield'][2,:,:,:]*=3.0
#
##extend model --> write data
#DB.extendModelData({DATABASE_FIBERORIENTATION:[0.8,0.1],DATABASE_FIBERVOLUMEFRACTION:1.0},{DATABASE_MODELTYPE:'LinearElastic'},ModelData)
#DB.extendModelData({DATABASE_FIBERORIENTATION:[0.8,0.1],DATABASE_FIBERVOLUMEFRACTION:1.0},{DATABASE_MODELTYPE:'LinearElastic'},ModeData,PODMODE=True)

#DB.extendModelData({DATABASE_FIBERORIENTATION:[0.8,0.1]},{DATABASE_MODELTYPE:'LinearElastic'},VolumeFieldData,VOLUMEFIELD=True)

#DB.updateModelData({DATABASE_FIBERORIENTATION:[0.8,0.10],DATABASE_FIBERVOLUMEFRACTION:0.991},{DATABASE_MODELTYPE:'LinearElastic'},ModelData)
#DB.updateModelData({DATABASE_FIBERORIENTATION:[0.8,0.1],DATABASE_FIBERVOLUMEFRACTION:1.0},{DATABASE_MODELTYPE:'LinearElastic'},ModeData,PODMODE=True)
#DB.updateModelData({DATABASE_FIBERORIENTATION:[0.8,0.1],DATABASE_FIBERVOLUMEFRACTION:1.0},{DATABASE_MODELTYPE:'LinearElastic'},VolumeFieldData,VOLUMEFIELD=True)
#
##write xdmf file --> paraview visualization
#DB.writeXDMF("/m/scratch/itwm/koebler/XDMF_J/t.xdmf")
#
##close database
#DB.close()
#
##-------------------------------------------------------------

##-------------------------------------------------------------
##read
#
##open database
#DB=Database("/m/bulk/sms/opt2/users/koebler/test.h5") #filename
#
##read single matrix
#Matrix1=DB.getData(DATABASE_MODEL_MATRIX,1,1,'Matrix1')# point 1, model 1
#
##or use
#Matrix1=DB.getData(DATABASE_MODEL_MATRIX,{DATABASE_FIBERORIENTATION:[0.8,0.1],DATABASE_FIBERVOLUMEFRACTION:1.0},{DATABASE_MODELTYPE:'LinearElastic'},'Matrix1')
#
##read all matrices
#
#DictContainMatrices= DB.getData(DATABASE_MODEL_MATRIX,1,1)
#
##read single PODMODE
#mode1=DB.getData(DATABASE_MODEL_PODMODES,{DATABASE_FIBERORIENTATION:[0.8,0.1],DATABASE_FIBERVOLUMEFRACTION:1.0},{DATABASE_MODELTYPE:'LinearElastic'},'mode1')
#
##close database
#DB.close()
#
##-------------------------------------------------------------

##-------------------------------------------------------------
##check
#
##open database
#DB=Database("/m/bulk/sms/opt2/users/koebler/t.h5") #filename, indicator, write access
#
##print attrs
#DB.DB.visititems(print_attrs)
#
##hdf5dump
#h5dump  /m/bulk/sms/opt2/users/koebler/test.h5
#
##close database
#
#DB.close()
#
##-------------------------------------------------------------
