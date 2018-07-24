# Biomass Supply Chain Network Data Extractor
# Author name: Nathanial Cooper
# Author email: nathanial.cooper AT imperial.ac.uk (preferred), nattiecooper AT gmail.com (alternate)
# Date: 24 July 2018
# Version: 1.0
#-------------------------------------------DESCRIPTION-----------------------------------------------------------
# This program extracts data relevant for biomass supply chain network optimization.
# This script requires: 
# 1) a raster with data about the land use (for example the CORINE Lane Use database raster)
# 2) the no data value of the raster
# 3) the value of the land use of interest for the land use raster
# 4) a raster of constraints on the usability of the land
# 5) a raster of the biomass yield values (up to 5)
# 6) the biomass name used in the optimisation associated with the yield raster
# 7) a shapefile or vector layer of some sort that represents the area of interest.
# 8) the desired number of grid cells
# 9) the minimum cell area fraction to be counted
# This script outputs: 
# 1) a CSV of the ID numbers of all cells, the area of all the cells 
# 2) a CSV of all linear distances between the centroids of all cells
# 3) a CSV of the area of the raster bands of interest in each cell
# 4) a CSV of the tortuosity matrix between cells
# 5) a CSV of the average biomass yield per cell
# 6) CSV of the slope and intercepts  for linear estimators of the biomass yield per cell
# 7) a shapefile of the grid used
#----------------------------------------------------------------------------------------------------------------------

##Supply Chain Network Data Aggregator and Linear Estimator=name
##Land_Use_Raster_Data=raster
##No_Data_Value=number 48
##Raster_Value_of_Interest_separate_values_by_commas=string 12
##Constraint_Layer=optional raster
##Bioyield_Map_A=optional raster
##Bioyield_Product_A=optional string Product_A
##Bioyield_Map_B=optional raster
##Bioyield_Product_B=optional string Product_B
##Bioyield_Map_C=optional raster
##Bioyield_Product_C=optional string Product_C
##Bioyield_Map_D=optional raster
##Bioyield_Product_D=optional string Product_D
##Bioyield_Map_E=optional raster
##Bioyield_Product_E=optional string Product_E
##Region_Shapefile=vector Polygon
##Desired_Number_of_Cells=number 100
##Only_count_cells_with_at_least_x_percent_of_the_area_of_the_largest_cell_where_0_includes_all_cells=number 5
##Lower_Breakpoint_Fraction=number 0.333
##Upper_Breakpoint_Fraction=number 0.666
##Output_Cell_Area=output table
##Output_Linear_Distance_Between_Cells=output table
##Output_Arable_Land_Area=output table
##Output_Tortuosity_Between_Cells=output table
##Output_Bioyield_Averages=output table
##Output_Bioyield_Low_Approximation_m=output table
##Output_Bioyield_Mid_Approximation_m=output table
##Output_Bioyield_High_Approximation_m=output table
##Output_Bioyield_Low_Approximation_b=output table
##Output_Bioyield_Mid_Approximation_b=output table
##Output_Bioyield_High_Approximation_b=output table
##Output_Grid_SHP=output vector
# ##Output_yieldmap1=output table
# ##Output_yieldmap2=output table

# Load libraries from QGIS and Python core
from qgis.core import *
from qgis.utils import *
from PyQt4.QtCore import *
import math
import csv
import time
import bisect


# ------------------------Variable Declarations and Initializations
t_0 = time.time()
numGridCells = max(Desired_Number_of_Cells,1) # approximate number of cells that will be created in the grid
minCellFrac = float(Only_count_cells_with_at_least_x_percent_of_the_area_of_the_largest_cell_where_0_includes_all_cells)/100 # the minimum cell size as a fraction of the largest cell size
noDataVal = No_Data_Value
tmpNum = 0
rastNumInterest = [int(s) for s in Raster_Value_of_Interest_separate_values_by_commas.split(',')]
outATcsv = Output_Cell_Area
outLDcsv = Output_Linear_Distance_Between_Cells
outLAcsv = Output_Arable_Land_Area
outTTcsv = Output_Tortuosity_Between_Cells
outBYcsv = Output_Bioyield_Averages
outBYLUmcsv = Output_Bioyield_Low_Approximation_m
outBYMUmcsv = Output_Bioyield_Mid_Approximation_m
outBYHUmcsv = Output_Bioyield_High_Approximation_m
outBYLUbcsv = Output_Bioyield_Low_Approximation_b
outBYMUbcsv = Output_Bioyield_Mid_Approximation_b
outBYHUbcsv = Output_Bioyield_High_Approximation_b
outSHP = Output_Grid_SHP
maxIter = 10 # maximum number of iterations to get correct number of cells
minUnder = 0.95
maxOver = 1.05
iter = 0
areaConv = 100/4 #converts hectares to km^2 and then annual production to monthly production
x_lo = Lower_Breakpoint_Fraction
x_hi = Upper_Breakpoint_Fraction
#yieldmapcount = 0

# new attribute names
fieldNameID = "IDnum" 
fieldNameArea = "Area-sq_km"
fieldNameX = "Centroid X"
fieldNameY = "Centroid Y"
fieldNameLen = "Road Len"
fieldNameCount = "Road Count"

# Copy the layers and leave the originals intact
fullRast_noConstr = processing.getObject(Land_Use_Raster_Data)
countrySHP = processing.getObject(Region_Shapefile)
constRast = processing.getObject(Constraint_Layer)
byARast = processing.getObject(Bioyield_Map_A)
byBRast = processing.getObject(Bioyield_Map_B)
byCRast = processing.getObject(Bioyield_Map_C)
byDRast = processing.getObject(Bioyield_Map_D)
byERast = processing.getObject(Bioyield_Map_E)


# ------------------------Reproject shapefiles to the CRS of the raster
def layerReproj(layer_unknownCRS, layer_destCRS):
    if not(layer_destCRS.crs().authid() == layer_unknownCRS.crs().authid()): # check to see if the CRS are different
        iface.mainWindow().statusBar().showMessage("Checking the layer CRS") # Update the user in the status bar
        tmpReproj = processing.runalg("qgis:reprojectlayer", layer_unknownCRS, layer_destCRS.crs().authid(),None) # reroject the country shapefile to the new CRS
        layer_newCRS = QgsVectorLayer(tmpReproj['OUTPUT'], "Reprojected Shapefile", "ogr") # add the reprojection to the project
        if layer_newCRS.isValid():
            iface.mainWindow().statusBar().showMessage("Reprojected the country shp CRS: %s" % layer_newCRS.crs().authid()) # Update the user in the status bar
        else:
            iface.messageBar().pushInfo("Error","Could not reproject vector layer, please check vector layer", level=QgsMessageBar.CRITICAL)
    else: # if the CRS are the same, we are done
        layer_newCRS = layer_unknownCRS
        iface.mainWindow().statusBar().showMessage("Shapefile already has correct CRS") # Update the user in the status bar
    return layer_newCRS

# ------------------------Reproject rasters to the CRS of the raster
def rastReproj(layer_unknownCRS, layer_destCRS,floatFlag):
    canvExtent = iface.mapCanvas().extent()
    exmin = canvExtent.xMinimum()
    exmax = canvExtent.xMaximum()
    eymin = canvExtent.yMinimum()
    eymax = canvExtent.yMaximum()
    extent = "%f,%f,%f,%f" %(exmin, exmax, eymin, eymax)
    
    if not(layer_destCRS.crs().authid() == layer_unknownCRS.crs().authid()): # check to see if the CRS are different
        iface.mainWindow().statusBar().showMessage("Checking the layer CRS") # Update the user in the status bar
        tmpReproj = processing.runalg("gdalogr:warpreproject",layer_unknownCRS,layer_unknownCRS.crs().authid(),layer_destCRS.crs().authid(),noDataVal,0.0,0,extent,'',(4+floatFlag),3,75,6,1,False,2,False,'',None)
        layer_newCRS = QgsRasterLayer(tmpReproj['OUTPUT'], "Reprojected Raster") # add the reprojection to the project
        if layer_newCRS.isValid():
            iface.mainWindow().statusBar().showMessage("Reprojected the country shp CRS: %s" % layer_newCRS.crs().authid()) # Update the user in the status bar
        else:
            iface.messageBar().pushInfo("Error","Could not reproject vector layer, please check vector layer", level=QgsMessageBar.CRITICAL)
    else: # if the CRS are the same, we are done
        layer_newCRS = layer_unknownCRS
        iface.mainWindow().statusBar().showMessage("Shapefile already has correct CRS") # Update the user in the status bar
    return layer_newCRS

# ------------------------Create a tight grid overlay for  the country
def gridGen(layer_countryProfile, cellsize):
    if (not layer_countryProfile.isValid()):
        return layer_countryProfile
    
    # Generate the overall grid
    iface.mainWindow().statusBar().showMessage("Generating Grid Overlay") # Update the user in the status bar
    # define the extent of the grid, based on the size of the country shapefile
    extent = "%f, %f, %f, %f" % (layer_countryProfile.extent().xMinimum(), layer_countryProfile.extent().xMaximum(), layer_countryProfile.extent().yMinimum(), layer_countryProfile.extent().yMaximum()) 
    tmpGrid = processing.runalg('qgis:vectorgrid', extent, cellsize, cellsize, 0, None) # create the grid and save it
    layer_newGrid = QgsVectorLayer(tmpGrid['OUTPUT'], "Full Grid", "ogr") # add the grid to the project
    if layer_newGrid.isValid():
        iface.mainWindow().statusBar().showMessage("Generating Grid Overlay - Complete") # Update the user in the status bar
    else:
        iface.messageBar().pushInfo("Error","Trouble generating grid - vectorgrid", level=QgsMessageBar.WARNING)
    
    # Tighten the grid to just the country
    iface.mainWindow().statusBar().showMessage("Clipping Grid to country profile") # Update the user in the status bar
    tmpFull = processing.runalg('qgis:clip',layer_newGrid, layer_countryProfile,None) # create the tight grid by clipping the grid with the country profile
    layer_fullGrid = QgsVectorLayer(tmpFull['OUTPUT'], "Full Grid", "ogr") # add the tight grid to the project
    if layer_fullGrid.isValid():
        iface.mainWindow().statusBar().showMessage("Clipping grid to country profile - complete") # Update the user in the status bar
    else:
        iface.messageBar().pushInfo("Error","Trouble generating grid - clip", level=QgsMessageBar.WARNING)
    
    return layer_fullGrid


# ------------------------Delete old fields. Add field for ID, and field for area. Add area, x and y coordinates to vectors for each
def replaceFields(layer_repFields):
    iface.mainWindow().statusBar().showMessage("Adding cell ID and area to vector layer") # Update the user in the status bar
    caps = layer_repFields.dataProvider().capabilities() # need to have access to layer capabilities in order to update fields
    if not caps: # if the capabilities are not accessible, this a problem.
        iface.messageBar().pushInfo("Error","loading layer capabilities failed", level=QgsMessageBar.WARNING) 
    elif ((not QgsVectorDataProvider.DeleteAttributes) and (not QgsVectorDataProvider.AddAttributes)): # see if attribute methods are available. if they are not, this is a problem
        iface.messageBar().pushInfo("Error","loading data provider methods failed", level=QgsMessageBar.WARNING)
    else: # if these are available, we may proceed
        fieldCount = 0 # set the field count at 0, then count the number of current fields (info needed to delete them)
        for field in layer_repFields.pendingFields():
            fieldCount += 1
        delFields = range(fieldCount) 
        layer_repFields.dataProvider().deleteAttributes(delFields) # delete all current attributes
        layer_repFields.updateFields()
        # add the attributes that we care about
        layer_repFields.dataProvider().addAttributes([QgsField(fieldNameID, QVariant.Int), QgsField(fieldNameArea, QVariant.Double)]) 
        layer_repFields.updateFields()
        
        layer_repFields.startEditing() # initiate editing of the layer attributes
        iter = layer_repFields.getFeatures()
        areaOp = QgsDistanceArea()
        areaOp.computeAreaInit()
        for feature in iter: # loop across all features in the layer
            if feature.geometry().wkbType() == QGis.WKBPolygon: # check to make sure the geometry is a valid choice
                feature[0] = feature.id()+1 # set the new feature ID value
                feature[1] = areaOp.measurePolygon(feature.geometry().asPolygon()[0])/(1000*1000)
                layer_repFields.updateFeature(feature) # update the feature with the field changes
            elif feature.geometry().wkbType() == QGis.WKBMultiPolygon: # check to make sure the geometry is a valid choice
                feature[0] = feature.id()+1 # set the new feature ID value
                inst = feature.geometry().asMultiPolygon()
                instArea = 0
                for geom in inst:
                    instArea = instArea + areaOp.measurePolygon(geom[0])/(1000*1000)
                feature[1] = instArea
                layer_repFields.updateFeature(feature) # update the feature with the field changes
            else: #can't / don't want to calculate for other geometry things
                iface.messageBar().pushInfo("Error","Unknown Geometry", level=QgsMessageBar.WARNING) # if the feature is not a polygon, this is a problem
        iface.mainWindow().statusBar().showMessage("Adding area and centroid location to vector layer - Complete") # Update the user in the status bar
        layer_repFields.commitChanges() # commit these changes 
    
    return layer_repFields


# ------------------------Number of grid cells that are too small
def countSmallFeatures(layer_smallFeat):
    if (not layer_smallFeat.isValid()):
        return -1
    
    minArea = minCellFrac*layer_smallFeat.maximumValue(layer_smallFeat.fieldNameIndex(fieldNameArea))
    strTooSmall = "\"%s\" < %f" % (fieldNameArea,minArea)
    requestSmallFeat = QgsFeatureRequest().setFilterExpression(strTooSmall)
    requestSmallFeat.setFlags(QgsFeatureRequest.NoGeometry)
    smallFeat = layer_smallFeat.getFeatures(requestSmallFeat)

    return smallFeat

# ------------------------Write attribute table to CSV
def writeAttributeTable(layer_attr):
    iface.mainWindow().statusBar().showMessage("Writing attribute table to CSV") # Update the user in the status bar
    if (not layer_attr.isValid()):
        return -1
        
    # will only keep the file open during the with statement, then closes
    with open(outATcsv,"wb") as ofileAT:
        writerAT = csv.writer(ofileAT, delimiter=',', quoting=csv.QUOTE_NONNUMERIC) # writer is necessary to create csv
    
        # loop across all features in the layer
        iter = layer_attr.getFeatures()
        for feature in iter:
            # check to make sure the geometry is a valid choice
            if feature.geometry().type() == QGis.Polygon:
                writerAT.writerow(feature.attributes()) # write all feature attributes to the csv
            else:
                iface.messageBar().pushInfo("Error","Unknown Geometry", level=QgsMessageBar.WARNING) # if the feature is not a polygon, this is a problem
        iface.mainWindow().statusBar().showMessage("Writing attribute table to CSV - complete") # Update the user in the status bar
    
    return 1

# ------------------------Write vector to file
def writeVec(vect,destFile):
    iface.mainWindow().statusBar().showMessage("Writing vector to CSV") # Update the user in the status bar
    if (not vect):
        return -1
        
    with open(destFile,"wb") as ofileDest:
        writerDest = csv.writer(ofileDest, delimiter=',', quoting=csv.QUOTE_NONNUMERIC) # writer is necessary to create csv
        
        isList = sum([isinstance(x,list) for x in vect])
        maxLen = 0
        
        if isList > 0:
            maxLen = max([len(x) for x in vect])+1
        for k in range(0,len(vect)):
            if isinstance(vect[k],list):
                row_k = [k+1] + vect[k]
            else:
                row_k = [k+1] + [vect[k]]
            row_k = row_k + [0]*(maxLen-len(row_k))
            writerDest.writerow(row_k)
        
    iface.mainWindow().statusBar().showMessage("Writing vector to CSV - Complete") # Update the user in the status bar
    
    return 1

# ------------------------Calulate linear distance between all centroid locations
def writeLinDist(layer_grid):
    iface.mainWindow().statusBar().showMessage("Writing linear distances to CSV") # Update the user in the status bar
    if (not layer_grid.isValid()):
        return -1
    
    # will only keep the file open during the with statement, then closes
    with open(outLDcsv,"wb") as ofileLD:
        writerLD = csv.writer(ofileLD, delimiter=',', quoting=csv.QUOTE_NONNUMERIC) # writer is necessary to create csv
        headerRow = (range(0,len(list(layer_grid.getFeatures()))+1)) # top line that lists all cell ID
        writerLD.writerow(headerRow)
        
        meas = QgsDistanceArea()
        meas.setEllipsoid(layer_grid.crs().authid())
        meas.setEllipsoidalMode(True)
        
        iterV = layer_grid.getFeatures()
        for featureV in iterV : # iterate over each cell, row-wise
            row_i = [featureV[0]] # starts the row with the ID of the current cell
            iterH = layer_grid.getFeatures()
            for featureH in iterH: # iterate over each cell, column-wise
                if (featureH[0] == featureV[0]):
                    dist = 0
                elif featureV.geometry().touches(featureH.geometry()):
                    dist = meas.measureLine(featureV.geometry().centroid().asPoint(),featureH.geometry().centroid().asPoint())/1000
                else:
                    dist = 999
                row_i.append(dist) # find the linear distance in km and append it to the current vector of distances
            writerLD.writerow(row_i) # write the vector of distances to the csv
            
    iface.mainWindow().statusBar().showMessage("Writing linear distances to CSV - Complete") # Update the user in the status bar
    
    return 1


# ------------------------Calculate Tortuosity Data
def writeTortuosity(layer_tort):
    if(not(layer_tort.isValid())):
        return -1

    iface.mainWindow().statusBar().showMessage("Writing Tortuosities to CSV")
    with open(outTTcsv,"wb") as ofileTT:
        writerTT = csv.writer(ofileTT, delimiter=',', quoting=csv.QUOTE_NONNUMERIC) # writer is necessary to create csv
        headerRow = range(0,len(list(layer_tort.getFeatures()))+1) # top line that lists all cell ID
        headerRow.insert(0,0)
        writerTT.writerow(headerRow)
        
        iterV = layer_tort.getFeatures()
            
        for featureV in iterV:
            row_i = [featureV[0]]
            row_i.append("truck")
            iterH = layer_tort.getFeatures()
                
            for featureH in iterH:
                if (featureH[0] == featureV[0]):
                    row_i.append(0)
                elif featureV.geometry().touches(featureH.geometry()):
                    row_i.append(1.4)
                else:
                    row_i.append(2.6)
            writerTT.writerow(row_i)
    iface.mainWindow().statusBar().showMessage("Writing linear distances to CSV - Complete") # Update the user in the status bar
    
    return 1


# ------------------------convert grid to raster w/ ID - rasterize
def gridToRast(layer_grid):
    if (not layer_grid.isValid()):
        return -1
        
    fGxmin = layer_grid.extent().xMinimum()
    fGxmax = layer_grid.extent().xMaximum()
    fGymin = layer_grid.extent().yMinimum()
    fGymax = layer_grid.extent().yMaximum()
    extent = "%f,%f,%f,%f" %(fGxmin, fGxmax, fGymin, fGymax)
    
    iface.mainWindow().statusBar().showMessage("Generating raster from country profile") # Update the user in the status bar
    
    tmpRast = processing.runalg("gdalogr:rasterize",layer_grid, fieldNameID, 1, 100, 100, extent, 0, 1, -999, 3, 75, 6, 1, 0, 2, '', None)
    layer_rast = QgsRasterLayer(tmpRast['OUTPUT'], "Full Grid Raster")
    
    if layer_rast.isValid():
        iface.mainWindow().statusBar().showMessage("Generating raster from country profile - Complete") # Update the user in the status bar
    else:
        iface.messageBar().pushInfo("Error","Trouble generating graster from country profile", level=QgsMessageBar.WARNING)
        
    return layer_rast


# ------------------------Perform Statistics on the raster, and write it into a temporary vector
def gridStats(layer_gridRast, layer_fullRast, layer_gridSource):
    if ((not layer_gridRast.isValid()) or (not layer_fullRast.isValid())):
        return -1
        
    fGxmin = layer_gridSource.extent().xMinimum()
    fGxmax = layer_gridSource.extent().xMaximum()
    fGymin = layer_gridSource.extent().yMinimum()
    fGymax = layer_gridSource.extent().yMaximum()
    extent = "%f,%f,%f,%f" %(fGxmin, fGxmax, fGymin, fGymax)

    iface.mainWindow().statusBar().showMessage("Performing statistics on Land Area Usage") # Update the user in the status bar
    rastPath = "%s;%s" %(layer_gridRast.dataProvider().dataSourceUri(), layer_fullRast.dataProvider().dataSourceUri())
    tmpVar = processing.runalg("grass7:r.stats",rastPath,"space","-999","255",False,False,True,False,False,False,False,False,False,True,True,False,False,extent,None,None)
    rastAreas = open(tmpVar['rawoutput'])

    # read raster stats data into vector
    currLine = rastAreas.readline()
    vector_land = [0]*len(list(layer_gridSource.getFeatures()))
    while currLine != "":
        cID,cLandUse,cArea = currLine.split()
        if int(cLandUse) in rastNumInterest:
            vector_land[int(cID)-1] = vector_land[int(cID)-1] + float(cArea)
        currLine = rastAreas.readline()
    iface.mainWindow().statusBar().showMessage("Performing statistics on Land Area Usage - Complete") # Update the user in the status bar
    
    return vector_land
    
# ------------------------Perform bioyield raster calculations
def bioyield(layer_constraint, layer_bio,layer_fullRast,layer_gridSource,layer_gridRast):
    if (not(layer_constraint) or not(layer_bio) or not(layer_constraint.isValid()) or not(layer_bio.isValid()) or not(layer_fullRast.isValid())):
        print "Something is wrong!"
        return -1
        
    fGxmin = layer_gridSource.extent().xMinimum()
    fGxmax = layer_gridSource.extent().xMaximum()
    fGymin = layer_gridSource.extent().yMinimum()
    fGymax = layer_gridSource.extent().yMaximum()
    extent = "%f,%f,%f,%f" %(fGxmin, fGxmax, fGymin, fGymax)
    
    iface.mainWindow().statusBar().showMessage("Performing Bioyield Raster Calculations") # Update the user in the status bar
    if len(rastNumInterest) == 1:
        eqn = '(A==%s)*B' % str(rastNumInterest[0])
    elif len(rastNumInterest) > 1:
        eqn = '((A==%s)' % str(rastNumInterest[0])
        for i in range(1,len(rastNumInterest)):
            eqn = eqn + ('+ (A==%s)' % str(rastNumInterest[i]))
        eqn = eqn + ')*B'
    output_tmpRast = processing.runalg("grass:r.mapcalculator",layer_fullRast.dataProvider().dataSourceUri(),layer_bio.dataProvider().dataSourceUri(),None,None,None,None,eqn,extent,0,None)
    bioRast = processing.getObject(output_tmpRast['outfile'])
    
    rastPath = "%s;%s" %(layer_gridRast.dataProvider().dataSourceUri(), bioRast.dataProvider().dataSourceUri())
    tmpVar = processing.runalg("grass7:r.stats",rastPath,"space","-999","255",False,True,True,False,False,False,False,False,False,True,True,False,False,extent,None,None)
    rastAreas = open(tmpVar['rawoutput'])
    
    currLine = rastAreas.readline()
    vector_BYland = [[] for i in range(len(list(layer_gridSource.getFeatures())))]
    vector_BYlandAvg = [[] for i in range(len(list(layer_gridSource.getFeatures())))]
    while currLine != "":
        cID,cBY,cArea = currLine.split()
        vector_BYland[int(cID)-1].append(float(cBY)*areaConv) # [tdm / (km^2 * months)]
        vector_BYland[int(cID)-1].append(float(cArea)/(1000.0*1000.0)) # [km^2]
        currLine = rastAreas.readline()
    iface.mainWindow().statusBar().showMessage("Performing Bioyield Raster Calculations - Complete") # Update the user in the status bar
    
#    global yieldmapcount
#    if (yieldmapcount == 0):
#        tmpflag = writeVec(vector_BYland,Output_yieldmap1)
#        yieldmapcount = 1
#    elif (yieldmapcount == 1):
#        tmpflag = writeVec(vector_BYland,Output_yieldmap2)
#        yieldmapcount = 2
    
    # Calculate the average bioyield for each cell
    for i in range(len(vector_BYland)):
        totArea = 0
        totGrowth = 0
        for j in range(len(vector_BYland[i])/2):
            totArea = totArea+vector_BYland[i][j*2+1] # [km^2]
            totGrowth = totGrowth+(vector_BYland[i][j*2+1]*vector_BYland[i][j*2]) # [tdm / month]
        if totGrowth>0:
            vector_BYlandAvg[i] = totArea/totGrowth # [km^2 month / tdm]
        else:
            vector_BYlandAvg[i] = 0
            
    # Calculate the total land area used and the total biomass produced by summing, going from best land to worst land
    landArea = [[] for i in range(len(vector_BYland))]
    bioPerYR = [[] for i in range(len(vector_BYland))]
    landReqd = [[] for i in range(len(vector_BYland))]
    bioProd= [[] for i in range(len(vector_BYland))]
    for i in range(len(vector_BYland)):
        for j in range(len(vector_BYland[i])/2):
            if not(vector_BYland[i][j*2] == 0):
                landArea[i] = [vector_BYland[i][j*2+1]] + landArea[i]
                bioPerYR[i] = [vector_BYland[i][j*2+1]*vector_BYland[i][j*2]] + bioPerYR[i]
        landReqd[i] = [sum(landArea[i][:k]) for k in range(0,len(landArea[i])+1)]
        bioProd[i] = [sum(bioPerYR[i][:k]) for k in range(0,len(bioPerYR[i])+1)]

    # Calculate linear estimators

    for i in range(len(landReqd)):
        if ((bioProd[i]) and (len(bioProd[i]) > 1) and not(bioProd[i][-1] == 0)):
            
            # Normalize the vectors
            y_s = [float(bioProd[i][j])/float(bioProd[i][-1]) for j in range(0,len(bioProd[i]))]
            x_s = [float(landReqd[i][j])/float(landReqd[i][-1]) for j in range(0,len(landReqd[i]))]
            
            #find bisection points (identify where the breakpoints area)
            x_lo_idx = [bisect.bisect(x_s,x_lo)-1,bisect.bisect(x_s,x_lo)]
            x_hi_idx = [bisect.bisect(x_s,x_hi)-1,bisect.bisect(x_s,x_hi)]
            
            m_interp_lo = (y_s[x_lo_idx[1]] - y_s[x_lo_idx[0]])/(x_s[x_lo_idx[1]] - x_s[x_lo_idx[0]])
            b_interp_lo = y_s[x_lo_idx[0]] - m_interp_lo*x_s[x_lo_idx[0]]
            
            y_lo = x_lo*m_interp_lo + b_interp_lo
            
            m_interp_hi = (y_s[x_hi_idx[1]] - y_s[x_hi_idx[0]])/(x_s[x_hi_idx[1]] - x_s[x_hi_idx[0]])
            b_interp_hi = y_s[x_hi_idx[0]] - m_interp_hi*x_s[x_hi_idx[0]]
            
            y_hi = x_hi*m_interp_hi + b_interp_hi
            
            # Get the slope and intersection of linear estimators
            [m_lo, b_lo] = linCalc(x_s[0],y_s[0],x_lo,y_lo,landReqd[i][-1],bioProd[i][-1])
            [m_med, b_med] = linCalc(x_lo,y_lo,x_hi,y_hi,landReqd[i][-1],bioProd[i][-1])
            [m_hi, b_hi] = linCalc(x_hi,y_hi,x_s[-1],y_s[-1],landReqd[i][-1],bioProd[i][-1])
            
            vector_BYland[i] = [m_lo,b_lo,m_med,b_med,m_hi,b_hi]
            
        else:
            vector_BYland[i] = [0,0,0,0,0,0]
        
    return [vector_BYland, vector_BYlandAvg]

# ------------------------Calculate the linear coefficients for the bioyield
def linCalc(x1,y1,x2,y2,xmax,ymax):
    if ((x2-x1) == 0) or ((y2-y1) == 0) or (xmax == 0) or (ymax == 0):
        return [0,0]
    
    m = ((y2-y1)*ymax)/((x2-x1)*xmax)
    b = ymax*(y1-((y2-y1)/(x2-x1))*x1)
    
    return [m,b]

# ------------------------Write bioyield results to GAMS-readable format
def writeBioyield(vectA,nameA,vectB,nameB,vectC,nameC,vectD,nameD,vectE,nameE,destFile):
    if (not(vectA) and not(vectB) and not(vectC) and not(vectD) and not(vectE)):
        return -1
        
    with open(destFile,"wb") as ofileDest:
        writerDest = csv.writer(ofileDest, delimiter=',', quoting=csv.QUOTE_NONNUMERIC) # writer is necessary to create csv
        headerRow = [0,0] + range(1,13) # top line that lists all cell ID
        writerDest.writerow(headerRow)
        
        for k in range(0,len(vectA)):
            row_k = [nameA, k+1] + [vectA[k]]*4 + [0]*8
            writerDest.writerow(row_k)
            
        for k in range(0,len(vectB)):
            row_k = [nameB, k+1] + [vectB[k]]*4 + [0]*8
            writerDest.writerow(row_k)
            
        for k in range(0,len(vectC)):
            row_k = [nameC, k+1] + [vectC[k]]*4 + [0]*8
            writerDest.writerow(row_k)
            
        for k in range(0,len(vectD)):
            row_k = [nameD, k+1] + [vectD[k]]*4 + [0]*8
            writerDest.writerow(row_k)
            
        for k in range(0,len(vectE)):
            row_k = [nameE, k+1] + [vectE[k]]*4 + [0]*8
            writerDest.writerow(row_k)
        
    iface.mainWindow().statusBar().showMessage("Writing Bioyield to CSV - Complete") # Update the user in the status bar
    
    return 1

# ------------------------Begin the work - call functions and subroutines


# reproject the layers
reprojSHP = layerReproj(countrySHP,fullRast_noConstr)
QgsMapLayerRegistry.instance().addMapLayer(reprojSHP)

if constRast:
    reproj_constRast = rastReproj(constRast,fullRast_noConstr,0)
    
if byARast:
    reproj_byARast = rastReproj(byARast,fullRast_noConstr,1)
if byBRast:
    reproj_byBRast = rastReproj(byBRast,fullRast_noConstr,1)
if byCRast:
    reproj_byCRast = rastReproj(byCRast,fullRast_noConstr,1)
if byDRast:
    reproj_byDRast = rastReproj(byDRast,fullRast_noConstr,1)
if byERast:
    reproj_byERast = rastReproj(byERast,fullRast_noConstr,1)

t_1 = time.time()
t_01 = t_1 - t_0
print "Time to reproject: %.2f s" % t_01 



# Apply the constraint map to the full rast
if constRast:
    fRxmin = fullRast_noConstr.extent().xMinimum()
    fRxmax = fullRast_noConstr.extent().xMaximum()
    fRymin = fullRast_noConstr.extent().yMinimum()
    fRymax = fullRast_noConstr.extent().yMaximum()
    extent = "%f,%f,%f,%f" %(fRxmin, fRxmax, fRymin, fRymax)
    
    eqn = 'A*B'
    output_tmpRast = processing.runalg("grass:r.mapcalculator",fullRast_noConstr.dataProvider().dataSourceUri(),reproj_constRast.dataProvider().dataSourceUri(),None,None,None,None,eqn,extent,0,None)
    fullRast = processing.getObject(output_tmpRast['outfile'])
    QgsMapLayerRegistry.instance().addMapLayer(fullRast)

t_2 = time.time()
t_12 = t_2 - t_1
print "Time to apply constraint to full raster: %.2f s" % t_12



# determine the necessary cell size to get approximately numGridCells of cells in the grid
reprojArea = (reprojSHP.extent().xMaximum()-reprojSHP.extent().xMinimum())*(reprojSHP.extent().yMaximum()-reprojSHP.extent().yMinimum())
cellDimension = math.sqrt(reprojArea/numGridCells)
gridSHP = gridGen(reprojSHP,cellDimension)

t_3 = time.time()
t_23 = t_3 - t_2
print "Time to create initial grid: %.2f s" % t_23



# replace the attribute fields with the ID and area
tightGridSHP = replaceFields(gridSHP)

t_4 = time.time()
t_34 = t_4 - t_3
print "Time to replace fields: %.2f s" % t_34



# get the right number of grid cells requested, excluding ones that are too small
numFeat = len(list(tightGridSHP.getFeatures())) 
prevNum = numGridCells
smallFeatList = countSmallFeatures(tightGridSHP)
numSmallFeat = len(list(smallFeatList))

while ((iter < maxIter) and (((numFeat-numSmallFeat) < numGridCells*minUnder) or ((numFeat-numSmallFeat) > numGridCells*maxOver))):
    iter = iter+1
    newNum = prevNum + ((numGridCells - numFeat)+numSmallFeat)*.9
    cellDimension = math.sqrt(reprojArea/newNum)
    gridSHP = gridGen(reprojSHP,cellDimension)
    tightGridSHP = replaceFields(gridSHP)
    prevNum = newNum
    numFeat = len(list(tightGridSHP.getFeatures()))
    smallFeatList = countSmallFeatures(tightGridSHP)
    numSmallFeat = len(list(smallFeatList))
    print "*** Iteration number: %d" % iter
    
t_5 = time.time()
t_45 = t_5 - t_4
print "Time to create good grid: %.2f s" % t_45



# Remove grid cells that are too small from the layer
with edit(tightGridSHP):
    deleteFeat = countSmallFeatures(tightGridSHP)
    for feature in deleteFeat:
        tightGridSHP.deleteFeature(feature.id())

fullGridSHP = replaceFields(tightGridSHP)

print "*** Remaining grid cells: %d" % len(list(fullGridSHP.getFeatures()))

QgsMapLayerRegistry.instance().addMapLayer(fullGridSHP)
QgsVectorFileWriter.writeAsVectorFormat(fullGridSHP,outSHP,"utf-8",fullGridSHP.crs(),"ESRI Shapefile",False)

t_6 = time.time()
t_56 = t_6 - t_5
print "Time to remove small grid cells: %.2f s" % t_56



# Export Attribute Table csv
writeAttributeTableFlag = writeAttributeTable(fullGridSHP)

t_7 = time.time()
t_67 = t_7 - t_6
print "Time to export attribute table: %.2f s" % t_67



# Export Linear Distance matrix csv
writeLinDistFlag = writeLinDist(fullGridSHP)

t_8 = time.time()
t_78 = t_8 - t_7
print "Time to export centroid distances: %.2f s" % t_78



# Export Tortuosity matrix csv
writeToruosityFlag = writeTortuosity(fullGridSHP)
    

t_9 = time.time()
t_89 = t_9 - t_8
print "Time to calculate and export tortuosity: %.2f s" % t_89



# Perform the raster cell counts for each feature
gridRast = gridToRast(fullGridSHP)
QgsMapLayerRegistry.instance().addMapLayer(gridRast)

t_10 = time.time()
t_910 = t_10 - t_9
print "Time to rasterize grid: %.2f s" % t_910



tmp_availLand = gridStats(gridRast,fullRast,fullGridSHP)
availLand = [x/(1000.0*1000.0) for x in tmp_availLand]

t_11 = time.time()
t_1011 = t_11 - t_10
print "Time to stats: %.2f s" % t_1011



# save raster stats data to csv
flagLA = writeVec(availLand,outLAcsv)

t_12 = time.time()
t_1112 = t_12 - t_11
print "Time to write stats to csv: %.2f s" % t_1112



# Calculate the bio-yield
proc_byARast = []
proc_byBRast = []
proc_byCRast = []
proc_byDRast = []
proc_byERast = []
proc_byARastAvg = []
proc_byBRastAvg = []
proc_byCRastAvg = []
proc_byDRastAvg = []
proc_byERastAvg = []

if constRast:
    if byARast:
        [proc_byARast, proc_byARastAvg] = bioyield(reproj_constRast,reproj_byARast,fullRast_noConstr,fullGridSHP,gridRast)
    if byBRast:
        [proc_byBRast, proc_byBRastAvg] = bioyield(reproj_constRast,reproj_byBRast,fullRast_noConstr,fullGridSHP,gridRast)
    if byCRast:
        [proc_byCRast, proc_byCRastAvg] = bioyield(reproj_constRast,reproj_byCRast,fullRast_noConstr,fullGridSHP,gridRast)
    if byDRast:
        [proc_byDRast, proc_byDRastAvg] = bioyield(reproj_constRast,reproj_byDRast,fullRast_noConstr,fullGridSHP,gridRast)
    if byERast:
        [proc_byERast, proc_byERastAvg] = bioyield(reproj_constRast,reproj_byERast,fullRast_noConstr,fullGridSHP,gridRast)

t_13 = time.time()
t_1213 = t_13 - t_12
print "Time to do bioyield raster calculations: %.2f s" % t_1213


# do this right for averages
flagBYavg = writeBioyield(proc_byARastAvg, Bioyield_Product_A,proc_byBRastAvg, Bioyield_Product_B,proc_byCRastAvg, Bioyield_Product_C,proc_byDRastAvg, Bioyield_Product_D,proc_byERastAvg, Bioyield_Product_E,outBYcsv)

# Everything down - only if there is lin estimates - what about these output otherwise?
if (proc_byARast):
    ALUm = [i[0] for i in proc_byARast]
    AMUm = [i[2] for i in proc_byARast]
    AHUm = [i[4] for i in proc_byARast]
    ALUb = [i[1] for i in proc_byARast]
    AMUb = [i[3] for i in proc_byARast]
    AHUb = [i[5] for i in proc_byARast]
else:
    ALUm = []
    AMUm = []
    AHUm = []
    ALUb = []
    AMUb = []
    AHUb = []
if (proc_byBRast):
    BLUm = [i[0] for i in proc_byBRast]
    BMUm = [i[2] for i in proc_byBRast]
    BHUm = [i[4] for i in proc_byBRast]
    BLUb = [i[1] for i in proc_byBRast]
    BMUb = [i[3] for i in proc_byBRast]
    BHUb = [i[5] for i in proc_byBRast]
else:
    BLUm = []
    BMUm = []
    BHUm = []
    BLUb = []
    BMUb = []
    BHUb = []
if (proc_byCRast):
    CLUm = [i[0] for i in proc_byCRast]
    CMUm = [i[2] for i in proc_byCRast]
    CHUm = [i[4] for i in proc_byCRast]
    CLUb = [i[1] for i in proc_byCRast]
    CMUb = [i[3] for i in proc_byCRast]
    CHUb = [i[5] for i in proc_byCRast]
else:
    CLUm = []
    CMUm = []
    CHUm = []
    CLUb = []
    CMUb = []
    CHUb = []
if (proc_byDRast):
    DLUm = [i[0] for i in proc_byDRast]
    DMUm = [i[2] for i in proc_byDRast]
    DHUm = [i[4] for i in proc_byDRast]
    DLUb = [i[1] for i in proc_byDRast]
    DMUb = [i[3] for i in proc_byDRast]
    DHUb = [i[5] for i in proc_byDRast]
else:
    DLUm = []
    DMUm = []
    DHUm = []
    DLUb = []
    DMUb = []
    DHUb = []
if (proc_byERast):
    ELUm = [i[0] for i in proc_byERast]
    EMUm = [i[2] for i in proc_byERast]
    EHUm = [i[4] for i in proc_byERast]
    ELUb = [i[1] for i in proc_byERast]
    EMUb = [i[3] for i in proc_byERast]
    EHUb = [i[5] for i in proc_byERast]
else:
    ELUm = []
    EMUm = []
    EHUm = []
    ELUb = []
    EMUb = []
    EHUb = []

flagLUm = writeBioyield(ALUm, Bioyield_Product_A,BLUm, Bioyield_Product_B,CLUm, Bioyield_Product_C,DLUm, Bioyield_Product_D,ELUm, Bioyield_Product_E,outBYLUmcsv)
flagMUm = writeBioyield(AMUm, Bioyield_Product_A,BMUm, Bioyield_Product_B,CMUm, Bioyield_Product_C,DMUm, Bioyield_Product_D,EMUm, Bioyield_Product_E,outBYMUmcsv)
flagHUm = writeBioyield(AHUm, Bioyield_Product_A,BHUm, Bioyield_Product_B,CHUm, Bioyield_Product_C,DHUm,Bioyield_Product_D,EHUm, Bioyield_Product_E,outBYHUmcsv)
flagLUb = writeBioyield(ALUb, Bioyield_Product_A,BLUb, Bioyield_Product_B,CLUb, Bioyield_Product_C,DLUb, Bioyield_Product_D,ELUb, Bioyield_Product_E,outBYLUbcsv)
flagMUb = writeBioyield(AMUb, Bioyield_Product_A,BMUb, Bioyield_Product_B,CMUb, Bioyield_Product_C,DMUb, Bioyield_Product_D,EMUb, Bioyield_Product_E,outBYMUbcsv)
flagHUb = writeBioyield(AHUb, Bioyield_Product_A,BHUb, Bioyield_Product_B,CHUb, Bioyield_Product_C,DHUb, Bioyield_Product_D,EHUb, Bioyield_Product_E,outBYHUbcsv)

t_14 = time.time()
t_1314 = t_14 - t_13
print "Time to write bioyield results: %.2f s" % t_1314


QgsMapLayerRegistry.instance().removeMapLayer(gridRast)
if not(countrySHP.crs().authid() == fullRast_noConstr.crs().authid()):
    QgsMapLayerRegistry.instance().removeMapLayer(reprojSHP)
QgsMapLayerRegistry.instance().removeMapLayer(fullGridSHP)
if constRast:
    QgsMapLayerRegistry.instance().removeMapLayer(fullRast)

t_total = t_13 - t_0
print "Total time:"
print t_total
