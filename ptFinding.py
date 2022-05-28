import time

from Data_Processing import *
from Base_Functions import *
import json

import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate

class SparLocation:
    def __init__(self, xc, zc, t):
        self.Xc = xc
        self.Zc = zc
        self.theta = t

#generates a filename string for any file needed for data management
#designed to make file structure more managable and easier to change
def gen_filename(USER, des_file, save):

    if USER == 'SHARED':
        user_dir = 'shared_data\\'
    else:
        user_dir = f'user_data\{USER}\\'

    if des_file == 'user_preferences':
        type_dir = ''
    elif des_file == 'active_file_list':
        type_dir = ''
    elif 'database' in des_file:
        type_dir = 'databases\\'
    elif 'log' in des_file:
        type_dir = 'logs\\'
    elif 'Geometry' in des_file:
        type_dir = f'files\{save}\\temp_save_data\geometry_files\\'
    elif 'Geo' in des_file:
        type_dir = f'files\{save}\\temp_save_data\\'
    elif 'Properties' in des_file:
        type_dir = f'files\{save}\\temp_save_data\\'
    elif des_file == 'active_body_list':
        type_dir = f'files\{save}\\temp_save_data\\'

    if '.jpg' in des_file:
        file_type = ''
    if '.json' in des_file:
        file_type = ''
    else:
        file_type = '.csv'

    filename = f'{user_dir}{type_dir}{des_file}{file_type}'

    return filename

def get_airfoil_database2():

    shared_airfoil_database_path = gen_filename('SHARED', 'airfoil_database.json', 0)

    with open(shared_airfoil_database_path) as json_file:
        data = json.load(json_file)

    airfoil_database_strings = data['name_list']
    airfoil_database_x = data['x_list']
    airfoil_database_y = data['y_list']

    return [airfoil_database_strings, airfoil_database_x, airfoil_database_y]

#seperate the upper and lower profiles of the airfoil and providing lists of points for plotting and calculations
def seperateUpperLower(airfoilptsX, airfoilptsY):
    # seperate upper and lower sections of airfoil
    upper_x = []
    upper_y = []
    lower_x = [0]  # give initial origin point
    lower_y = [0]  # give initial origin point
    on_up = 1
    for i in range(len(airfoilptsX)):

        x = airfoilptsX[i]
        y = airfoilptsY[i]

        if on_up == 1:
            upper_x.append(x)
            upper_y.append(y)

        elif on_up == 0:
            lower_x.append(x)
            lower_y.append(y)

        if i > 0 and x == 0 and y == 0:
            on_up = 0

    return [upper_x, upper_y, lower_x, lower_y]

#from Paul Draper (https://stackoverflow.com/questions/20677795/how-do-i-compute-the-intersection-point-of-two-lines)
def line_intersection(line1, line2):
    xdiff = (line1[0][0] - line1[1][0], line2[0][0] - line2[1][0])
    ydiff = (line1[0][1] - line1[1][1], line2[0][1] - line2[1][1])

    def det(a, b):
        return a[0] * b[1] - a[1] * b[0]

    div = det(xdiff, ydiff)
    if div == 0:
       raise Exception('lines do not intersect')

    d = (det(*line1), det(*line2))
    x = det(d, xdiff) / div
    y = det(d, ydiff) / div
    return [x, y]

#spar theta must be in radians
def findSparIntersection(sparXc, sparZc, sparTheta, splineX, splineZ):

    #find upper index closest to intersection

    #find the delta in X and Z axes of all airfoil pts from the spar center
    Xdiff = [x - sparXc for x in splineX]
    Zdiff = [z - sparZc for z in splineZ]

    theta = []

    #find the angles of all the airfoil pts from the spar center pt
    for i in range(len(Zdiff)):

        thetaN = np.arctan2(Zdiff[i], Xdiff[i])

        theta.append(thetaN)

    #calculate the difference in angle to each point from the spar angle
    thetaDiff = [t - sparTheta for t in theta]

    #find the 2 airfoil pts that the spar angle passes between
    thetaIndexPos = thetaDiff.index(min([i for i in thetaDiff if i > 0]))
    thetaIndexNeg = thetaDiff.index(max([i for i in thetaDiff if i < 0]))

    #plot bounding spline points of intersection
    #plt.scatter(splineX[thetaIndexNeg], splineZ[thetaIndexNeg], s=2)
    #plt.scatter(splineX[thetaIndexPos], splineZ[thetaIndexPos], s=2)

    #get the definition of the line that the spar intersects
    boundingLine = [[splineX[thetaIndexNeg], splineZ[thetaIndexNeg]], [splineX[thetaIndexPos], splineZ[thetaIndexPos]]]

    #output the bounding line and the indeces of surrounding pts
    return boundingLine, thetaIndexNeg, thetaIndexPos


def findSparIntersections(sparXc, sparZc, sparTheta):
    sparThetaOpp = sparTheta - np.pi  # the opposite angle of the spar to be used for the lower section of the airfoil

    # find points to draw large spar line
    legLen = 0.2
    sZmax = np.sin(sparTheta) * legLen + sparZc
    sXmax = np.cos(sparTheta) * legLen + sparXc
    sZmin = -np.sin(sparTheta) * legLen + sparZc
    sXmin = -np.cos(sparTheta) * legLen + sparXc

    [upperX, upperZ, lowerX, lowerZ] = seperateUpperLower(airfoilX, airfoilZ)

    # find intersection lines of spar on upper and lower profiles
    upperLine, UnegIndex, UposIndex = findSparIntersection(sparXc, sparZc, sparTheta, upperX, upperZ)
    lowerLine, LnegIndex, LposIndex = findSparIntersection(sparXc, sparZc, sparThetaOpp, lowerX, lowerZ)

    # find intersection points on upper and lower profiles
    upperIntersectionX, upperIntersectionZ = line_intersection(upperLine, [[sXmin, sZmin], [sXmax, sZmax]])
    lowerIntersectionX, lowerIntersectionZ = line_intersection(lowerLine, [[sXmin, sZmin], [sXmax, sZmax]])


    #plt.scatter(upperIntersectionX, upperIntersectionZ, s=3)
    #plt.scatter(lowerIntersectionX, lowerIntersectionZ, s=3)

    plt.scatter(sparXc, sparZc, s=3)
    plt.plot([sXmin, sXmax], [sZmin, sZmax], linewidth=0.25, c="blue")

    return [[upperIntersectionX, upperIntersectionZ], [lowerIntersectionX, lowerIntersectionZ]], [[UnegIndex, UposIndex], [LnegIndex, LposIndex]]

#calculate the points and plot wingbox based on spar locations
def calculateWingbox(Spar1, Spar2, shapeInd, T1dir, T2dir):
    sparXc_1 = Spar1.Xc
    sparZc_1 = Spar1.Zc
    sparTheta_1 = Spar1.theta

    sparXc_2 = Spar2.Xc
    sparZc_2 = Spar2.Zc
    sparTheta_2 = Spar2.theta

    upShapeList = [1,1,1,1,0,0,0,1,1,1,0]
    lowShapeList = [1,1,1,0,1,1,1,0,0,0,1]
    frontShapeList = [1,1,0,1,1,1,0,1,0,0,0]
    backShapeList = [1,0,1,1,1,0,1,0,1,0,0]

    T1ptList = [-1,0,2,3,1,1,2,0,3,0,2]
    T2ptList = [-1,3,1,2,0,3,0,2,1,1,3]

    T1dirList = [-1, 1, 1, T1dir, T1dir, T1dir, 1, 1, T1dir, 1, 1]
    T2dirList = [-1, 0, 0, T2dir, T2dir, 0, T2dir, T2dir, 0, 0, 0]


    # define what to draw
    drawUpperSpan = upShapeList[shapeInd]
    drawLowerSpan = lowShapeList[shapeInd]
    drawInitialSpar = frontShapeList[shapeInd]
    drawSecondSpar = backShapeList[shapeInd]

    #find intersection of 1st and second spars
    #also pull in the indeces of surrounding pts to draw spans and tabs
    [[UInterX_1, UInterZ_1],[LInterX_1, LInterZ_1]],[[UnegIndex_1, UposIndex_1], [LnegIndex_1, LposIndex_1]] = findSparIntersections(sparXc_1, sparZc_1, sparTheta_1)
    [[UInterX_2, UInterZ_2],[LInterX_2, LInterZ_2]],[[UnegIndex_2, UposIndex_2], [LnegIndex_2, LposIndex_2]] = findSparIntersections(sparXc_2, sparZc_2, sparTheta_2)

    #organize pts to find the list of points between intersection of upper and lower surfaces
    upperSpanX = upperX[UposIndex_2:UnegIndex_1+1]
    upperSpanZ = upperZ[UposIndex_2:UnegIndex_1+1]

    upperSpanX = [UInterX_2] + upperSpanX + [UInterX_1]
    upperSpanZ = [UInterZ_2] + upperSpanZ + [UInterZ_1]

    lowerSpanX = lowerX[LposIndex_1:LnegIndex_2+1]
    lowerSpanZ = lowerZ[LposIndex_1:LnegIndex_2+1]

    lowerSpanX = [LInterX_1] + lowerSpanX + [LInterX_2]
    lowerSpanZ = [LInterZ_1] + lowerSpanZ + [LInterZ_2]

    frontSparX = [UInterX_1, LInterX_1]
    frontSparZ = [UInterZ_1, LInterZ_1]
    rearSparX = [UInterX_2, LInterX_2]
    rearSparZ = [UInterZ_2, LInterZ_2]

    # Plot Sections
    if drawInitialSpar == 1:
        plt.plot([UInterX_1, LInterX_1], [UInterZ_1, LInterZ_1], c="blue")  # plot initial Spar
    if drawSecondSpar == 1:
        plt.plot([UInterX_2, LInterX_2], [UInterZ_2, LInterZ_2], c="blue")  # plot initial Spar
    if drawUpperSpan == 1:
        plt.plot(upperSpanX, upperSpanZ, c="blue")
    if drawLowerSpan == 1:
        plt.plot(lowerSpanX, lowerSpanZ, c="blue")

    # ToDo: figure out where to place tabs based on what type of wingbox is being used (also apply to LE and TE)

    #calculuate and add tabs
    TabXptOrder = [UInterX_2,UInterX_1,LInterX_1,LInterX_2]
    TabZptOrder = [UInterZ_2,UInterZ_1,LInterZ_1,LInterZ_2]

    TabNegIndOrder = [UnegIndex_2,UnegIndex_1,LnegIndex_1,LnegIndex_2]
    TabPosIndOrder = [UposIndex_2, UposIndex_1, LposIndex_1, LposIndex_2]

    Tab1SurfList = [-1,0,1,1,0,0,1,0,1,0,1]
    Tab2SurfList = [-1,1,0,1,0,1,0,1,0,0,1]

    TabXSurfList = [upperX, lowerX]
    TabZSurfList = [upperZ, lowerZ]

    T1 = T1ptList[shapeInd]
    T2 = T2ptList[shapeInd]

    T1Surf = Tab1SurfList[shapeInd]
    T2Surf = Tab2SurfList[shapeInd]

    if shapeInd != 0:
        Tab_1X, Tab_1Z = calculateTab([TabXptOrder[T1], TabZptOrder[T1]], [TabNegIndOrder[T1], TabPosIndOrder[T1]], T1dirList[shapeInd], 0.05, TabXSurfList[T1Surf], TabZSurfList[T1Surf])
        Tab_2X, Tab_2Z = calculateTab([TabXptOrder[T2], TabZptOrder[T2]], [TabNegIndOrder[T2], TabPosIndOrder[T2]], T2dirList[shapeInd], 0.05, TabXSurfList[T2Surf], TabZSurfList[T2Surf])

        plt.plot(Tab_1X, Tab_1Z, 'red')
        plt.plot(Tab_2X, Tab_2Z, 'red')

    else:
        Tab_1X = []
        Tab_1Z = []
        Tab_2X = []
        Tab_2Z = []

    memberStartList = [0,0,2,3,1,1,2,0,3,0,2]

    if memberStartList[shapeInd] == 0:
        surfListX = upperSpanX * (drawUpperSpan) + frontSparX * (drawInitialSpar) + lowerSpanX * (drawLowerSpan) + rearSparX * (
        drawSecondSpar)

        surfListZ = upperSpanZ * (drawUpperSpan) + frontSparZ * (drawInitialSpar) + lowerSpanZ * (
            drawLowerSpan) + rearSparZ * (drawSecondSpar)
    elif memberStartList[shapeInd] == 1:
        surfListX = frontSparX * (drawInitialSpar) + lowerSpanX * (drawLowerSpan) + rearSparX * (
        drawSecondSpar) + upperSpanX * (drawUpperSpan)

        surfListZ = frontSparZ * (drawInitialSpar) + lowerSpanZ * (drawLowerSpan) + rearSparZ * (
            drawSecondSpar) + upperSpanZ * (drawUpperSpan)

    elif memberStartList[shapeInd] == 2:
        surfListX = lowerSpanX * (drawLowerSpan) + rearSparX * (
        drawSecondSpar) + upperSpanX * (drawUpperSpan) + frontSparX * (drawInitialSpar)

        surfListZ = lowerSpanZ * (drawLowerSpan) + rearSparZ * (
            drawSecondSpar) + upperSpanZ * (drawUpperSpan) + frontSparZ * (drawInitialSpar)

    elif memberStartList[shapeInd] == 3:
        surfListX = rearSparX * (
            drawSecondSpar) + upperSpanX * (drawUpperSpan) + frontSparX * (drawInitialSpar) + lowerSpanX * (drawLowerSpan)

        surfListZ = rearSparZ * (
            drawSecondSpar) + upperSpanZ * (drawUpperSpan) + frontSparZ * (drawInitialSpar) + lowerSpanZ * (drawLowerSpan)

    Tab_1X.reverse()
    Tab_1Z.reverse()

    #create a total list of points for the member being created
    memberX = Tab_1X + surfListX + Tab_2X
    memberZ = Tab_1Z + surfListZ + Tab_2Z

    #plot of the full member in order
    #plt.plot(memberX,memberZ)

    calculateCentroid(memberX, memberZ)

#calculate the points and plot wingbox based on spar locations
def calculateSingleSpar(Spar1, sparShape):
    sparXc_1 = Spar1.Xc
    sparZc_1 = Spar1.Zc
    sparTheta_1 = Spar1.theta

    sparShapes = ['C','invC','Z','invZ']

    sparShapeInd = sparShapes.index(sparShape)

    T1dirList = [1,0,0,1]
    T2dirList = [0,1,0,1]

    T1dir = T1dirList[sparShapeInd]
    T2dir = T2dirList[sparShapeInd]

    # define what to draw
    drawInitialSpar = 0

    #find intersection of 1st and second spars
    #also pull in the indeces of surrounding pts to draw spans and tabs
    [[UInterX_1, UInterZ_1],[LInterX_1, LInterZ_1]],[[UnegIndex_1, UposIndex_1], [LnegIndex_1, LposIndex_1]] = findSparIntersections(sparXc_1, sparZc_1, sparTheta_1)

    # Plot Sections
    if drawInitialSpar == 1:
        plt.plot([UInterX_1, LInterX_1], [UInterZ_1, LInterZ_1], c="blue")  # plot initial Spar

    #calculuate and add tabs for spar
    UTab_1X, UTab_1Z = calculateTab([UInterX_1, UInterZ_1], [UnegIndex_1, UposIndex_1], T1dir, 0.05, upperX, upperZ)
    LTab_1X, LTab_1Z = calculateTab([LInterX_1, LInterZ_1], [LnegIndex_1, LposIndex_1], T2dir, 0.05, lowerX, lowerZ)

    #link points for inital spar
    UTab_1X.reverse()
    UTab_1Z.reverse()
    SparPtsX_1 = UTab_1X + LTab_1X
    SparPtsZ_1 = UTab_1Z + LTab_1Z

    plt.plot(SparPtsX_1, SparPtsZ_1)
    calculateCentroid(SparPtsX_1, SparPtsZ_1)

#find the centroid of a 2d list of points with even line weights
def calculateCentroid(memberX, memberZ):

    xSum = 0
    zSum = 0

    lenSum = 0

    for n in range(len(memberX) - 1):

        x1 = memberX[n]
        x2 = memberX[n+1]

        z1 = memberZ[n]
        z2 = memberZ[n+1]

        xMid = (x1+x2)/2
        zMid = (z1+z2)/2

        segLen = np.sqrt((x2-x1)**2 + (z2-z1)**2)

        xSum += xMid * segLen
        zSum += zMid * segLen

        lenSum +=segLen

    xCentroid = xSum / lenSum
    zCentroid = zSum / lenSum

    plt.scatter(xCentroid, zCentroid)

#calculate the wingbox that goes around the leading edge
def calculateLeadingEdge(Spar):

    sparXc_1 = Spar.Xc
    sparZc_1 = Spar.Zc
    sparTheta_1 = Spar.theta

    #define what to draw
    drawLE = 0
    drawSpar = 0

    # find intersection of spar
    # also pull in the indeces of surrounding pts to draw spans and tabs
    [[UInterX_1, UInterZ_1], [LInterX_1, LInterZ_1]], [[UnegIndex_1, UposIndex_1],[LnegIndex_1, LposIndex_1]] = findSparIntersections(sparXc_1,
                                                                                                           sparZc_1,
                                                                                                           sparTheta_1)

    # organize pts to find the list of points between intersection of upper and lower surfaces
    upperSpanX = upperX[UposIndex_1:-1]
    upperSpanZ = upperZ[UposIndex_1:-1]

    upperSpanX = [UInterX_1] + upperSpanX
    upperSpanZ = [UInterZ_1] + upperSpanZ

    lowerSpanX = lowerX[0:LposIndex_1]
    lowerSpanZ = lowerZ[0:LposIndex_1]

    lowerSpanX = lowerSpanX + [LInterX_1]
    lowerSpanZ = lowerSpanZ + [LInterZ_1]

    LEspanX = upperSpanX + lowerSpanX
    LEspanZ = upperSpanZ + lowerSpanZ

    if drawSpar == 1:
        plt.plot([UInterX_1, LInterX_1], [UInterZ_1, LInterZ_1], c="red")  # plot initial Spar
    if drawLE == 1:
        plt.plot(LEspanX,LEspanZ, c="red")

#calculate the wingbox that goes around the trailing edge
def calculateTrailingEdge(Spar):

    sparXc_1 = Spar.Xc
    sparZc_1 = Spar.Zc
    sparTheta_1 = Spar.theta

    #define what to draw
    drawTE = 0
    drawSpar = 0

    # find intersection of spar
    # also pull in the indeces of surrounding pts to draw spans and tabs
    [[UInterX_1, UInterZ_1], [LInterX_1, LInterZ_1]], [[UnegIndex_1, UposIndex_1],[LnegIndex_1, LposIndex_1]] = findSparIntersections(sparXc_1,
                                                                                                           sparZc_1,
                                                                                                           sparTheta_1)

    # organize pts to find the list of points between intersection of upper and lower surfaces
    upperSpanX = upperX[0:UposIndex_1]
    upperSpanZ = upperZ[0:UposIndex_1]

    upperSpanX = upperSpanX + [UInterX_1]
    upperSpanZ = upperSpanZ + [UInterZ_1]

    lowerSpanX = lowerX[LposIndex_1:-1]
    lowerSpanZ = lowerZ[LposIndex_1:-1]

    lowerSpanX = [LInterX_1] + lowerSpanX
    lowerSpanZ = [LInterZ_1] + lowerSpanZ

    LEspanX = lowerSpanX + upperSpanX
    LEspanZ = lowerSpanZ + upperSpanZ

    if drawSpar == 1:
        plt.plot([UInterX_1, LInterX_1], [UInterZ_1, LInterZ_1], c="red")  # plot initial Spar
    if drawTE == 1:
        plt.plot(LEspanX,LEspanZ, c="red")

#TabDirection: 0 = neg, 1 = pos
def calculateTab(IntersectionPt, IntersectionBoundingIndeces, TabDirection, TabLength, profileX, profileZ):

    #for loop that increments the index in the Tab direction from the intersection pt, until the total length of all additions equals that of the TabLength


    if TabDirection == 1:
        initialIndex = IntersectionBoundingIndeces[TabDirection]
        initialTabX = profileX[0:initialIndex] + [IntersectionPt[0]] #add intersection point to tab list
        initialTabX.reverse() #reverse so that intersection point is the first point
        initialTabZ = profileZ[0:initialIndex] + [IntersectionPt[1]]
        initialTabZ.reverse()

    elif TabDirection == 0:
        initialIndex = IntersectionBoundingIndeces[TabDirection] + 1
        initialTabX = [IntersectionPt[0]] + profileX[initialIndex:]
        initialTabZ = [IntersectionPt[1]] + profileZ[initialIndex:]

    #chop tab to the proper length
    TabX = [initialTabX[0]]
    TabZ = [initialTabZ[0]]
    totalLen = 0
    lastX = initialTabX[0]
    lastZ = initialTabZ[0]
    for i in range(len(initialTabX)-1):

        currX = initialTabX[i+1]
        currZ = initialTabZ[i+1]

        displacement = np.sqrt((currX-lastX)**2 + (currZ-lastZ)**2)

        totalLen += displacement

        if totalLen <= TabLength:
            TabX.append(currX)
            TabZ.append(currZ)

            lastX = currX
            lastZ = currZ
        else:
            diffX = (currX-lastX)
            diffZ = (currZ-lastZ)

            diffPerc = (displacement - (totalLen - TabLength)) / displacement

            newX = lastX + (diffX * diffPerc)
            newZ = lastZ + (diffZ * diffPerc)

            TabX.append(newX)
            TabZ.append(newZ)

            break

    return TabX, TabZ


#find the angle opposite of l3
def lawOfCos(l1, l2, l3):

    angle = (np.arccos((l1**2 + l2**2 - l3**2) / (2*l1*l2)))

    return angle

#find the hinge profile based on an optimal hinge radius to be tangent to the upper surface of the wing
def findHingeIntersection(hingeX, hingeZ, upperSplineX, upperSplineZ, lowerSplineX, lowerSplineZ):

    # find the delta in X and Z axes of all airfoil pts from the spar center
    Xdiff = [x - hingeX for x in upperSplineX]
    Zdiff = [z - hingeZ for z in upperSplineZ]

    linearDiff=[]
    minusMinDiff=[]

    for i in range(len(Zdiff)):

        xdiff = Xdiff[i]
        zdiff = Zdiff[i]

        linearDiff.append(np.sqrt(xdiff**2 + zdiff**2))
        minusMinDiff.append(np.sqrt(xdiff ** 2 + zdiff ** 2))

    minDiffInd = linearDiff.index(min([n for n in linearDiff if n>0]))
    del minusMinDiff[minDiffInd]
    secminDiffInd = minusMinDiff.index(min([i for i in minusMinDiff if i>0]))


    plt.scatter(hingeX,hingeZ,s=3) #plot hinge pt

    l3Xdiff = upperSplineX[secminDiffInd] - upperSplineX[minDiffInd]
    l3Zdiff = upperSplineZ[secminDiffInd] - upperSplineZ[minDiffInd]

    diff1 = linearDiff[minDiffInd]
    diff2 = linearDiff[secminDiffInd]

    diff3 = np.sqrt(l3Xdiff**2 + l3Zdiff**2)

    alpha = lawOfCos(diff1,diff3,diff2)

    phi = (np.pi/2) - alpha

    percdiff = diff1*np.sin(phi)
    perc = percdiff / diff3

    tangentPtX = upperSplineX[minDiffInd] + (l3Xdiff * perc)
    tangentPtZ = upperSplineZ[minDiffInd] + (l3Zdiff * perc)

    plt.scatter(tangentPtX, tangentPtZ)

    hingeType = 0 #may be used in future to specify how to shape the leading edge of the control surface

    if hingeType == 0:
        hingeR = diff1 * np.cos(phi)

    # find the delta in X and Z axes of all airfoil pts from the spar center
    Xdiff = [x - hingeX for x in lowerSplineX]
    Zdiff = [z - hingeZ for z in lowerSplineZ]

    linearDiff = []
    minusMinDiff = []
    linearDiffpastHinge = []

    for i in range(len(Zdiff)):
        xdiff = Xdiff[i]
        zdiff = Zdiff[i]

        if lowerSplineX[i] < hingeX:

            linearDiff.append(np.sqrt(xdiff ** 2 + zdiff ** 2))
            minusMinDiff.append(np.sqrt(xdiff ** 2 + zdiff ** 2))

        linearDiffpastHinge.append(np.sqrt(xdiff ** 2 + zdiff ** 2))

    minDiffInd = linearDiff.index(min([n for n in linearDiff if n > 0]))

    if linearDiff[minDiffInd] > hingeR:
        secminDiffInd = linearDiffpastHinge.index(min([i for i in linearDiffpastHinge if i > 0]))
        l2 = linearDiffpastHinge[secminDiffInd]
    else:
        del minusMinDiff[minDiffInd]
        secminDiffInd = minusMinDiff.index(min([i for i in minusMinDiff if i > 0]))
        l2 = linearDiff[secminDiffInd]

    plt.scatter([lowerSplineX[minDiffInd], lowerSplineX[secminDiffInd]],
                [lowerSplineZ[minDiffInd], lowerSplineZ[secminDiffInd]])

    l1 = linearDiff[minDiffInd]

    l3Xdiff = lowerSplineX[secminDiffInd] - lowerSplineX[minDiffInd]
    l3Zdiff = lowerSplineZ[secminDiffInd] - lowerSplineZ[minDiffInd]

    l3 = np.sqrt(l3Xdiff**2 + l3Zdiff**2)

    phi = lawOfCos(l1,l2,l3)

    alpha = lawOfCos(l2,l3,l1)

    beta = lawOfCos(l1,l3,l2)

    omega = (l1 / hingeR) * np.sin(beta)

    theta = np.pi - beta - omega

    x = hingeR * (np.sin(theta) / np.sin(beta))

    perc = x / l3

    lowerIntPtX = lowerSplineX[minDiffInd] + (l3Xdiff * perc)
    lowerIntPtZ = lowerSplineZ[minDiffInd] + (l3Zdiff * perc)

    plt.scatter(lowerIntPtX, lowerIntPtZ)

    intDeltaX = (lowerIntPtX - hingeX)
    intDeltaZ = (lowerIntPtZ - hingeZ)

    if intDeltaZ < 0 and intDeltaX < 0:
        lowerIntAngle = np.arctan(intDeltaZ / intDeltaX)
        lowerIntAngle += np.pi
    else:
        lowerIntAngle = np.arccos((lowerIntPtX - hingeX) / hingeR)

    upperTangentAngle = np.arccos((tangentPtX - hingeX) / hingeR)



    circPtsX, circPtsZ = arc_pts(hingeX,hingeZ,hingeR, lowerIntAngle,upperTangentAngle,50)

    plt.plot(circPtsX, circPtsZ)

    #ToDo: create list of points that represents everythign infront of the hinge (including the hinge), and a seperate list of all points that would be on the control surface


######################

# ToDo: Figure out how to prevent error when initial point lies outside of airfoil boundaries

# User should define a list of spars
# User should then define wingboxes based on 1 or 2 spars, then define which parts of the box will actually be applied

t0 = time.time()

#define the airfoil to use from the database
airfoil = "NLF(1)-0215F" #"NACA2414"

#define the initial spar
sparXc_1 = .2
sparZc_1 = 0
sparTheta_1 = (1 / 2) * np.pi

Spar1 = SparLocation(sparXc_1,sparZc_1,sparTheta_1)


#define second spar
sparXc_2 = .7
sparZc_2 = 0
sparTheta_2 = (1 / 3) * np.pi

Spar2 = SparLocation(sparXc_2,sparZc_2,sparTheta_2)

[airfoil_database_strings, airfoil_database_x, airfoil_database_y] = get_airfoil_database2()

airfoil_index = airfoil_database_strings.index(airfoil)


hingeX = 0.8
hingeZ = 0

#pull airfoil pts list
airfoilX = airfoil_database_x[airfoil_index]
airfoilZ = airfoil_database_y[airfoil_index]

[upperX, upperZ, lowerX, lowerZ] = seperateUpperLower(airfoilX, airfoilZ)

#plot the base airfoil in the back
plt.plot(airfoilX, airfoilZ, linewidth=0.5, c='red')


findHingeIntersection(hingeX,hingeZ,upperX,upperZ,lowerX,lowerZ)

#calculate the center wingbox based on 2 spar locations
calculateWingbox(Spar1, Spar2, 5, 1, 0) #prob at 3, 8
#calculateSingleSpar(Spar1,'Z')

calculateLeadingEdge(Spar1)
calculateTrailingEdge(Spar2)

t1 = time.time()

print("Time to process: " + str(t1-t0))

plt.axis('equal')
plt.show()


