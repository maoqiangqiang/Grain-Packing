#!/usr/bin/env python
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
##################################################################################################
import numpy as np
import matplotlib.pyplot as plt
import math
import random
from scipy.stats import truncnorm
import time
from awcdf import*


##############################################################Function##############################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Creat one Region~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#was not used actually
def CreateMatrix(row, line):
    # initial matrix
    fMatrix = np.zeros((row, line))
    # first column and last columnn
    for ii in np.arange(0, row):
        fMatrix[ii, 0] = 1
        fMatrix[ii, line-1] = 1
    # middle
    for i in np.arange(1, row):
        for j in np.arange(1, line-2):
            fMatrix[i, j] = 0
    # bottom row
    for jj in np.arange(1, line-1):
        fMatrix[0, jj] = 1
    return fMatrix



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def ballLocal(x, y, radius,row,line,visual):
    for i in np.arange(math.floor(x-radius),math.ceil(x+radius)):
        for j in np.arange(math.floor(y-radius),math.ceil(y+radius)):
            #if i < line:
            distance = math.sqrt((i-x)**2+(j-y)**2)
            if distance <= radius:
                visual[int(j),int(i%line)]=1
    return visual
'''
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~generate one ball within region~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def ball(x, y, radius):
    global row, line
    
    for i in np.arange(0, row-1):
        for j in np.arange(0, line-1):
            distance = math.sqrt(((j-x))**2+((i-y))**2)
            distanceNegative = math.sqrt(((j-(x-line))**2+(i-y)**2))                #To avoid generate ball beyond boundry
            distancePositive = math.sqrt(((j-(x+line))**2+(i-y)**2))                
            if distance <= radius or  distanceNegative <= radius or distancePositive <= radius:
                region[i, j] = 1      
    return region
'''

#~~~~~~~~~~~~~~~~~~~~~~~create another ball with random x , 300 and random radius~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def createBall(frow,fline,fminRad,fmaxRad,fgsdFunction,fprobability,fx,fmean,fstd,finterval,fa,fb):
    # print('\ncreating next ball now !!!')
    fGrainSize = {
        2: awGrainSize(fprobability,fx,fmean,fstd,finterval,fa,fb),
        1: truncnorm.ppf(random.random(),fa,fb,fmean,fstd),
        0: random.uniform(fminRad,fmaxRad)
    }[fgsdFunction]
    return random.uniform(0, fline),frow+fGrainSize,fGrainSize#random.uniform(fminRad,fmaxRad)

#~~~~~Traverse each ball and find the first ball (KissBall) as well as store each hitted ball in (ContactList)~~~~~~~~~~~~~~ 
def TraverseBall(newball, ballList,kissball,line):                                #here newball refers NextBall; ballList means BallList
    fContactList = []                                        #ContactList: store all the ball that would be hitted 
    kiss = 0
    if kissball:
        kissball=(kissball[0]%line,kissball[1],kissball[2])
        ballList.remove(kissball)

    for existedBall in ballList:
        if existedBall[1] >= newball[1]:            #ignore oldball which has same y and in a straight line
            continue
        
        else:
            distanceXNE1 = (newball[0]-existedBall[0])%line             #(a-b)%500
            distanceXNE2 = (existedBall[0]-newball[0])%line             #(b-a)%500
            MinDistanceXNE = min(distanceXNE1,distanceXNE2)             #FIND minimum distanc    
            #print(distanceXNE1,distanceXNE2,MinDistanceXNE)
            radiaSumNE = newball[2]+existedBall[2]

            if round(MinDistanceXNE) < round(radiaSumNE):
                kiss = 1
                fContactList.append(existedBall)
    
    #find highest y which be regarded as KissBall   / find maxium element in a tuple list 
    if fContactList:
        fKissBall= max(fContactList, key=lambda x:(x[1]+x[2])) 
        fKissBall = (round(fKissBall[0],10),round(fKissBall[1],10),round(fKissBall[2],10)) 
    else:
        fKissBall = None    
    #print('Traverse ball')
    return fKissBall, fContactList, kiss                   #return KissBall, ContactList and interact flag

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~determination of Stacking~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Calulate minimun distance of two points        ---------used for period boundry
def PeriodNegMinDistance(ball1,ball2,fline):                              #part of ball is on the left 
    #kissball is ball1, oldball is ball 2
    x1,y1,r1 = ball1[0],ball1[1],ball1[2]
    x2,y2,r2 = ball2[0],ball2[1],ball2[2]

    distanceNegative = math.sqrt(((x1-fline)-x2)**2+(y1-y2)**2)          
    return distanceNegative

def PeriodPosiMinDistance(ball1,ball2,fline):
    #kissball is ball1, oldball is ball 2
    x1,y1,r1 = ball1[0],ball1[1],ball1[2]
    x2,y2,r2 = ball2[0],ball2[1],ball2[2]

    distancePositive = math.sqrt(((x1+fline)-x2)**2+(y1-y2)**2)           
    return distancePositive
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#calculate 3 balls system and  get the top balls' position 
def insec(fKissBall,fStackBall,fNextBall):
   
    x,y,R = fKissBall[0],fKissBall[1],fKissBall[2]+fNextBall[2]
    a,b,S = fStackBall[0],fStackBall[1],fStackBall[2]+fNextBall[2]
   
    d = math.sqrt((abs(a-x))**2 + (abs(b-y))**2)

    A = (R**2 - S**2 + d**2) / (2 * d)
    h = math.sqrt(R**2 - A**2)
    x2 = x + A * (a-x)/d
    y2 = y + A * (b-y)/d
    x3 = (x2 - h * (b - y) / d)
    y3 = (y2 + h * (a - x) / d)
    x4 = (x2 + h * (b - y) / d)
    y4 = (y2 - h * (a - x) / d)
  
    insecPoint=[(x3,y3),(x4,y4)]
    final=max(insecPoint,key=lambda x:x[1])
    return final[0],final[1]

#find the highest y of StackBall:
def CloseBallIntersection(fKissBall,fCloseBall,fNextBall,fStackBallList,fline):
    ##here StackBallList means final position of nextball after stack
    #fStackballhighesr meas final position
    
    for oldball in fCloseBall:
        #print(oldball,fKissBall)
        fx,fy = insec(fKissBall,oldball,fNextBall)
        finalstackball = ((fx%fline),fy,fNextBall[2])
        fStackBallList.append(finalstackball)

    fStackBallHighest = max(fStackBallList,key=lambda x:x[1])
    for oldball in fCloseBall:
        fx,fy = insec(fKissBall,oldball,fNextBall)
        finalstackball = ((fx%fline),fy,fNextBall[2])
        if finalstackball == fStackBallHighest:
            fStackBall = oldball    
    return fStackBallHighest,fStackBall
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~end~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#chech and avoid overlapping when set nextball[1] equal to KissBall[1]
def CheckOverlapping(nextball,CorrectList,KissBall,fStackBall,fline):
    # print('now come to check overlapping')
    # print('see here stackball '+str(fStackBall))
    #print('fCheckBallList:'+str(fcheckBallList))
    unfinish =0
    disAll=[]
    OverlapBallList=[]
    OverlapBall=()
    fStackBall=(round(fStackBall[0]%fline,10),fStackBall[1],fStackBall[2])
    CorrectList.remove(KissBall)
    CorrectList.remove(fStackBall)
    #print(CorrectList)
    for oldball in CorrectList:
        #if round(oldball[1],10)>=round(nextball[1],10):           ##==failed
        #    continue
        #else:   
            #print(oldball)
        dis1=math.sqrt((nextball[0]-oldball[0])**2+(nextball[1]-oldball[1])**2)
        dis2=math.sqrt((nextball[0]+fline-oldball[0])**2+(nextball[1]-oldball[1])**2)
        dis3=math.sqrt((nextball[0]-fline-oldball[0])**2+(nextball[1]-oldball[1])**2)
        disAll=[dis1,dis2,dis3]
        dis = min(disAll)
        if round(dis) < round(nextball[2]+oldball[2]):
            # print(round(dis),round(nextball[2]+oldball[2]))
            
            # print(oldball)
            # print('this Stack is useless')
            unfinish = 1
            OverlapBallList.append(oldball)
        if OverlapBallList:
            OverlapBall = max(OverlapBallList, key=lambda x:(x[1]+x[2]))
            # print('see here'+str(OverlapBall))
    return unfinish,OverlapBall


















#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Ending~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
