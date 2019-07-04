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
#import matplotlib.cm as cm
import math
import random
import time
from scipy.stats import truncnorm
from GrainPacking2D_Function import*
from awcdf import*

def GrainPacking(fgsdFunction,fmean,fstd,fminrad,fmaxrad):
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    gsdFunction=fgsdFunction #0 is uniform, 1 is normal distribution, 2 is area weighted normal distribution
    random.seed(1)
    """initial parameters"""
    row = 1300#750#1300
    line = 1000#600#1000                  #model legnth
    minRad = fminrad #10
    maxRad = fmaxrad #25

    mu=fmean#random.uniform(10,25)
    sigma=fstd#random.uniform(10,25)

    a1,a2=mu-2*sigma,mu+2*sigma
    x1,x2,interval= a1,a2,1000
    x= np.linspace(x1,x2,interval)
    a,b=(a1-mu)/sigma,(a2-mu)/sigma

    proList=[]
    probability = random.random()
    proList.append(probability)
    ###########################################################
    GrainSize = {
        2: awGrainSize(probability,x,mu,sigma,interval,a,b),
        1: truncnorm.ppf(random.random(),a,b,mu,sigma),
        0: random.uniform(minRad,maxRad)
    }[gsdFunction]

    # print('first Grain Size = '+str(GrainSize))
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Creat Region ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    radius = GrainSize
    ballFirst=(random.uniform(0,line),radius,radius)
    # print('first ball '+str(ballFirst))
    BallList=[]
    ballFirst=(round(ballFirst[0],10),round(ballFirst[1],10),round(ballFirst[2],10))
    BallList.append(ballFirst)

    ################################################### generate ball one by one and check for hit only #######################################################
    # BallNumbers=250
    # for ii in range(0,BallNumbers):
    genGrain = 0
    while genGrain == 0:
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~create a new one and traverse contact~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        probability = random.random()
        proList.append(probability)
        NextBall=createBall(row,line,minRad,maxRad,gsdFunction,probability,x,mu,sigma,interval,a,b)                                           #generate a new next ball
 
        # print(NextBall)
        finish =0
        KissBall=()
        rerotate =0
        # print(len(BallList))
        while finish != 1:
            NextBall = (round(NextBall[0],10),round(NextBall[1],10),round(NextBall[2],10))
            # print('nextball'+str(NextBall))
            if KissBall:
                KissBall = (round(KissBall[0],10),round(KissBall[1],10),round(KissBall[2],10))
                if rerotate !=1:
                    TraverseKissBall=KissBall
                    TraverseBallList=BallList[:]
                    KissBall, ContactList, interact = TraverseBall(NextBall,TraverseBallList,TraverseKissBall,line)

            else:
                TraverseKissBall=()
                TraverseBallList=BallList[:]
                KissBall, ContactList, interact = TraverseBall(NextBall,TraverseBallList,TraverseKissBall,line)
                
            # print ('kissball '+str(KissBall))

            CloseBall=[]
            CloseBallSpecial =[]
            StackBallList=[]
            CheckBallList = BallList[:]
            TBallList=BallList[:]
            OverlapBall=()
            if len(ContactList)==0:
                #no interaction
                NextBall=(NextBall[0]%line,NextBall[2],NextBall[2])
                finish = 1
                # print('no kissball and then generating this ball on ground\n')
                
            elif len(ContactList)>0:
                #will be hitted
                #Case 1 :this situation is that KissBall has part of leftside and NextBall hit on the Left Period Boundry
                if (KissBall[0]+KissBall[2]>line and NextBall[0]<=(KissBall[2]+NextBall[2])) or (NextBall[0]-NextBall[2]<0 and KissBall[0]>=(line-NextBall[2]-KissBall[2])):    
                    TBallList.remove(KissBall)
                    #print(TBallList)
                    for oldball in TBallList:
                        if oldball[0]<NextBall[0]:      
                            continue
                        else:
                            dis1 = PeriodNegMinDistance(KissBall,oldball,line)
                            dis2 = math.sqrt((KissBall[0]-oldball[0])**2+(KissBall[1]-oldball[1])**2)
                            if dis2<dis1:
                                if oldball[0]<KissBall[0]:
                                    continue
                                else:
                                    if dis2 <= KissBall[2]+oldball[2]+2*NextBall[2] and dis2 >0:
                                        if oldball[0]<(KissBall[2]+2*NextBall[2]):
                                            CloseBall.append(oldball)
                            else:
                                if dis1 <= KissBall[2]+oldball[2]+2*NextBall[2] and dis1 >0:
                                    oldball=(oldball[0]+line,oldball[1],oldball[2])
                                    CloseBall.append(oldball)
                    #CloseBall.remove(KissBall)

                    if CloseBall:               #if CloseBall is not Null
                        #KissBall=(KissBall[0]-line,KissBall[1],KissBall[2])         #chage into x1= x1 -l
                        NextBall,StackBall =CloseBallIntersection(KissBall,CloseBall,NextBall,StackBallList,line)
                        NextBall = (round(NextBall[0],10),round(NextBall[1],10),round(NextBall[2],10))
                        StackBall = (round(StackBall[0],10),round(StackBall[1],10),round(StackBall[2],10))
                        # print('final position after stacking'+str(NextBall))
                        # print('Case 1 and there is closeball and Stacking \n')
                        if KissBall[0]<NextBall[0]+line<StackBall[0]:
                            rerotate = 0
                            anoKissBall = KissBall
                            unfinish,OverlapBall = CheckOverlapping(NextBall,CheckBallList,anoKissBall,StackBall,line)
                            if unfinish==1:
                                # NextBall=((KissBall[0]+NextBall[2]+KissBall[2])%line,KissBall[1],NextBall[2])
                                KissBall = OverlapBall
                                rerotate = 1
                                # print('Case 1 and there is closeball and Stable but overlapping \n')
                                finish = 0
                            else:
                                NextBall = (round(NextBall[0],10),round(NextBall[1],10),round(NextBall[2],10))
                                if (NextBall[1]-NextBall[2]<=0):
                                    NextBall = (round(NextBall[0],10),round(NextBall[2],10),round(NextBall[2],10))
                                # print('Case 1 and there is closeball and Stable no overlapping then finished \n')
                                finish =1
                        else:
                            # NextBall=((KissBall[0]+NextBall[2]+KissBall[2])%line,KissBall[1],NextBall[2])
                            #x of stackball is possible to be negative. so ensure it within (0,line)
                            KissBall = (StackBall[0]%line,StackBall[1],StackBall[2])
                            rerotate = 1
                            # print('case 1 Unstable and rerotate')
                            finish = 0

                    else:
                        rerotate = 0
                        StackBall = None
                        NextBall=((KissBall[0]+NextBall[2]+KissBall[2])%line,KissBall[1],NextBall[2])      #no need %line
                        if (NextBall[1]-NextBall[2]<=0):
                            NextBall=(NextBall[0]%line,NextBall[2],NextBall[2])
                            NextBall = (round(NextBall[0],10),round(NextBall[1],10),round(NextBall[2],10))
                            finish =1
                            # print('Case 1 No Closeball stand on ground')
                        else:
                            # print('Case 1 No Closeball and no stack do it again\n')
                            finish =0
     
                #Case 2 :this situation is that KissBall has part of leftside and NextBall hit on the Left Period Boundry
                elif (KissBall[0]-KissBall[2]<0 and NextBall[0]>=(line-KissBall[2]-NextBall[2])) or (NextBall[0]+NextBall[2]>line and KissBall[0]<(NextBall[2]+KissBall[2])):
                    TBallList.remove(KissBall)
                    for oldball in TBallList:
                        #print(oldball)
                        if oldball[0]>NextBall[0]:
                            continue
                        else:
                            #print('here '+str(oldball))
                            dis1 = PeriodPosiMinDistance(KissBall,oldball,line)
                            dis2 = math.sqrt((KissBall[0]-oldball[0])**2+(KissBall[1]-oldball[1])**2)
                            if dis2<dis1:
                                if oldball[0]>KissBall[0]:
                                    continue
                                else: 
                                    if dis2 <= KissBall[2]+oldball[2]+2*NextBall[2] and dis2>0:
                                        if oldball[0]>(line-KissBall[2]-2*NextBall[2]):
                                            CloseBall.append(oldball)
                            else:
                                if dis1 <= KissBall[2]+oldball[2]+2*NextBall[2] and dis1>0:
                                    oldball=(oldball[0]-line,oldball[1],oldball[2])
                                    CloseBall.append(oldball)
                    # if CloseBall:
                    #     for closeball in CloseBall:
                    #         if closeball[0] < (line-KissBall[2]-2*NextBall[2]):
                    #             CloseBall.remove(closeball)
                    #     print('Now after remove possible closeball'+str(CloseBall))
                    if CloseBall:               #if CloseBall is not Null
                        # print('possible Closeball'+str(CloseBall))
                        #KissBall=(KissBall[0]+line,KissBall[1],KissBall[2])                 #chage into x1= x1 -l

                        NextBall,StackBall =CloseBallIntersection(KissBall,CloseBall,NextBall,StackBallList,line)
                        NextBall = (round(NextBall[0],10),round(NextBall[1],10),round(NextBall[2],10))
                        StackBall = (round(StackBall[0],10),round(StackBall[1],10),round(StackBall[2],10))
                        # print('final position after stacking'+str(NextBall))               
                        # print('nextball, stackball'+str(NextBall)+str(StackBall))
                        # print('Case 2 there is closeball and Stacking\n')
                        if StackBall[0]<NextBall[0]-line<KissBall[0]:
                            rerotate = 0
                            anoKissBall = KissBall
                            unfinish,OverlapBall = CheckOverlapping(NextBall,CheckBallList,anoKissBall,StackBall,line)
                            if unfinish == 1:
                                # NextBall=((KissBall[0]-NextBall[2]-KissBall[2])%line,KissBall[1],NextBall[2])
                                KissBall=OverlapBall
                                rerotate =1
                                # print('Case 2 there is closeball and Stacking. Stable but overlapping\n')
                                finish =0
                            else:
                                NextBall = (round(NextBall[0],10),round(NextBall[1],10),round(NextBall[2],10))
                                if (NextBall[1]-NextBall[2]<=0):
                                    NextBall = (round(NextBall[0],10),round(NextBall[2],10),round(NextBall[2],10))
                                # print('Case 2 there is closeball and Stacking. Stable No overlapping finished\n')
                                finish = 1     
                        else:
                            # NextBall=((KissBall[0]-NextBall[2]-KissBall[2])%line,KissBall[1],NextBall[2])
                            KissBall = (StackBall[0]%line,StackBall[1],StackBall[2])
                            rerotate = 1
                            # print('Case 2 Unstable and rerotate')
                            finish =0
                    else:
                        rerotate = 0
                        StackBall = None
                        NextBall=((KissBall[0]-NextBall[2]-KissBall[2])%line,KissBall[1],NextBall[2])           #no need %line
                        if (NextBall[1]-NextBall[2]<=0):
                            NextBall=(NextBall[0]%line,NextBall[2],NextBall[2])
                            NextBall = (round(NextBall[0],10),round(NextBall[1],10),round(NextBall[2],10))
                            finish =1
                            # print('Case 2 No Closeball stand on ground')
                        else:
                            # print('Case 3 No Closeball and no stack do it again\n')
                            finish =0

                #Case 3 :this situation is that KissBall is a normal ball with line length. but it divided into 4 sub-situations
                else:
                    #one situation x3<x1 which would rotate to left

                    if NextBall[0]<KissBall[0]:
                        #first and second Situation. they are same
                        # print('Case 3 x3<x1')
                        TBallList.remove(KissBall)
                        for oldball in TBallList:
                            if (oldball[0]+oldball[2]>line and NextBall[0]<oldball[2]+NextBall[2]) or KissBall[0]-KissBall[2]-2*NextBall[2]<0:

                                dis1 = PeriodPosiMinDistance(KissBall,oldball,line) 
                                dis2 = math.sqrt((KissBall[0]-oldball[0])**2+(KissBall[1]-oldball[1])**2)
                                if dis2 < dis1:
                                    if oldball[0]>KissBall[0]:
                                        continue
                                    else:
                                        if dis2 <= KissBall[2]+oldball[2]+2*NextBall[2] and dis2>=0:
                                            CloseBall.append(oldball)
                                else:        
                                    if dis1 <=  KissBall[2]+oldball[2]+2*NextBall[2] and dis1 > 0:
                                        oldball=(oldball[0]-line,oldball[1],oldball[2])
                                        CloseBallSpecial.append(oldball)              

                            elif oldball[0]+oldball[2]<=line and KissBall[0]-KissBall[2]-2*NextBall[2]>=0:
                                #print(oldball)
                                if oldball[0]>KissBall[0]:
                                    continue
                                    #print(oldball)
                                else:
                                    dis = math.sqrt((KissBall[0]-oldball[0])**2+(KissBall[1]-oldball[1])**2)
                                    if dis <= KissBall[2]+oldball[2]+2*NextBall[2] and dis>=0:
                                        CloseBall.append(oldball)
                        #print('closeballspecial'+str(CloseBallSpecial))
                        #print('closeball'+str(CloseBall))
                        CloseBall.extend(CloseBallSpecial)

                        # print('all the possible closeball'+str(CloseBall))

                        CloseBall1 = CloseBall[:]
                        if CloseBall1:
                            NextBall,StackBall = CloseBallIntersection(KissBall,CloseBall,NextBall,StackBallList,line)
                            NextBall = (round(NextBall[0],10),round(NextBall[1],10),round(NextBall[2],10))
                            StackBall = (round(StackBall[0],10),round(StackBall[1],10),round(StackBall[2],10))
                            # print('final position after stacking'+str(NextBall))
                            # print('final Stackball'+str(StackBall))      
                            # print('Case 3 and x3<x1  then Stacking\n')

                            if StackBall[0]<NextBall[0]<KissBall[0] or StackBall[0]<NextBall[0]-line<KissBall[0]:
                                # print('there is Closeball. stable and within range\n')
                                rerotate = 0
                                anoKissBall = KissBall

                                unfinish,OverlapBall = CheckOverlapping(NextBall,CheckBallList,anoKissBall,StackBall,line)
                                if unfinish == 1:
                                    # NextBall=((KissBall[0]-NextBall[2]-KissBall[2])%line,KissBall[1],NextBall[2])
                                    KissBall = OverlapBall
                                    rerotate = 1
                                    # print('there is Closeball. stable and within range but overlapping. have to do it again\n')
                                    finish =0
                                else:
                                    NextBall = (round(NextBall[0],10),round(NextBall[1],10),round(NextBall[2],10))
                                    if (NextBall[1]-NextBall[2]<=0):
                                        NextBall = (round(NextBall[0],10),round(NextBall[2],10),round(NextBall[2],10))
                                    finish =1    

                            else:
                                # NextBall=((KissBall[0]-NextBall[2]-KissBall[2])%line,KissBall[1],NextBall[2])   
                                KissBall = (StackBall[0]%line,StackBall[1],StackBall[2])
                                rerotate = 1
                                # print('Case 3 x3<x1 there is Closeball.Unstable -Rerotate')
                                finish =0
                        else:
                            rerotate = 0
                            StackBall = None
                            NextBall=((KissBall[0]-NextBall[2]-KissBall[2])%line,KissBall[1],NextBall[2])                #i have deleted this %line
                            if (NextBall[1]-NextBall[2]<=0):
                                NextBall=(NextBall[0]%line,NextBall[2],NextBall[2])
                                NextBall = (round(NextBall[0],10),round(NextBall[1],10),round(NextBall[2],10))
                                finish =1
                                # print('Case 3 and x3<x1 No Closeball stand on ground')
                            else:
                                # print('Case 3 and x3<x1 No Closeball and no stack do it again\n')
                                finish =0

                    #another situation x3>x1which rotate to right
                    elif NextBall[0]>KissBall[0]:
                        #first and second situation. they are same
                        # print('Case 3 x3>x1')
                        TBallList.remove(KissBall)
                        for oldball in TBallList:
                            if (oldball[0]-oldball[2]<0 and NextBall[0]>(line-oldball[2]-NextBall[2])) or (KissBall[0]+KissBall[2]+2*NextBall[2])>line:
                                dis1 = PeriodNegMinDistance(KissBall,oldball,line)
                                dis2 = math.sqrt((KissBall[0]-oldball[0])**2+(KissBall[1]-oldball[1])**2)
                                if dis2<dis1:
                                    if oldball[0]<KissBall[0]:
                                        continue
                                    else:
                                        if dis2 <= KissBall[2]+oldball[2]+2*NextBall[2] and dis2 > 0:
                                            CloseBall.append(oldball)
                                else:
                                    if dis1 <=  KissBall[2]+oldball[2]+2*NextBall[2] and dis1 >0:
                                        oldball=(oldball[0]+line,oldball[1],oldball[2])
                                        CloseBallSpecial.append(oldball)

                            elif oldball[0]-oldball[2]>=0 and (KissBall[0]+KissBall[2]+2*NextBall[2])<=line:
                                if oldball[0]<KissBall[0]:
                                    continue
                                else:
                                    #normal one 
                                    dis = math.sqrt((KissBall[0]-oldball[0])**2+(KissBall[1]-oldball[1])**2)
                                    if dis <= KissBall[2]+oldball[2]+2*NextBall[2] and dis > 0:
                                        CloseBall.append(oldball)

                        CloseBall.extend(CloseBallSpecial)
    
                        CloseBall1 = CloseBall[:]
                        if CloseBall1:
                            NextBall,StackBall = CloseBallIntersection(KissBall,CloseBall,NextBall,StackBallList,line)
                            NextBall = (round(NextBall[0],10),round(NextBall[1],10),round(NextBall[2],10))
                            StackBall = (round(StackBall[0],10),round(StackBall[1],10),round(StackBall[2],10))
                            # print ('there is CloseBall and nextball, stackball'+str(NextBall)+str(StackBall))
                            
                            if KissBall[0]<NextBall[0]<StackBall[0] or KissBall[0]<NextBall[0]+line<StackBall[0]:
                                rerotate = 0
                                anoKissBall =KissBall
                                unfinish,OverlapBall = CheckOverlapping(NextBall,CheckBallList,anoKissBall,StackBall,line)
                                if unfinish == 1:
                                    # NextBall=((KissBall[0]+NextBall[2]+KissBall[2])%line,KissBall[1],NextBall[2])
                                    KissBall=OverlapBall
                                    rerotate = 1
                                    # print('there is Closeball.stable but overlapping and unfinished')
                                    finish =0
                                else:
                                    NextBall = (round(NextBall[0],10),round(NextBall[1],10),round(NextBall[2],10))
                                    if (NextBall[1]-NextBall[2]<=0):
                                        NextBall = (round(NextBall[0],10),round(NextBall[2],10),round(NextBall[2],10))
                                    # print('Case 3 and x3>x1 there is Closeball.stable Stacking finished\n')
                                    finish =1

                            else:
                                # NextBall=((KissBall[0]+NextBall[2]+KissBall[2])%line,KissBall[1],NextBall[2])                    # i have deleted this %line
                                KissBall = (StackBall[0]%line,StackBall[1],StackBall[2])
                                rerotate = 1
                                # print('Case 3 and x3>x1 there is Closeball.unstable After remove None -rerotate\n')
                                finish =0

                        else:
                            rerotate = 0
                            StackBall = None
                            NextBall=((KissBall[0]+NextBall[2]+KissBall[2])%line,KissBall[1],NextBall[2])                    # i have deleted this %line
                            if (NextBall[1]-NextBall[2]<=0):
                                NextBall=(NextBall[0]%line,NextBall[2],NextBall[2])
                                NextBall = (round(NextBall[0],10),round(NextBall[1],10),round(NextBall[2],10))
                                finish =1
                                # print('Case 3 and x3>x1 no CloseBall and stand on ground')
                            else:
                                # print('Case 3 and x3>x1 NoCloseBall and no stack do it again\n')
                                finish =0
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Packing Finished~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~                

        NextBall=(NextBall[0],NextBall[1],NextBall[2])                  #set the "next ball " to the final position 
        NextBall = (round(NextBall[0],10),round(NextBall[1],10),round(NextBall[2],10))
        BallList.append(NextBall)
        if NextBall[1]>1000:
            genGrain = 1
        

    # print('print all balls inforamtion with (x,y,r)') 
    # print(BallList)
    # print('\nFirst ball info: ' + str(ballFirst))
    # print('\nnow visulizing current simulation\n')


    ######################################visulization############################################

    visual=np.zeros((row,line))     
    #visual= CreateMatrix(row, line)
    for currentBall in BallList:
        #print(currentBall)
        visual= ballLocal(currentBall[0],currentBall[1],currentBall[2],row,line,visual)

    plt.imshow(visual, interpolation='nearest', origin='lower')

    # plt.savefig("/home/maoq/Desktop/debug/figure/"+ 'GrainPacking_with_'+str(BallNumbers)+'_Balls.png')
    
    # plt.savefig('GrainPacking_Balls.png')
    plt.show()

    GrainMatrix = [row[:] for row in visual]
    return GrainMatrix
    

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~here you run by itself~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#run in this file 
if __name__ == "__main__":
    
    start = time.time()
    # initial parameter essentially 
    gsdFunction = 2
    mean = 25
    std =  3
    minRad = 10
    maxRad = 25

    #Grain Packing
    GrainMatrix = GrainPacking(gsdFunction,mean,std,minRad,maxRad)
    # print(GrainMatrix)


    end = time.time()
    print('running time ' +str(end - start))