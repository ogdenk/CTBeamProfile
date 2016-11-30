
import numpy
import math
import array
from numpy import sin, cos, tan, pi, arctan, sqrt, arccos

#DONE -  First take the CTDI measurements as a function of x
# Done - use a best fit to get a continuous function, get the CTDI for each degree, this will be Vector "CTDI" with 360 enries
# find the time per degree for each location x that CTDI measured at (Create a function that gives the weighting factor for
#                  each theta in the beam profile, says ho much that beam profile angle contributes to the total CTDI) This is Matrix W    
# Make a Vector named "BP" for each degree of the BP, can calculate, if radius is 57 then arctan((57/2)/57) is the max angle and can go 2x this
# Linear algebra solve for W*BP = CTDI  [ looks like CTDI vector and BP vector need to be the same in length, could do this by making
#           CTDI vector 360 and the BP sub degree to make 360)
# This would give the BP

# First take the CTDI measurements as a function of x
DIM = 21 # used 21 at first
CTDI = [1,0.936,0.76,0.576,0.444,0.372,0.292]
dist = [0,5,10,15,20,25,30]

# use a best fit to get a continuous function, get as many CTDI values as will have Beam Profile Values (Goal is one for each degree), this will be Vector "CTDI" with 360 enries
CTDIc = numpy.polyfit(dist, CTDI,4)
#d_locs = numpy.linspace(0,57,180)
# numpy.polyval(CTDIc,step*(57/180))
CTDIv = list()
distances = list()
step = 0
while step < DIM:
    d_loc = step*(30.0/(DIM - 1))
    distances.insert(step,step*(30.0/(DIM - 1)))
    CTDIv.insert(step, (d_loc*d_loc*d_loc*d_loc)*CTDIc[0] + (d_loc*d_loc*d_loc)*CTDIc[1] + (d_loc*d_loc)*CTDIc[2] + (d_loc)*CTDIc[3] + CTDIc[4])
    #print d_loc
    #print CTDIv[step]
    step = step +1

# CTDIv is a Vector with 21 components
# Now to make a 2D matrix that will be multiplied by a BP (Beam Profile) vector to yield the CTDIv Vector
#          Each row will consist of the contributions that each angle (definded in the BP vector) of the BP contributes to the CTDIv at a sepefic distnace d
# The BP vector will only depend on angle theta so the matrix will have to account for distance and relative time that angle contributes to the CTDI

#Creating distance correction

W = numpy.zeros((DIM,DIM)) # can call eliments by W[1][2], can make this as large as desired
index = 0
index1 = 0
index2 = 0
max_theta = 0.0
max_phi = 0.0
theta_arr_long = list()
Phi_arr_long = list() 
  
while index < 1000: # 1) Find the total theta change
    Phi = index*((pi)/1000)
    Rx = 57*cos(Phi)
    Ry = 57*sin(Phi)
    r_primex =  57*cos(Phi) - dist[6]
    theta = arccos((Ry*Ry+(Rx*r_primex))/(sqrt(Rx*Rx + Ry*Ry)*sqrt(Ry*Ry+r_primex*r_primex)))
    Phi_arr_long.insert(index, Phi)
    theta_arr_long.insert(index, theta)
    if theta > max_theta:
        max_theta = theta
        max_phi = Phi
    index = index + 1

print 'max theta', max_theta*(180/pi)
print 'max Phi', max_phi*(180/pi)

index = 0
index_sub = 0
delta = 1000
theta = 0.0
Phi = 0.0
theta_arr = list()
Phi_arr = list() 
while index < DIM: # go through each theta (to be used to make the matrix W columns) and calc Phi
    theta = index*(max_theta)/(DIM - 1) # equally break theta into 21 parts
    theta_arr.insert(index,theta)
    while index_sub < 1000:
        Phi = index_sub*((max_phi)/1000)
        Rx = 57*cos(Phi)
        Ry = 57*sin(Phi)
        r_primex =  57*cos(Phi) - dist[6]
        if numpy.absolute(theta - arccos((Ry*Ry+(Rx*r_primex))/(sqrt(Rx*Rx + Ry*Ry)*sqrt(Ry*Ry+r_primex*r_primex)))) < delta:
            delta = numpy.absolute(theta - arccos((Ry*Ry+(Rx*r_primex))/(sqrt(Rx*Rx + Ry*Ry)*sqrt(Ry*Ry+r_primex*r_primex))))        
            Phi_arr.insert(index,Phi) # find the Phi that goes with each theta
        index_sub = index_sub +1
    print 'theta', theta_arr[index]*180/pi
    print 'Phi', Phi_arr[index]*180/pi
    delta = 1000
    index_sub = 0
    index = index + 1

Phi = 0
Phi_new = 0.0
theta_new = 0.0
delta = 1000
index1 = 0
index2 = 0
index_sub = 0

while index1 < DIM: #180: # index1 is for each row where we consider a fixed d, We are considering d to be going in the positive x axis
    d = distances[index1]
    
    while index2 < DIM: # index2 is for each column where we consider a theta and look at the contributions from different theta's (angles of the BP),
                        #       21 was used aribrarily, this will give us 21 points of the BP. 
                        #       We are considering d to be going in the positive x axis and Phi of zero to be in the positive x axis.
                        #       Because of symetrry we will go from 0 to max_theta 
        while index_sub <100:
            Phi = index_sub*(pi/2)/100 # Angle between (vector from the origen to the focal spot) and (Vector from origin to positive X axis also knonw as 90 degree position)
            Rx = 57*cos(Phi)   # x component of R, this is the same as r_prime_x
            Ry = 57*sin(Phi) 
            r_primex =  57*cos(Phi) - d
            theta = arccos((Ry*Ry+(Rx*r_primex))/(sqrt(Rx*Rx + Ry*Ry)*sqrt(Ry*Ry+r_primex*r_primex))) # Angle between (vector from the origen to the focal spot) and (Vector from location d on the y axis where CTDI measurement made and the focal spot), this angle is the BP angle
            #print index2
            #print theta_arr[index2]
            #print arccos((Ry*Ry+(Rx*r_primex))/(sqrt(Rx*Rx + Ry*Ry)*sqrt(Ry*Ry+r_primex*r_primex))  )          
            if numpy.absolute(theta_arr[index2] - arccos((Ry*Ry+(Rx*r_primex))/(sqrt(Rx*Rx + Ry*Ry)*sqrt(Ry*Ry+r_primex*r_primex)))) < delta:
                delta = numpy.absolute(theta_arr[index2] - arccos((Ry*Ry+(Rx*r_primex))/(sqrt(Rx*Rx + Ry*Ry)*sqrt(Ry*Ry+r_primex*r_primex))))        
                Phi_new = Phi # We now have the Phi used for the theta column position listed by index2 and diatnce in index1
                theta_new = theta
            index_sub = index_sub + 1

        #if index1 == 10:
        #    print 'index2 = ' , index2, ', Phi = ', Phi_new, ', theta = ', theta_new
        if index1 == 0:
            print 'theta =' , theta, ' theta arr =', theta_arr[index2]
            
        
        if theta_arr[index2] <= theta:           
            W[index1][index2] = (57*57)/(math.pow((57*sin(Phi_new)),2) +math.pow((57*cos(Phi_new)-d),2))
            if index1 == 0:
                print 'assigned and theta =' , theta, ' theta arr =', theta_arr[index2]            
        else:
            W[index1][index2] = 0

        delta = 1000
        Phi_new = 0
        theta_new = 0
        index2 = index2 + 1
        index_sub =0

# New plan needed, BP vector goes from +arctan(30/57) to -arctan(30/57), however each row of W will not, each row will go from +arctan(d/57) to -arctan(d/57)
#      Maybe when makeing each row W, have each column be for a different degree (or subdegree of theta), this way would not have to calculate theta
#      for each column would need theta (known), dtheta/dt (can calculate), what to do with Phi then (Could calculate Phi based off theta with dot product Cos(90-Phi) = (Rv dot dv)/(R^2 +d^2) 

    index1 = index1 + 1
    index2=0
    print 'index1 = ', index1
    #print ""
    #print index1


print('\n'.join([''.join(['{:4}'.format(round(item,4)) for item in row]) 
      for row in W]))

BP = numpy.linalg.solve(W,CTDIv)

index = 0
while index < len(BP):
    print BP[index]
    index = index + 1



# Next to do is use linear algebra to solve
# then take into account the dtheta/dt

    



    



    


