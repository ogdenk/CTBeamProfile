
import numpy
import math
from numpy import sin, cos, tan, pi, arctan, sqrt, arccos

# First take the CTDI measurements as a function of x
DIM = 21 # used 21 at first but that be any value
CTDI = [1,0.936,0.76,0.576,0.444,0.372,0.292]
dist = [0,5,10,15,20,25,30]


CTDIc = numpy.polyfit(dist, CTDI,4) # use a best fit to get a continuous function, get as many CTDI values as will have Beam Profile Values if the matrix approach is used
#d_locs = numpy.linspace(0,57,180)
# numpy.polyval(CTDIc,step*(57/180))
CTDIv = list() # CTDIv is a Vector with length DIM which can be set above
distances = list() # distances is a vector with length DIM and evenly divides the max dist into smaller distances
step = 0
while step < DIM:
    d_loc = step*(30.0/(DIM - 1)) # This value is set to 30 but could be anything, chose 30 to get the full range of the input distances
    distances.insert(step,step*(30.0/(DIM - 1)))
    CTDIv.insert(step, d_loc**4.0*CTDIc[0] + d_loc**3.0*CTDIc[1] + d_loc**2.0*CTDIc[2] + d_loc*CTDIc[3] + CTDIc[4])
    #print d_loc
    #print CTDIv[step]
    step = step +1

W = numpy.zeros((DIM,DIM)) # can call elements by W[1][2], can make this matrix as large as desired
dist_ad = numpy.zeros((DIM,DIM))
Theta_ad = numpy.zeros((DIM,DIM))

index0 = 0
index = 0
index1 = 0
index2 = 0
max_theta = 0.0
max_phi = 0.0
theta_arr_long = list()
Phi_arr_long = list() 
max_theta_v = list()

while index0 < DIM :
    while index < 1000: # Find the total theta change, angle between central beam and furthest out distance
        Phi = index*((pi)/1000)
        Rx = 57*cos(Phi)
        Ry = 57*sin(Phi)
        r_primex =  57*cos(Phi) - distances[index0]
        theta = arccos((Ry*Ry+(Rx*r_primex))/(sqrt(Rx*Rx + Ry*Ry)*sqrt(Ry*Ry+r_primex*r_primex)))
        Phi_arr_long.insert(index, Phi)
        theta_arr_long.insert(index, theta)
        if theta > max_theta: # From this we will find the maximum theta value as well as the Phi value which is associated with it
            max_theta = theta
            max_theta_v.insert(index0,theta) 
            max_phi = Phi
        index = index + 1

    #print 'max theta for a distance ', distances[index0], ' is ' , max_theta_v[index0]*180/pi 
    index = 0
    index0 = index0 + 1


#print 'max theta', max_theta*(180/pi)
#print 'max Phi', max_phi*(180/pi)

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
        if numpy.absolute(theta - arccos((Ry**2.0+(Rx*r_primex))/(sqrt(Rx**2.0 + Ry**2.0)*sqrt(Ry**2.0+r_primex*r_primex)))) < delta:
            delta = numpy.absolute(theta - arccos((Ry**2.0+(Rx*r_primex))/(sqrt(Rx**2.0 + Ry**2.0)*sqrt(Ry**2.0+r_primex*r_primex))))
            Phi_arr.insert(index,Phi) # find the Phi that goes with each theta in "theta_arr" which evenly samples the value of max_theta 
        index_sub = index_sub +1
    #print 'theta', theta_arr[index]*180/pi
    #print 'Phi', Phi_arr[index]*180/pi
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

while index1 < DIM: # index1 is for each row where we consider a fixed distance d, We are considering d to be zero at the origin and going in the positive x axis
    d = distances[index1]    
    while index2 < DIM: # index2 is for each column where we consider a fixed theta. Contributions from different theta's (angles of the BP) will be multiplied by the BP vector to yield the CTDIv vector,
                        #       DIM values was used, this will give us DIM points of the BP. 
                        #       We are considering d to be zero at the origin and going in the positive x axis, for a location on the positive x-axis a Phi of zero is used. 
                        #       Because of symmetry the BP and W matrix will only use theta from 0 to max_theta, max_theta is the angle which contributes to the furthest d value from the origin.  
        while index_sub <100: 
            Phi = index_sub*(pi/2)/100 # Angle between (vector from the origin to the focal spot) and (Vector from origin to positive X axis also knonw as zero degree position), chose 90 degrees because the max_theta should be less than this so any other theta would be as well
            Rx = 57*cos(Phi)   # x component of R
            Ry = 57*sin(Phi)   # y component of R, this is the same as r_primey
            r_primex =  57*cos(Phi) - d # x component of the vector from the origin to the focal spot, r_primey would be the same as Ry
            theta = arccos((Ry**2.0+(Rx*r_primex))/(sqrt(Rx**2.0 + Ry**2.0)*sqrt(Ry**2.0+r_primex*r_primex)))
            # theat is the angle between (vector from the origin to the focal spot) and (Vector from location d on the x axis where CTDI measurement was made and the focal spot), this angle is the BP angle
            #print index2
            #print theta_arr[index2]
            #print arccos((Ry*Ry+(Rx*r_primex))/(sqrt(Rx*Rx + Ry*Ry)*sqrt(Ry*Ry+r_primex*r_primex))  )          
            if numpy.absolute(theta_arr[index2] - arccos((Ry**2.0+(Rx*r_primex))/(sqrt(Rx**2.0 + Ry**2.0)*sqrt(Ry**2.0+r_primex*r_primex)))) < delta:
                # The purpose of this If statement and this 3rd while loop it to find the Phi associated with each theta (column) and distance d (row), brute force method used as can be seen with minimizing delta
                delta = numpy.absolute(theta_arr[index2] - arccos((Ry**2.0+(Rx*r_primex))/(sqrt(Rx**2.0 + Ry**2.0)*sqrt(Ry**2.0+r_primex*r_primex))))
                Phi_new = Phi # We now have the Phi used for the theta column position listed by index2 and distance in index1
                theta_new = theta # Theta new is the theta associated with the Phi above
            index_sub = index_sub + 1

        rprime = sqrt(57*sin(Phi_new)*57*sin(Phi_new) + (57*cos(Phi_new) - d)*(57*cos(Phi_new) - d)) # vector from the CTDI measurement to the focal spot
        rprime_dot = (2*57**2*sin(Phi_new)*cos(Phi_new)-2*57*sin(Phi_new)*(57*cos(Phi_new)-d))/(2*rprime) # r prime time derivative, rate at which r prime changes with time
        theta_dot = (d*cos(Phi_new) - rprime_dot*sin(theta_new))/(rprime*cos(theta_new)) # rate at which theta changes with time 

        #if index1 == 10:
        #    print 'index2 = ' , index2, ', Phi = ', Phi_new, ', theta = ', theta_new 'theta =' , theta,
        if index1 == 0:
            print ' theta arr =', theta_arr[index2], ' theta arr in degree =', theta_arr[index2]*180/pi # these are the 21 theta's used

        if theta_arr[index2] <= max_theta_v[index1]: #Each row goes to the max possible theta, this is only achievable for the largest d or the last row of the matrix, each other row will have a max theta that differs and dependa on the distance d            
            dist_ad[index1][index2] = (57**2)/(math.pow((57*sin(Phi_new)),2) +math.pow((57*cos(Phi_new)-d),2)) #weighting factor that takes into account inverse square to give a BP at a fixed distance R
            # the contributions wrt theta_dot will not be linear but inverse since we want to know how long the focal spot stays at a theta position, also although there will be negative values of theta_dot to show change in
            #   direction we only care about the magnitude
            Theta_ad[index1][index2] = abs(1/theta_dot) # the weighting factor based on how long a angle theta contributes 
            W[index1][index2] = ((57**2)/(math.pow((57*sin(Phi_new)),2) +math.pow((57*cos(Phi_new)-d),2)))*abs(1/theta_dot)
        else:
            W[index1][index2] = 0.0

        delta = 1000
        Phi_new = 0
        theta_new = 0
        index2 = index2 + 1
        index_sub =0

    index1 = index1 + 1
    index2=0
    print 'index1 = ', index1

print 'W matrix'
print('\n'.join([''.join(['{:4}'.format(round(item,4)) for item in row]) 
      for row in W]))
print  ' Adjust for dist'
print('\n'.join([''.join(['{:4}'.format(round(item,4)) for item in row]) 
      for row in dist_ad]))
print  ' Adjust for theta'
print('\n'.join([''.join(['{:4}'.format(round(item,4)) for item in row]) 
      for row in Theta_ad]))

# Go through each row and normalize
index1 = 0
index2 = 0
Sum = 0
W[0][0] = 1
while index1 < DIM:
   
    while index2 < DIM:
        Sum = Sum + W[index1][index2]
        index2 = index2 +1

    index2 = 0
    while index2 < DIM:
        W[index1][index2] = W[index1][index2]/Sum
        index2 = index2 +1
        
    index2=0
    Sum = 0
    index1 = index1 + 1

index1 = 0
index2 = 0
while index1 < DIM:
    while index2 < DIM:
        if W[index1][index2] == 0.0:
            if index1 == 0.0:
                W[index1][index2] = 0.0 # 0.0001
        index2 = index2+1
    index2 = 0
    index1 = index1 + 1

print 'Normalized W'
print('\n'.join([''.join(['{:4}'.format(round(item,4)) for item in row]) 
      for row in W]))

#BP = numpy.linalg.solve(W,CTDIv)
BP = numpy.linalg.lstsq(W,CTDIv, 0.001) # http://docs.scipy.org/doc/numpy/reference/generated/numpy.linalg.lstsq.html 

print 'Beam profile using least squares'
index = 0
while index < len(BP):
    print BP[index]
    index = index + 1

print 'CTDI interpolated'
index = 0
while index < len(CTDIv):
    print CTDIv[index]
    index = index + 1


    



    



    


