
import numpy
from numpy import sin, cos, tan, pi, arctan
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


# Use a List for the measured CTDI values, in next version allow variable distance for now used with fixed 5 cm values

print ("Opening and closing file")
f = open("Num32.txt","r")
lines = f.readlines() # Put all the lines into an array
f.close()

index = 0
place = 0
num = list()
time = list()
rate = list()

for line in lines:
    #print line
    #print line.find('\t')
    num.insert(index, line[0:line.find('\t')-1])
    #print num[index]
    line = line[line.find('\t')+1:]
    #print line.find('\t')
    time.insert(index, line[0:line.find('\t')-1])
    #print time[index]
    line = line[line.find('\t')+1:]
    #print line.find('\t')
    rate.insert(index, line)
    #print rate[index]
    index = index +1
    
# Data imported

#print('Please enter the radius of the path that the detector takes to make one rotation, in cm')
#R = float(input())
#print('Please enter the distance from the center of rotation that the data was measured at, in cm ')
#r = float(input())
R = 60.0
#R = float(R)
r = 15.0
#r = float(r)

rot_int = 2*pi/len(time) # In radians
rot_ang = pi/2 # Because of geometry the largest exposure rate starts at t=0, Phi = 90, Theta=0
Phi = list()
index = 0
for i in time:
    Phi.insert(index,rot_ang)
    rot_ang = rot_ang - rot_int
    #print Phi[index]*180/pi, arctan((r*cos(pi/2-Phi[index])/(R-r*cos(pi/2-Phi[index])))*tan(pi/2-Phi[index]))*180/pi
    index = index +1
#print "The above is Phi of length: "
#print len(Phi)
# Now we have a Phi for each list value, this is the angle (in radians) that the exposure is measures at wrt the center of rotation

Theta = list()
index = 0
calc = 0.0
for i in Phi:
    calc = arctan((r*cos(pi/2-Phi[index])/(R-r*cos(pi/2-Phi[index])))*tan(pi/2-Phi[index]))
    Theta.insert(index,calc)
    #print Theta[index]*180/pi
    index = index + 1
#print "The above is Theta of length: "
#print len(Phi)
# Now we have a Theta for each list value, this is the angle (in radians) wrt the central beam that contributes to location "r" 

Beam_profile = list()
index = 0
calc = 0.0
for i in Theta:
     rate.insert(index, float(rate[index]))
     calc = rate[index]*(R**2+r**2-2*r*R*sin(Phi[index]))/(R**2)
     Beam_profile.insert(index,calc)
     #print Beam_profile[index]
     index = index + 1
#print "The above is beam profile of length: "
#print len(Beam_profile)
# Beam Profile calculated for each list value

results = open("Results2.txt","w")
index = 0
res_out = ""
for i in Theta:
    res_out = ""+str(Theta[index])+ "\t"+ str(Beam_profile[index])+ "\n"
    results.write(res_out)
    index = index +1
results.close()
# return results in terms of Theta and Beam Profile

red_patch = mpatches.Patch(color='red', label='Beam Profile')
plt.legend(handles=[red_patch])

plt.plot(Theta,Beam_profile, 'r', ms=1)
plt.show()


#Dist = [0,5,10,15,20,25]
#CTDI = [1,0.92,0.68,0.37,0.23,0.11]
#Beam_profile = [1,0.87,0.59,0.29,0.17,0.12]
#R = 30 # R in cm
