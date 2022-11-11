"""
Calculating the ephemeris of a celestial body

"""
#############################################################
## IMPORTS

from Vector import Vector
from Matrix import Matrix
import math
import os
import time
from NR import new_rap

#############################################################
## GLOBAL CONSTANTS

# Meta usage
bar = "*************************************************************************"

# Celestial instance constants
# 23.439 at J2000
# 23.436 at 2022
eps = 23.439*math.pi/180    #rad
J2000 = 2451545.0

# Physical constants
G = 6.674e-11 # m**3/kg/s**2
solar_mass = 1.98847e30 # kg
k = 0.01720209895   # rad/day
au = 1.495978707e11 # m
c =  2.99792458e8 # m/s

# Lambda funcs
clear = lambda: os.system('cls')

#############################################################
## FUNCTIONS

############################################
## PRINT FORMATTING

# encloses text with a bar
def enclose(text, f = 0):
    if f == 0:
        print(text)
        print(bar)
    else:
        f.write(text)
        f.write("\n")
        f.write(bar)
        f.write("\n")
    
""" 
converts radians to degrees, arcmins, and arcsecs
60' = 1deg, 60" = 1'
"""
def rad_to_deg(ang):
    
    ang = ang*180/math.pi   # converts to degrees
    
    deg = int(ang)          # floors the degree
    res = abs(ang - deg)
    
    min = int(res*60)
    
    res = abs(res-min/60)

    sec = int(res*60*60)
    
    return f"{deg}◦ {min}\' {sec}\""

"""
converts radians to hours, mins, and secs
(one revolution is 24 hours)
"""
def rad_to_time(ang):

    ang = ang*180/math.pi# converts to degrees
    
    # positive angles
    if ang < 0:
        ang = 360 + ang
    
    h = int(ang/360*24)
    
    res = ang - h*360/24
    
    m = int(res/360*24*60)
    
    res = res - m*360/24/60
    
    s = int(res/360*24*3600+0.5)
    
    return f"{h}h {m}m {s}s"

############################################
## INPUT

""" 
handles input of celestial object
data for 6918 Manaslu is saved
"""
def choose_object(celestial_info):
    while True:
        print(bar)
        enclose("use data for 6918 Manaslu?\nyes: y\tno: n")
        choice = input("input: ")
        
        if choice == "y":
            elements = [2.40726, 0.1410730,\
                math.radians(1.86650),math.radians(346.94588),\
                    math.radians(183.43911)]
            name = "6918 Manaslu"
            celestial_info.append(name)
            celestial_info.append(elements)
            
            # since the closest passage is after 22 feb 2022
            # we deduct a full period from it
            celestial_info.append(2459930.4545718)
            
            celestial_info.append("20220324")
            break
        
        elif choice == "n":
            name = input("name of object: ")
            clear()
            enclose(f"{name}-properties")
            elements = define_elements()
            celestial_info.append(name)
            celestial_info.append(elements)
            celestial_info.append(float(input("time since last passage [JD]: ")))
            celestial_info.append(input("time at calc [yyyymmdd]: "))
            break
        clear()
    
    print(bar)
    print("orbital elements set")
    time.sleep(1)
    clear()
        
""" 
sets the position of the Earth
if exact position is not known, approximate
position is calculated using approx_earth_pos()
"""
def set_earth(celestial_info, t):
    pos_earth = 0
    
    # position of the earth (heliocentric ecliptic)
    # at 2022/03/24
    if celestial_info[3] == "20220324":
        X = -9.953944023483982e-1
        Y = -5.260671184560772e-2
        Z = 9.130349223770484e-6
        pos_earth = Vector(X,Y,Z)
        
    else:
        while True:
            print("provide accurate position for Earth?")
            enclose("yes: y\t\tno: n")
            choice = input("input: ")
            
            if choice == "y":
                pos_earth = define_earth_pos()
                break
            elif choice == "n":
                pos_earth = approx_earth_pos(t)
                break
            
            clear()
        clear()
        print("Earth position set")
        time.sleep(1)
        
    return pos_earth

############################################
## PRINTING / WRITING

""" 
prints the data in an attractive manner
"""
def print_data(celestial_info, positions, corr_pos, calcs,\
    corr_calcs, angles, corr_angles, corr):
    print(bar)
    enclose(celestial_info[0])
    print("Ephemeris")
    enclose(f"{celestial_info[3][:4]}/{celestial_info[3][4:6]}/{celestial_info[3][6:]}\
        \t<-->\t{nrml_to_JDT(celestial_info[3])} [JD]")
    print("Topocentric correction was not made")
    enclose(f"Aberration correction in minutes: {corr}")
    print("EQUATORIAL ANGLES")
    print(f"Angle\t\t\tUncorrected\tCorrected")
    print("-------------------------------------------------------------------------")
    print(f"Declination: \t\t{rad_to_deg(angles[0])} \t{rad_to_deg(corr_angles[0])}")
    print(f"Right ascension: \t{rad_to_time(angles[1])} \t{rad_to_time(corr_angles[1])}")
    
    print(bar)
    print("POSITION")
    print(f"{'System':15}{'Uncorrected':<30}{'Corrected':<25}")
    print("-------------------------------------------------------------------------")    
    
    for pos, cpos, labl in zip(positions, corr_pos, ["orbital:", "ecliptic:", "equatorial:"]):
        print(f"|{labl:15}{str(pos.rnd()):<30}{str(cpos.rnd()):<25} | [au]")
        
    print(bar)
    print("CALCULATED PROPERTIES")
    print(f"{'Property':17}{'Uncorrected':<25}{'Corrected':<25}")
    print("-------------------------------------------------------------------------")
    for elm, celm, labl in zip(calcs, corr_calcs, ["P [yrs]:","n [deg]:", "t [JD]:",\
        "M [deg]:", "E [deg]:", "r [au]:", "Delta [au]:"]):
        print(f"|{labl:15} {elm:<25} {celm:<25}|")
    print(bar)
    
""" 
writes the data in an attractive manner
(not attractive code though)
"""
def write_data(celestial_info, positions, corr_pos, calcs,\
    corr_calcs, angles, corr_angles, corr, f):
    space = lambda : f.write("\n")
    f.write(bar)
    space()
    enclose(celestial_info[0], f)
    f.write("Ephemeris")
    space()
    enclose(f"{celestial_info[3][:4]}/{celestial_info[3][4:6]}/{celestial_info[3][6:]}\
        \t<-->\t{nrml_to_JDT(celestial_info[3])} [JD]", f)
    f.write("Topocentric correction was not made")
    space()
    enclose(f"Aberration correction in minutes: {corr}", f)
    f.write("EQUATORIAL ANGLES")
    space()
    f.write(f"Angle\t\t\tUncorrected\tCorrected")
    space()
    f.write("-------------------------------------------------------------------------")
    space()
    f.write(f"Declination: \t\t{rad_to_deg(angles[0])} \t{rad_to_deg(corr_angles[0])}")
    space()
    f.write(f"Right ascension: \t{rad_to_time(angles[1])} \t{rad_to_time(corr_angles[1])}")
    space()
    
    f.write(bar)
    space()
    f.write("POSITION")
    space()
    f.write(f"{'System':15}{'Uncorrected':<30}{'Corrected':<25}")
    space()
    f.write("-------------------------------------------------------------------------")    
    space()
    
    for pos, cpos, labl in zip(positions, corr_pos, ["orbital:", "ecliptic:", "equatorial:"]):
        f.write(f"|{labl:15}{str(pos.rnd()):<30}{str(cpos.rnd()):<25} | [au]")
        space()
        
    f.write(bar)
    space()
    f.write("CALCULATED PROPERTIES")
    space()
    f.write(f"{'Property':17}{'Uncorrected':<25}{'Corrected':<25}")
    space()
    f.write("-------------------------------------------------------------------------")
    space()
    for elm, celm, labl in zip(calcs, corr_calcs, ["P [yrs]:","n [deg]:", "t [JD]:",\
        "M [deg]:", "E [deg]:", "r [au]:", "Delta [au]:"]):
        f.write(f"|{labl:15} {elm:<25} {celm:<25}|")
        space()
    f.write(bar)
    

############################################
## DEFINING ELEMENTS

""" 
handles the input for the orbital elements
"""
def define_elements():
    elements = []
    for elm in ["semimajor axis [au]", "eccentricity [N/A]", \
        "inclination [deg]", "argument of the periapsis [deg]",\
            "longitude of the ascending node [deg]"]:
        elements.append(float(input(f"{elm}: ")))
    
    elements[2] = math.radians(elements[2])
    elements[3] = math.radians(elements[3])
    elements[4] = math.radians(elements[4])
    
    return elements

""" 
handles the input for the position of the earth
(if accurate position is avaliable)
"""
def define_earth_pos():
    pos = []
    for p in ["x-position [au]", "y-position [au]", "z-position [au]"]:
        pos.append(float(input(f"{p}: ")))
    return Vector(pos[0], pos[1], pos[2])

""" 
defines the coordinate transformation matrix
from orbital coordinates to ecliptic
(both heliocentric)
"""
def R_orb2ec(elements):
    i = elements[2]
    w = elements[3]
    W = elements[4]
    cos = math.cos
    sin = math.sin
    R11 = cos(W)*cos(w)-sin(W)*cos(i)*sin(w)
    R12 = -cos(W)*sin(w)-sin(W)*cos(i)*cos(w)
    R13 = sin(W)*sin(i)
    R21 = sin(W)*cos(w)+cos(W)*cos(i)*sin(w)
    R22 = -sin(W)*sin(w)+cos(W)*cos(i)*cos(w)
    R23 = -cos(W)*sin(i)
    R31 = sin(i)*sin(w)
    R32 = sin(i)*cos(w)
    R33 = cos(i)
    R = Matrix("R_orb2ec", R11, R12, R13, R21, R22, R23,\
        R31, R32, R33)
    return R

""" 
defines the coordinate transformation matrix
from ecliptic coordinates to equatorial
(both heliocentric, rotation with eps)
"""
def R_ec2eq(deg):
    return Matrix("R_ec2eq", 1, 0, 0, 0, math.cos(deg),\
        -math.sin(deg), 0, math.sin(deg), math.cos(deg))

############################################
## CALCULATIONS

""" 
calculates the Julian date using a string
of the date, format 'yyyymmdd'
credits to Duffett-Smith
"""
def nrml_to_JDT(date, corr=0):
    Y = int(date[0:4])
    M = int(date[4:6])
    D = int(date[6:].strip())
    # TT is assumed to be 00 hrs
    
    JD = 367*Y - int(7/4*(Y + int((M+9)/12)))
    JD += int(275*M/9)
    JD += D + 1721013.5
    JD -=corr/24

    return JD

def calc_period(elements):
    return 2*math.pi*math.sqrt(elements[0]**3)/k

def mean_motion(period):
    return 2*math.pi/period

def mean_anamoly(n, t, T):
    return n*(t-T)

""" 
calculates the orbital coordinates (heliocentric)
z = 0 due to plane of motion
"""
def calc_orb_pos(celestial_info, calcs):
    x = celestial_info[1][0]*(math.cos(calcs[4])-celestial_info[1][1])
    
    y = celestial_info[1][0]*math.sqrt(1-celestial_info[1][1]**2)*math.sin(calcs[4])
    
    return Vector(x,y,0)

""" 
approximates the position of the Earth
t in JD after the standard epoch
"""
def approx_earth_pos(t):
    Tstar = t - J2000
    
    gamma = 357.528 + 0.9856003*Tstar
    
    lamb = 280.460 + 0.9856474*Tstar 
    lamb += 1.915*math.sin(math.radians(gamma)) + 0.020*math.sin(2*math.radians(gamma))
    
    lamb = math.radians(lamb)
    
    gamma = math.radians(gamma)
    
    r = 1.00014 - 0.01671*math.cos(gamma)
    r -= 0.00014*math.cos(2*gamma)
    
    return Vector(-r*math.cos(lamb),-r*math.sin(lamb),0)

""" 
handles the calculation of the orbital position,
ecliptic position, as well as the period P, average motion n,
mean anomaly M, eccentric anomaly E
"""
def calc_pos(celestial_info, postitions, calcs):
    # P
    calcs[0] = calc_period(celestial_info[1])

    # n
    calcs[1] = mean_motion(calcs[0])
    
    # M
    calcs[3] = mean_anamoly(calcs[1], calcs[2], celestial_info[2])

    # E, using Newton-Raphson
    calcs[4] = new_rap(10,calcs[3],celestial_info[1][1])
    
    postitions[0] = calc_orb_pos(celestial_info, calcs)
    
    R = R_orb2ec(celestial_info[1])

    postitions[1] = R*postitions[0]
    
""" 
calculates the equatorial angles
"""
def calc_eq_ang(pos_eq):
    # declination
    decl = math.asin(pos_eq.z/abs(pos_eq))
    
    # right ascension
    alpha = math.atan2(pos_eq.y, pos_eq.x)
    
    # positive angles
    if alpha < 0:
        alpha = 2*math.pi + alpha

    return [decl, alpha]

""" 
handles the calculation of 
"""
def handling_calcs(celestial_info, positions, pos_earth, calcs):
    
    calc_pos(celestial_info, positions, calcs)
    
    # rotational matrix, ecliptic to equatorial
    R2 = R_ec2eq(eps)
    
    # Delta
    calcs[6] = abs(positions[1]-pos_earth)
    
    # equatorial coordinates (geocentric)
    positions[2] = R2*(positions[1]-pos_earth)
    
    calcs[5] = abs(positions[0])
    
    return calc_eq_ang(positions[2])

""" 
converts the calculated elements to degrees
and years
"""
def conv_to_disp(calcs):
    temp = []
    temp.append(calcs[0]/365)
    temp.append(math.degrees(calcs[1]))
    temp.append(calcs[2])
    temp.append(math.degrees(calcs[3]))
    temp.append(math.degrees(calcs[4]))
    temp.append(calcs[5])
    temp.append(calcs[6])
    
    return temp

#############################################################
## MAIN
    
def main():
    # [name, orbital elements, closest perihelion passage,
    # date of calculation]
    celestial_info = [] 
    
    # [orb, ec, eq]
    positions = [0,0,0]
    
    # [orb, ec, eq]
    corr_pos = [0,0,0]
    
    # [P, n, t, M, E, r, Delta]
    calcs = [0,0,0,0,0,0,0]
    
    # [P, n, t, M, E, r, Delta]
    corr_calcs = [0,0,0,0,0,0,0]

    choose_object(celestial_info)
    
    # t
    calcs[2] = nrml_to_JDT(celestial_info[3])    
    
    # if the time T is after the observation date
    # remove an entire period
    if calcs[2] < celestial_info[2]:
        celestial_info[2] -= calc_period(celestial_info[1])
    
    pos_earth = set_earth(celestial_info, calcs[2])
    
    angles = handling_calcs(celestial_info, positions, pos_earth, calcs)
    
    print("correcting for planetary abberation...")
    
    # mins
    corr = calcs[6]*au/c/60
    
    # redefines the observation date due to
    # finite light propogation
    # planetary aberration
    corr_calcs[2] = calcs[2] - corr/60/24
    
    corr_angles = handling_calcs(celestial_info, corr_pos, pos_earth, corr_calcs)
    
    clear()
    
    calc_conv = conv_to_disp(calcs)
    corr_conv = conv_to_disp(corr_calcs)
    
    print_data(celestial_info, positions, corr_pos, calc_conv, corr_conv, angles, corr_angles, corr)
    
    input("press [ENTER] to exit...")
    
    print(f"saving to {celestial_info[0]}.txt")
    time.sleep(1)
    
    with open(f"{celestial_info[0]}.txt", "w") as f:
        write_data(celestial_info, positions, corr_pos, calc_conv, corr_conv, angles, corr_angles, corr, f)
        
#############################################################    
## RUN CODE

if __name__ == "__main__":
    main()
