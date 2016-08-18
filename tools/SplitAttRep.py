#!/usr/bin/python
##############################################################################
#                                                                            #
# Simple script to separate attractive and repulsive parts of the NCI plot   #
# in VMD                                                                     #
#                                                                            #
# Usage:                                                                     #
#   user:$ python SplitAttRep.py filename                                    #
#                                                                            #
# Here 'filename' should match filename-dens.cube and filename-grad.cube     #
#                                                                            #
##############################################################################
# Author: Eric G. Kratz, Wayne State University

#Import modules
import sys

#Check for extra arguments
if (len(sys.argv) != 2):
  line = '\n'
  line += "Usage:"
  line += '\n'
  line += '\n'
  line += "  user:$ python SplitAttRep.py filename"
  line += '\n'
  line += '\n'
  line += "NB: 'filename' should match "
  line += "filename-dens.cube and filename-grad.cube"
  line += '\n'
  print(line)
  exit(0)

#Extract filename label for cube files
try:
  cubename = sys.argv[1]
  cubename = cubename.strip()
except:
  line = '\n'
  line += "Usage:"
  line += '\n'
  line += '\n'
  line += "  user:$ python SplitAttRep.py filename"
  line += '\n'
  line += '\n'
  line += "NB: 'filename' should match "
  line += "filename-dens.cube and filename-grad.cube"
  line += '\n'
  print(line)
  exit(0)

#Print initial output
line = '\n'
line += "Separating attractive and repulsive NCI surfaces..."
line += '\n'
line += '\n'
line += "Input files:"
line += '\n'
line += "  "
line += cubename+"-dens.cube "
line += cubename+"-grad.cube "
line += cubename+".vmd"
line += '\n'
print(line)

#Safely open cube files
try:
  densfile = open((cubename+"-dens.cube"),"r")
except:
  line = ""
  line += "Error: File "+cubename+"-dens.cube "
  line += "not found..."
  line += '\n'
  print(line)
  exit(0)
try:
  gradfile = open((cubename+"-grad.cube"),"r")
except:
  line = ""
  line += "Error: File "+cubename+"-grad.cube "
  line += "not found..."
  line += '\n'
  print(line)
  exit(0)

#Safely open VMD file
try:
  vmdfile = open((cubename+".vmd"),"r")
except:
  line = ""
  line += "Error: File "+cubename+".vmd "
  line += "not found..."
  line += '\n'
  print(line)
  exit(0)

#Create new cube files
attgradfile = open((cubename+"-attract-grad.cube"),"w")
repgradfile = open((cubename+"-repulse-grad.cube"),"w")

#Create new VMD files
attvmdfile = open((cubename+"-attract.vmd"),"w")
repvmdfile = open((cubename+"-repulse.vmd"),"w")

#Check number of atoms and clear junk
densfile.readline()
densfile.readline()
Natoms = densfile.readline()
Natoms = Natoms.strip()
Natoms = Natoms.split()
Natoms = Natoms[0]
Natoms = int(Natoms)
for i in range(Natoms+3):
  #Skip atomic positions
  densfile.readline()

#Save cube header
for i in range(Natoms+6):
  #Copy cube data and atomic positions
  line = gradfile.readline()
  attgradfile.write(line)
  repgradfile.write(line)

#Process cube files
line = ""
line += "Writing new cube files..."
line += '\n'
print(line)
densline = densfile.readline()
gradline = gradfile.readline()
while (densline and gradline):
  #Ugly, but it loops over the remaining lines
  densline = densline.split()
  gradline = gradline.split()
  #Error check
  if (len(densline) != len(gradline)):
    line = '\n'
    line += "Error: "+cubename+"-dens.cube does not match "
    line += cubename+"-grad.cube"
    line += '\n'
    print(line)
    exit(0)
  #Write data to new grad.cube files
  for i in range(len(densline)):
    #Write attractive interactions
    if (float(densline[i]) < 0):
      #Keep point
      extraspace = " "
      Nchar = len(gradline[i])
      for j in range(12-Nchar):
        #Make sure the number is 12 characters long
        extraspace += " "
      attgradfile.write((extraspace+gradline[i]))
    else:
      #Remove point
      attgradfile.write("  0.10000E+03")
    #Write repulsive interactions
    if (float(densline[i]) > 0):
      #Keep point
      extraspace = " "
      Nchar = len(gradline[i])
      for j in range(12-Nchar):
        #Make sure the number is 12 characters long
        extraspace += " "
      repgradfile.write((extraspace+gradline[i]))
    else:
      #Remove point
      repgradfile.write("  0.10000E+03")
  #Add newline characters
  attgradfile.write('\n')
  repgradfile.write('\n')
  #Read new data for the EOF check
  densline = densfile.readline()
  gradline = gradfile.readline()

#Create new cubename.vmd files
line = ""
line += "Writing new VMD files..."
line += '\n'
print(line)
vmdinfo = vmdfile.readlines()
for vmdline in vmdinfo:
  line = vmdline #Copy line
  vmdline = vmdline.split()
  if (len(vmdline) > 3):
    if (vmdline[2] == (cubename+"-grad.cube")):
      #Write new attractive line
      vmdline[2] = cubename+"-attract-grad.cube"
      line = ""
      for i in range(len(vmdline)):
        line += vmdline[i]
        if (i != (len(vmdline)-1)):
          line += " "
        else:
          line += '\n'
      attvmdfile.write(line)
      #Write new repulsive line
      vmdline[2] = cubename+"-repulse-grad.cube"
      line = ""
      for i in range(len(vmdline)):
        line += vmdline[i]
        if (i != (len(vmdline)-1)):
          line += " "
        else:
          line += '\n'
      repvmdfile.write(line)
    else:
      #No modifications are required
      attvmdfile.write(line)
      repvmdfile.write(line)
  else:
    #No modifications are required
    attvmdfile.write(line)
    repvmdfile.write(line)

#Tell users about the files before exiting
line = ""
line += "Attractive surfaces:"
line += '\n'
line += "  "
line += cubename+"-dens.cube "
line += cubename+"-attract-grad.cube "
line += cubename+"-attract.vmd"
line += '\n'
line += '\n'
line += "Repulsive surfaces:"
line += '\n'
line += "  "
line += cubename+"-dens.cube "
line += cubename+"-repulse-grad.cube "
line += cubename+"-repulse.vmd"
line += '\n'
line += '\n'
line += "Done."
line += '\n'
print(line)

#Close files and exit
vmdfile.close()
densfile.close()
gradfile.close()
attvmdfile.close()
repvmdfile.close()
attgradfile.close()
repgradfile.close()
exit(0)
