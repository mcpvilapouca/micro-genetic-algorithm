print("start_rpt...")
#Run abaqus
import subprocess
import os
pathh=os.getcwd()
#
#eliminate the report if the file is already there
if os.path.isfile("report.rpt"):
    os.remove("report.rpt")

#Create the stretch-stress data from abaqus
ger = subprocess.Popen('abaqus viewer noGui=get_rpt.py ',
stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
ger.wait()

#Eliminate the first three lines that have text
#be aware that depending on the model the number of variables
#to eliminate might change
print("open")
lines = open('report.rpt').readlines()
open('report.rpt', 'w').writelines(lines[3:])
print("...end_rpt")
