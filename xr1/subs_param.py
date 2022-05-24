#Replace parameters in inp files
import io
import os
data = []
#Read variables from file par.txt
with io.open('par.txt', encoding='latin-1') as myfile:
    for i in myfile.readlines():
        data.append(i)
#Define each variable as a string
var="C10=" + str(data[0])
var1="K11=" + str(data[1])
var2="K12=" + str(data[2])
#Replace the variables and create a new files
#var
with open("param_orig.inp","r") as f:
    newline=[]
    for word in f.readlines():
        newline.append(word.replace("C10=1.0",var))  ## Replace the keyword while you copy.
with open("auxfile.inp","w") as f:
    for line in newline:
        f.writelines(line)
#var1
with open("auxfile.inp","r") as f:
    newline=[]
    for word in f.readlines():
        newline.append(word.replace("K11=1.0",var1))  ## Replace the keyword while you copy.
with open("auxfile1.inp","w") as f:
    for line in newline:
        f.writelines(line)
# #var2
with open("auxfile1.inp","r") as f:
    newline=[]
    for word in f.readlines():
        newline.append(word.replace("K12=1.0",var2))  ## Replace the keyword while you copy.
with open("parameters.inp","w") as f:
    for line in newline:
        f.writelines(line)
# Remove auxiliar files
os.remove("auxfile.inp")
os.remove("auxfile1.inp")
#Remove blank lines from the file
with open("parameters.inp","r") as f:
        lines=f.readlines()
with open("parameters.inp","w") as f:
    [f.write(line) for line in lines if line.strip() ]
