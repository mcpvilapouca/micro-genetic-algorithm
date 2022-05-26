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
var3="BETAM1=" + str(data[3])
var4="TAUM1=" + str(data[4])
var5="BETAM2=" + str(data[5])
var6="TAUM2=" + str(data[6])
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
with open("auxfile2.inp","w") as f:
    for line in newline:
        f.writelines(line)
# #var2
with open("auxfile2.inp","r") as f:
    newline=[]
    for word in f.readlines():
        newline.append(word.replace("BETAM1=1.0",var3))  ## Replace the keyword while you copy.
with open("auxfile3.inp","w") as f:
    for line in newline:
        f.writelines(line)
# #var2
with open("auxfile3.inp","r") as f:
    newline=[]
    for word in f.readlines():
        newline.append(word.replace("TAUM1=1.0",var4))  ## Replace the keyword while you copy.
with open("auxfile4.inp","w") as f:
    for line in newline:
        f.writelines(line)
# #var2
with open("auxfile4.inp","r") as f:
    newline=[]
    for word in f.readlines():
        newline.append(word.replace("BETAM2=1.0",var5))  ## Replace the keyword while you copy.
with open("auxfile5.inp","w") as f:
    for line in newline:
        f.writelines(line)
# #var2
with open("auxfile5.inp","r") as f:
    newline=[]
    for word in f.readlines():
        newline.append(word.replace("TAUM2=1.0",var6))  ## Replace the keyword while you copy.
with open("parameters.inp","w") as f:
    for line in newline:
        f.writelines(line)
# Remove auxiliar files
os.remove("auxfile.inp")
os.remove("auxfile1.inp")
os.remove("auxfile2.inp")
os.remove("auxfile3.inp")
os.remove("auxfile4.inp")
os.remove("auxfile5.inp")
#Remove blank lines from the file
with open("parameters.inp","r") as f:
        lines=f.readlines()
with open("parameters.inp","w") as f:
    [f.write(line) for line in lines if line.strip() ]
