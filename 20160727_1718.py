import os
os.chdir("C:/Users/mzy/Desktop")
file1=open('qwe.csv','r')
file2=open("001","a")
list1=[]
list2=[]
l=[]
x=[]
for i,line in enumerate(file1):
    if i==0:
        fields=line.strip("\n").split(",")
        for j in range(len(fields)):
            l.append(fields[j])
    if i > 0:
        fields=line.strip("\n").split(",")
        for j in range(i+1,2328):
             x.append(l[i])
             x.append(l[j])
             x.append(fields[j])
             list2.append(x)
             x=[]
             file2.write(l[i]+'\t'+l[j]+'\t'+fields[j]+'\n')
file2.close()


