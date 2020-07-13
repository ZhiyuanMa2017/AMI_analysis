import os
import math
import numpy
from scipy import stats 
os.chdir("C:/Users/mzy/Desktop")

file1=open('finalgene.csv','r')
list1=[]
for i,line in enumerate(file1):
    if i>0:
        fields=line.strip("\n").split(",")
        list1.append(fields[1:])
list2=zip(*list1)
a,b=stats.spearmanr(list2)
numpy.savetxt('pvalue.txt',b)