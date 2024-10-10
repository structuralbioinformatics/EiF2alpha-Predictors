import pandas as pd 
from math import exp

file = open("UTRfiltered_ALL.txt", "r")
matrix = file.readlines()

output=pd.DataFrame(columns=["Gene","Score","ORFs","Atf4like","Length","%GC"])

counter=0
DFin=0

for i in matrix:
    
    row=i.split("\t")
    counter+=1
    previousAUG=-99
    gc=0
    orf0=0 
    n0transl="closed"
    orf1=0
    n1transl="closed"
    orf2=0
    n2transl="closed"
    closedAUG=0
    cand0=0
    cand1=0
    cand2=0
    stop0=[]
    stop1=[]
    stop2=[]
    Atf4like=0
       
    for n in range(0,(len(row[1])-1)):
        
        if row[1][n] == "G" or row[1][n] == "C":
            
            gc+=1
        
        if n < len(row[1])-4:
            
            if ((row[1][n-3] == "A" or row[1][n-3] == "G") or (row[1][n+3] == "G")) and ((row[1][n] == "A") and (row[1][n+1] == "T") and (row[1][n+2] == "G")):
                
                if n-previousAUG>=30:
                    orf=n%3
                    
                    if orf == 0:
                        
                        if n0transl == "closed":
                            n0transl="open"
                            orf0+=1
                            cand0=n
                            
                    elif orf == 1:
                        
                        if n1transl == "closed":
                            n1transl="open"
                            orf1+=1
                            cand1=n
                            
                    elif orf == 2:
                        
                        if n2transl == "closed":
                            n2transl="open"
                            orf2+=1
                            cand2=n
                            
                previousAUG=n           
            
            if ((row[1][n] == "T") and (row[1][n+1] == "A") and (row[1][n+2] == "A")) or ((row[1][n] == "T") and (row[1][n+1] == "A") and (row[1][n+2] == "G")) or ((row[1][n] == "T") and (row[1][n+1] == "G") and (row[1][n+2] == "A")):
                
                orf=n%3
                
                if orf == 0:
                    if orf0 > 0:
                        
                        orf0=0
                        n0transl="closed"
                        closedAUG+=1
                        stop0.append(n)
                        cand0=0
                        
                elif orf == 1:
                    
                    if orf1 > 0:
                        
                        orf1=0
                        n1transl="closed"
                        closedAUG+=1
                        stop1.append(n)
                        cand1=0
                        
                elif orf == 2:
                    
                    if orf2 > 0:
                        
                        orf2=0
                        n2transl="closed"
                        closedAUG+=1
                        stop2.append(n)
                        cand2=0
            
    openAUG=(orf0+orf1+orf2)
    ORFs=openAUG+closedAUG
    length=len(row[1])-1
    pergc=gc/(len(row[1])-1)*100
    
    candDEF=[]
    candDEF.append(cand0)
    candDEF.append(cand1)
    candDEF.append(cand2)
    stopDEF=stop0+stop1+stop2
    
    for i in candDEF:
        
        w=0
        
        for j in stopDEF:
            
            if i>j and w==0:
                
                Atf4like+=1
                w+=1

    coef_list=[-0.0684291477691681,-0.242290543856012,-0.0634125358003748,-0.00815035995448194,0.00187094356680248,-0.170471273131394,0.000301400719870572,0.00849572542167972,0.000850864786405846,0.0112296578307007,0.000096082448197154]
    prescore=coef_list[0]+coef_list[1]*ORFs+coef_list[2]*Atf4like+coef_list[3]*length+coef_list[4]*pergc+coef_list[5]*ORFs*Atf4like+coef_list[6]*ORFs*length+coef_list[7]*ORFs*pergc+coef_list[8]*Atf4like*length+coef_list[9]*Atf4like*pergc+coef_list[10]*length*pergc    
    score=(exp(prescore))/(1+exp(prescore))        
    
    if score>0.7 and ORFs>0:
        
        output.loc[DFin,"Gene"] = row[0]
        output.loc[DFin,"Length"] = len(row[1])-1
        output.loc[DFin,"%GC"] = gc/(len(row[1])-1)*100
        output.loc[DFin,"Score"] = score
        output.loc[DFin,"ORFs"] = ORFs
        output.loc[DFin,"Atf4like"] = Atf4like
        DFin+=1

    if counter%100==0:
        
        print("Working on it... Transcript #"+str(counter))
    
output.to_excel("positives.xls")
