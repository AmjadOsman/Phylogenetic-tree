import numpy as np
import math
from heapq import nsmallest
from ete3 import Tree



def compare(string1, string2, no_match_c=' ', match_c='|'): #gets the  of difference between sequence
    if len(string2) < len(string1):
        string1, string2 = string2, string1
    result = ''
    n_diff = 0
    for c1, c2 in zip(string1, string2):
        if c1 == c2:
            result += match_c
        else:
            result += no_match_c
            n_diff += 1
    delta = len(string2) - len(string1)
    result += delta * no_match_c
    n_diff += delta
    return  n_diff


def matrix(x,y): # makes a matrix 
    return [[0]*y for i in range(x)]


def printMatrix(matrix):
	#function to print the matrix 
	for i in range(len(matrix)):
		print(matrix[i])
	print()


x={}
with open('msa.fasta','r') as file:# extracts sequences and saves them in a dictionary 
    for lines in file:
        line=lines.strip()
        if line.startswith('>'):
            seqname=line[1:]
            if seqname not in x: 
                x[seqname]=[]
            continue 
        sequence=line
        x[seqname].append(sequence)
        fasta={k: ''.join(v) for k, v in x.items()}
    
    
# print(fasta
# 


match_score=1
mismatch_score=-1
gap=-2


def np_global_alignment(x,y): #global alignment acting as a MSA 


    n_matrix=np.zeros((len(x)+1,len(y)+1)) # add +1 here 
    main_matrix= np.zeros((len(x)+1,len(y)+1))
    trace_back=np.zeros((len(x)+1,len(y)+1),dtype=str)
    # trace_back[0][0]=(0,0)

    for i in range(len(x)): # match mismatch matrix 
        for j in range(len(y)):
            if x[i]==y[j]:
                n_matrix[i][j]=match_score
 
            else:
                n_matrix[i][j]=mismatch_score


    for i in range(len(x)+1):
        main_matrix[i][0]=gap*i
    for j in range(len(y)+1):
        main_matrix[0][j]=gap*j

    for i in range(1,len(x)+1):
        for j in range(1,len(y)+1):
            left=main_matrix[i][j-1]+gap
            up=main_matrix[i-1][j]+gap
            diagonal= main_matrix[i-1][j-1]+n_matrix[i-1][j-1]
            main_matrix[i][j]=max(left,diagonal,up)


    align1 =[]
    align2=[]
    xi=len(x)    
    yj=len(y)  
    while(xi > 0 and yj > 0):
        if (xi>0 and main_matrix[xi][yj] == main_matrix[xi-1][yj-1]+n_matrix[xi-1][yj-1]):
			# Diag 
            align1.append(x[xi-1])
            align2.append(y[yj-1])
            xi=xi-1
            yj=yj-1
        elif (yj>0 and main_matrix[xi][yj] == main_matrix[xi][yj-1]+gap):
			# Left 
            align2.append(y[yj-1])
            align1.append('-')
            yj=yj-1

        else:
			# Up 
            align1.append(x[xi-1])
            align2.append('-')
            xi = xi-1
		

    alignment1=align1[::-1]
    alignment2=align2[::-1]
    return alignment1,alignment2


distancedict={}


def p_distancematrix(function1,function2):# calulates the p distance 

    for key,value in fasta.items():
        for key1 ,value1 in fasta.items():
            if key+'-'+key1 not in distancedict and value != value1:
                align1,align2=function1(value,value1) #function1 = global alignment =aligns sequences 
                string1=''.join(align1)
                string2=''.join(align2)
                length1=len(string1)
                length2=len(string2)
                num_difference=function2(string1,string2)# function2 = compare: gets the # of difference between sequence
                if length1==length2:
                    pdistance= num_difference/length1
                else:
                    print('length arent equal')

                prespecies= key+'-'+key1
                distancedict[prespecies]=[]
                species=prespecies    
                distancescore =pdistance
                distancedict[species].append(distancescore) 
    print(distancedict)
    return distancedict
   




              
matrix=p_distancematrix(np_global_alignment,compare)
model_dict={}
filtered_modeldict= {}

def juke_and_cantor_model(): # nucleotide subsitiution model- formula
    for key,value  in distancedict.items():
            if key not in model_dict:
                pvalue=float(''.join(str(i)for i in value ))
                jcvalue=-3/4*math.log(1-4/3*pvalue) 
                # model_dict[key]=[]

                if jcvalue not in model_dict.values():
                   
                    model_dict[key]=jcvalue
                else:
                    model_dict[key]=None

                    continue 
                 
    print(model_dict)
    return model_dict

model=juke_and_cantor_model()



tup_dict={}

def keytuple(): # turns dictionary key to tuple 
    filtered_modeldict= {k:v for k,v in model_dict.items()  if v != None} #removes duplicates
    for k,v in filtered_modeldict.items(): # seperates the key into tuple 
        key1= tuple(k.split('-'))
        if key1 not in tup_dict:
            tup_dict[key1]=v
    print(tup_dict)

tuplfinal=keytuple()


def Pre_newick_format(): # switching dictionary to pre processed form 
    global tree
    tree =[]
    for k,v in tup_dict.items():
        for key in k:
            smallestvalue= min(tup_dict.values())
            if v==smallestvalue:
                first_tuple=tuple(k)
                large_tuple=(first_tuple,)
    
                for i in large_tuple:
                    for key1,value1 in tup_dict.items():
                        for key2,value2 in tup_dict.items():

                            if (smallestvalue!=value1) and  (key1[0]==key2[0] ) and (key2[1]!= key1[1] )and (key1[0]or key1[1]or key2[0]or key1[1]in i or large_tuple):
                                secondvalue= min(value1,value2)
                                if value2 == secondvalue:
                                    s=tuple(key2)
                                else: 
                                    s= tuple(key1)
    
                                if (value1 or value2==secondvalue )and (smallestvalue!=secondvalue) and (s[0] not in large_tuple):
                                        first_specie=(s[0],)
                                        tuple_1=(large_tuple+first_specie)
                                        largest_tuple=(tuple_1,)
                                        
                                else:
                                        first_specie=(s[1],)
                                        tuple_1=(large_tuple+first_specie)
                                        largest_tuple=(tuple_1,) 

                                secondsmallest_value=nsmallest(2,tup_dict.values())[-1]
        
                                for key4,value4 in tup_dict.items():
                                        for key5,value5 in tup_dict.items():
                                                for i in largest_tuple:
                                                    for j in i:

                                                            if (key5 and  key4 and (key4[0] or  key4[1] or key5[0] or  key5[1]) not in largest_tuple , j, i) and  (value4 and value5 != secondsmallest_value) and (value5 and value4 != secondvalue )and (value5 != value1) and (value5 !=value2) and (value4 != value2 ) and (value4 != value1) and (value4!=value5):
                                                                savevalue1=value5
                                                                savevalue2=value4
                                            
                                                                if (key4[0] or key4[1] and  key5[0] or key5[1] not in tree) and (value5!= value4) and (value5 and value4 != savevalue1 and savevalue2):
                                                                    finalvalue =min(savevalue1,savevalue2,value5,value4)
                                                                    if value4 ==finalvalue and key4[1] in tree:
                                                                        tree.append((key4[0],)+largest_tuple)
                                                                    else:
                                                                        tree.append ((key4[1],)+largest_tuple)
                                                            
                                                                for key3,value3 in tup_dict.items():

                                                            
                                                                    if value3==secondsmallest_value:
                                                                        tree.append(tuple(key3))
                                                                        print (tree)
                                                                        return tree
                                                        

new=Pre_newick_format()

def final_newick_format(): # finalizes the newick format
    y= str(tree).strip('[]')
    newy= y.replace("'","")
    final_tree='(' + newy + ');'
    # final_tree=
    t= Tree(final_tree,format=9)
    t.render('phylogenetic_Tree.png')

    print(t)
    

news= final_newick_format()
