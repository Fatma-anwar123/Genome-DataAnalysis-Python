# -*- coding: utf-8 -*-
"""
Created on Sun Apr  3 22:42:57 2022

@author: DELL
"""
import networkx, matplotlib.pyplot as plot 

def K_mer(seq,k): #this function take one read(string), length of k-mer
    kmer=[]
    nkemr=len(seq)-k+1  #from eqaution calc num. of k-mers can we have from this read
    for i in range(0,nkemr,1):
        kmer.append(seq[i:i+k])  #append to kmer list all substrings
    return kmer
  


def DebruijnGraph(read,k):  #take list of reads and k
    kmers=[]
    graph={}
    f=[]
    prefix=[]
    suffix=[]
    Edge=[]
    res=[]      
    for i in range(len(read)):# for each read return kmers of it and append in Kmers list
        kmers+=K_mer(read[i], k)
    for s in range (len(kmers)): #for each k-mer return kmers of k-1
         f+=K_mer(kmers[s], k-1)
    for item in range(0,len(f),2):
        prefix.append(f[item]) # appened prefix
        suffix.append(f[item+1]) # appened suffix
    for i in f:  # for remove duplicate in Nodes
         if i not in res:
           res.append(i)
    graph["nodes"]=res  #add it to key nodes in graph dictionary
    for j in range (len(prefix)):
         Edge.append([prefix[j],suffix[j]]) 
    graph['edges']=Edge # append edges to key edge in graph dictionary
    return graph
        
def visualizeDBGraph(graph): # to visulaize  graph
    dbGraph = networkx.DiGraph()
    dbGraph.add_nodes_from(graph['nodes']) #Add the nodes to the graph 
    dbGraph.add_edges_from(graph['edges']) #Add the edges to the graph
    networkx.draw(dbGraph, with_labels=True, node_size=1000) 
    plot.show()

graph = DebruijnGraph(['TTACGTT','CCGTTA','GTTAC','GTTCGA','CGTTC'],5)
visualizeDBGraph(graph) 


       
  
    
  
#print(K_mer("TTACGTT",5))   
    

    
#print(DebruijnGraph(['TTACGTT','CCGTTA','GTTAC','GTTCGA','CGTTC'],5))
        
        