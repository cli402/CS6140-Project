#!/usr/bin/python
##  CSE6140 HW1
##  This assignment requires installation of networkx package if you want to make use of available graph data structures or you can write your own!!
##  Please feel free to modify this code or write your own
import time
import sys
import os
import bisect
global weight_count
class RunExperiments:
    '''<executable>  <graph_file.gr> <change_file.extra> <output_file>.'''

    def read_graph(self, filename):
        G=Graph(filename)
        return G

    def computeMST(self,graph):
        ## containner to hold the info
        heap_queue=[]
        nodes=[]
        node_set=set([])
        edge_set=[]
        global weight_count
        weight_count=0

        ##initialize the set of the inidividual vertices
        nodes=graph.nodes()
        uni_find=graph.uni_found()

        ##initialize the priority queue
        heap_queue=graph.edges()
        count=0
        index=0
        # while the spanning tree is not complete
        #while len(graph.edge_dict.keys())+1< len(nodes) and count+1<len(heap_queue):
        while count+1<=len(heap_queue) :
            edge=heap_queue[count]
            #if the nodes are not in the same set
            if uni_find.find_root(edge[1]) !=  uni_find.find_root(edge[2]):
                uni_find.add(edge[1],edge[2])
                #swap for convinience
                if edge[1] > edge[2]:
                    edge[2],edge[1] = edge[1],edge[2]
                graph.edge_dict[(edge[1],edge[2])]=edge[0]
                weight_count+=int(edge[0])
                index=count
                edge_set.append(edge)
                node_set.add(edge[1])
                node_set.add(edge[2])
            count+=1

        graph.max_weight=heap_queue[index][0]
        graph.list=edge_set
        return weight_count

    def recomputeMST(self,u, v, weight,G):
        G.add_edge(weight,u,v)
        global weight_count
        if u > v:
            u,v = v,u
        #if the edge exists in the MST,just update the mst
        if (u,v) in G.edge_dict:
            if G.edge_dict[(u,v)]> weight:
                weight_count-=(G.edge_dict[(u,v)]-weight)
                G.edge_dict[(u,v)]=weight
            return weight_count
        #if the edge not in the current MST:
        #case 1,edge not in the MST but in the graph:
        elif u in G.node_set and v in G.node_set:
             #case 1.1 node not in one MST:
                 if G.uni_find.find_root(u) != G.uni_find.find_root(v):
                     weight_count+=weight
                     G.edge_dict[(u,v)]=weight
                     G.uni_find.add(u,v)
                 #case 1.2 edge in the same MST,find the cycle and get the maximum weight in that cycle and compare
                 else:
                     #else rerun the computeMST to get the weight
                     if  weight < G.max_weight:
                             G.rerun_unifind()
                             return self.computeMST(G)
                 return weight_count

             #case 2: if one or two node not in the graph,just add the edge and weight
        elif u in G.node_set or v in G.node_set:
                 weight_count+=weight
                 if u in G.node_set:
                     G.add_node(v)
                 else:
                     G.add_node(u)
                 G.uni_find.add(u,v)
                 G.edge_dict[(u,v)]=weight
                 return weight_count
            #case 3: if none of the node in the node set
        else:
         G.add_node(u)
         G.add_node(v)
         G.rerun_unifind()
         return self.computeMST(G)





    def main(self):
        '''<executable>  <graph_file.gr> <change_file.extra> <output_file>.'''
        num_args = len(sys.argv)
        
        if num_args < 4:
            print "error: not enough input arguments"
            exit(1)

        graph_file = sys.argv[1]
        change_file = sys.argv[2]
        output_file = sys.argv[3]

        #Construct graph
        G =self.read_graph(graph_file)

        start_MST = time.time() #time in seconds
        MSTweight =self.computeMST(G) #call MST function to return total weight of MST
        total_time = (time.time() - start_MST) * 1000 #to convert to milliseconds

        #Write initial MST weight and time to file
        output = open(output_file, 'w')
        output.write(str(MSTweight) + " " + str(total_time)+'\n')


        #Changes file
        change_file = os.path.dirname(os.path.dirname(__file__))+'/data/'+change_file
        with open(change_file, 'r') as changes:
            num_changes = changes.readline()

            for line in changes:
                #parse edge and weight
                edge_data = list(map(lambda x: int(x), line.split()))
                assert(len(edge_data) == 3)

                u,v,weight = edge_data[0], edge_data[1], edge_data[2]

                #call recomputeMST function 
                start_recompute = time.time()
                new_weight = self.recomputeMST(u, v, weight, G)
                total_recompute = (time.time() - start_recompute) * 1000 # to convert to milliseconds

                #write new weight and time to output file
                output.write(str(new_weight) + " " + str(total_recompute)+'\n')



class unioin_find():
    def __init__(self,list):
        self.mapping = {} # maps a member to the group's mapping
        self.group = {} # maps a group mapping to the group (which is a set)
        for node in list:
            self.group[node]=set([node])
            self.mapping[node]=node

    def find_root(self,node):
        if self.mapping[node] == node:
            return node
        else:
            return self.find_root(self.mapping[node])
    def add(self, a, b):
        mappinga = self.find_root(a)
        mappingb = self.find_root(b)
        if mappinga == mappingb: return # nothing to do
        groupa = self.group[mappinga]
        groupb = self.group[mappingb]
        if len(groupa) < len(groupb):
            a, mappinga, groupa, b, mappingb, groupb = b, mappingb, groupb, a, mappinga, groupa

        elif len(groupa) == len(groupb) and int(mappinga) > int(mappingb):
            a, mappinga, groupa, b, mappingb, groupb = b, mappingb, groupb, a, mappinga, groupa
        groupa |= groupb
        del self.group[mappingb]
        self.mapping[b]=a
        self.mapping[mappingb]=mappinga

    def mapping_index(self,node):
        return self.mapping[node]


class Graph(object):
    def __init__(self,filename):
        a= os.path.dirname(os.path.dirname(__file__))
        a=a+'/data/'+filename
        self.list=[]
        self.node_set=set([])
        self.edge_dict={}
        self.max_weight=0

        with open(a) as f:
             lines = f.readlines()
             print [int(e.strip()) for e in lines[0].split(' ')]
             for node in lines[1:]:
                 node=[int(e.strip()) for e in node.split(' ')]
                 node.insert(0,node.pop(2))
                 self.list.append(node)

                 self.node_set.add(node[1])
                 self.node_set.add(node[2])
                 # self.add_node(node[1])
                 # self.add_node(node[2])
       #initialize the priority queue
        #pq.heapify(self.list)
        self.list.sort(key=lambda x: x[0])
        self.uni_find=unioin_find(self.node_set)

    def add_edge(self,weight,node1,node2):
        node=[weight,node1,node2]
        bisect.insort(self.list,node)

    def add_node(self,node):
        self.node_set.add(node)
        self.uni_find.mapping[node] = node
        self.uni_find.group[node] = set([node])
    def nodes(self):
        return self.node_set
    def edges(self):
        return self.list
    def uni_found(self):
        return self.uni_find
    def rerun_unifind(self):
        self.uni_find= unioin_find(self.node_set)
        self.edge_dict={}


if __name__ == '__main__':
    # run the experiments
    runexp = RunExperiments()
    runexp.main()
