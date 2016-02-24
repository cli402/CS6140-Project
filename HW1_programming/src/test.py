__author__ = 'Amras'
import os
import copy
import heapq as pq

class Graph(object):
    def __init__(self,filename):
        a= os.path.dirname(os.path.dirname(__file__))
        a=a+'/data/'+filename
        self.list=[]
        self.node_set=set([])
        with open(a) as f:
             lines = f.readlines()
             print [int(e.strip()) for e in lines[0].split(' ')]
             for node in lines[1:]:
                 node=[int(e.strip()) for e in node.split(' ')]
                 node.insert(0,node.pop(2))
                 self.list.append(node)
                 if node[1] not in self.node_set:
                     self.node_set.add(node[1])
                 if node[2] not in self.node_set:
                     self.node_set.add(node[2])
       #initialize the priority queue
        pq.heapify(self.list)
        print self.list


    def add_edge(self,weight,node1,node2):
        node=[weight,node1,node2]
        pq.heappush(list,node)
        if node1 not in self.node_set:
            self.node_set.add(node1)
        if node2 not in self.node_set:
            self.node_set.add(node2)
    def n_largest(self,n):
        data = pq.nlargest(n, enumerate(self.list), key=lambda x:x[1])
    def nodes(self):
        return self.node_set
    def edges(self):
        return copy.deepcopy(self.list)

G=Graph('rmat0406.gr')



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



c=G.nodes()
e=G.edges()
while e:
    a=pq.heappop(e)
    print a
unifi=unioin_find(c)
print e
b=G.edges()
print b

# for ele in e:
#     unifi.add(ele[1],ele[2])

print unifi.group
print unifi.mapping


def computeMST(graph):
        ## containner to hold the info
        heap_queue=[]
        nodes=[]
        node_set=set([])
        weight_count=0

        ##initialize the set of the inidividual vertices
        nodes=graph.nodes()
        uni_find=unioin_find(nodes)

        ##initialize the priority queue
        heap_queue=graph.edges()

        # while the spanning tree is not complete
        while len(node_set) < len(nodes):
            edge=pq.heappop(heap_queue)
            #if the nodes are not in the same set
            if uni_find.find_root(edge[1]) !=  uni_find.find_root(edge[2]):
                uni_find.add(edge[1],edge[2])
                weight_count+=int(edge[0])
                if edge[1] not in node_set:
                    node_set.add(edge[1])
                if edge[2] not in node_set:
                    node_set.add(edge[2])
        return weight_count

import time
start_recompute = time.time()
a=computeMST(G)
total_recompute = (time.time() - start_recompute) * 1000
print total_recompute





    # def find_path(self,u,v,G):
    #     parent_u=copy.deepcopy(u)
    #     parent_v=copy.deepcopy(v)
    #     path_u = [u]
    #     path_v = [v]
    #     result=[]
    #     root=G.find_root[u]
    #     state= False
    #     while G.mapping[parent_u] != root and G.mapping[parent_v]!=root :
    #         if G.mapping[parent_u] == v:
    #             result= path_u.append(v)
    #             state = True
    #             break
    #         if G.mapping[parent_v] == u:
    #             result = path_v.append(u)
    #             state = True
    #             break
    #         path_u = path_u.append(G.mapping[parent_u])
    #         path_v = path_v.append(G.mapping[parent_v])
    #
    #         parent_v=G.mapping[parent_v]
    #         parent_u=G.mapping[parent_u]
    #
    #     if not state:
    #         state=True
    #         while G.mapping[parent_v] != root:
    #             if G.mapping[parent_v] == u:
    #                 result = path_v.append(u)
    #                 state =False
    #                 break
    #             path_v = path_v.append(G.mapping[parent_v])
    #             parent_v=G.mapping[parent_v]
    #
    #     if not state:
    #         state=True
    #         while G.mapping[parent_u] != root:
    #             if G.mapping[parent_u] == v:
    #                 result = path_u.append(v)
    #                 state =False
    #                 break
    #             path_u = path_u.append(G.mapping[parent_v])
    #             parent_v=G.mapping[parent_v]
    #     if not state:
    #         path_u.reverse()
    #
    #         result = path_u[1:]+path_v
    #     return result