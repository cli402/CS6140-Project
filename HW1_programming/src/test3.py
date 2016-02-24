from run_experiments import *
G=Graph('rmat0406.gr')
print G.edges()
print G.edge_dict
print G.max_weight
print G.uni_find.mapping
print G.uni_find.group
def computeMST(graph):
        ## containner to hold the info
        heap_queue=[]
        nodes=[]
        node_set=set([])
        edge_set={}

        global weight_count
        weight_count=0

        ##initialize the set of the inidividual vertices
        nodes=graph.nodes()
        uni_find=graph.uni_found()

        ##initialize the priority queue
        heap_queue=graph.edges()
        count=0
        # while the spanning tree is not complete
        while len(node_set) < len(nodes):
            edge=heap_queue[count]
            count+=1
            #if the nodes are not in the same set
            if uni_find.find_root(edge[1]) !=  uni_find.find_root(edge[2]):
                uni_find.add(edge[1],edge[2])
                #swap for convinience
                if edge[1] > edge[2]:
                    edge[2],edge[1] = edge[1],edge[2]
                graph.edge_dict[(edge[1],edge[2])]=edge[0]
                weight_count+=int(edge[0])
                graph.max_weight=max(edge[0],graph.max_weight)
                if edge[1] not in node_set:
                    node_set.add(edge[1])
                if edge[2] not in node_set:
                    node_set.add(edge[2])
        print G.uni_find.mapping
        print G.uni_find.group
        return weight_count
a=computeMST(G)
print G.uni_find.find_root(12)
print G.uni_find.find_root(5)
print a