import sys


class GraphNode:
    def __init__(self, value, incoming, outcoming, is_branch, is_visited, count):
        self.value = value
        self.count = count
        self.is_branch = is_branch
        self.is_visited = is_visited
        self.incoming = incoming
        self.outcoming = outcoming


class DeBruijnGraph:
    def __init__(self, reads, k):
        self.reads = reads
        self.nodes = dict()


    # assuming reads list is populated, builds a de bruijn graph from the reads
    def build_graph(self):
        for i in range(len(self.reads)):
            prev_node = ""
            for j in range(len(self.reads[i])-k+1):
                cur_kmer = self.reads[i][j:j+k]
                self.add_node(cur_kmer)
                if prev_node != "":
                    self.add_in_edge(cur_kmer, prev_node)
                    self.add_out_edge(prev_node, cur_kmer)
                prev_node = cur_kmer 
    

    # graph printing function for debugging purposes
    def printGraph(self):
        print("Printing Graph:")
        for i in self.nodes:
            x = self.nodes[i]
            print("Data: " + x.value)
            print("Incoming: " + str(x.incoming))
            print("Outgoing: " + str(x.outcoming))
        print()


    # adds the given kmer to the graph, creats GraphNode item
    def add_node(self, kmer):
        if kmer in self.nodes:
            self.nodes[kmer].count = self.nodes[kmer].count + 1
        else:
            self.nodes[kmer] = GraphNode(kmer, [], [], False, False, 1)
    

    # adds an in edge from kmer to out_kmer
    def add_in_edge(self, kmer, in_kmer):
        if in_kmer not in self.nodes[kmer].incoming:
             self.nodes[kmer].incoming.append(in_kmer)
        if (len(self.nodes[kmer].incoming) > 1):
            self.nodes[kmer].is_branch = True
    

    # adds an out edge from kmer to out_kmer
    def add_out_edge(self, kmer, out_kmer):
        if out_kmer not in self.nodes[kmer].outcoming:
            self.nodes[kmer].outcoming.append(out_kmer)
        if (len(self.nodes[kmer].outcoming) > 1):
            self.nodes[kmer].is_branch = True
    

    # removes all branch nodes from the graph
    def remove_branch_nodes(self):
        branching_nodes = []
        for name in self.nodes:
            if self.nodes[name].is_branch:
                branching_nodes.append(name)
        
        for i in range(len(branching_nodes)):
            self.remove_node(branching_nodes[i])


    # removes the given kmer from the graph and all corresponding edges
    def remove_node(self, kmer):
        self.nodes.pop(kmer)
        for x in self.nodes:
            node = self.nodes[x]
            for inc in node.incoming:
                if inc == kmer:
                    node.incoming.remove(kmer)
            for outc in node.outcoming:
                if outc == kmer:
                    node.outcoming.remove(kmer)

    
    #filters out graph to return list of nodes that appear more than once
    def get_good_reads(self):
        good_reads = []
        for i in range(len(self.reads)):
            if self.check_good_read(self.reads[i]) == True:
                good_reads.append(self.reads[i])
        return good_reads

    # given an individual read, checks all kmers in the read to see if it has
    # a count of at least 2
    def check_good_read(self, read):
        for i in range(len(read)-k+1):
            if self.nodes[read[i:i+k]].count == 1:
                return False
        return True

    # removes branch nodes and returns a tuple
    # where the first item is a list of contigs, and the second is a list of 
    # contig lengths
    def get_contigs(self, min_len=100):
        self.remove_branch_nodes()
        source_nodes = self.get_source_nodes()
        contigs, contig_lengths = self.traverse_graph(source_nodes, min_len)
        return contigs, contig_lengths

        
    #returns a list of all the source nodes (nodes with an indegree of 0)
    def get_source_nodes(self):
        source_nodes = []
        for key in self.nodes:
            node = self.nodes[key]
            if len(node.incoming) == 0:
                source_nodes.append(key)
        return source_nodes

    # given a list of source_node and a minimum contig lengt, returns a tuple
    # where the first item is a list of contigs, and the second is a list of 
    # contig lengths
    def traverse_graph(self, source_nodes, min):
        contigs = []
        contig_lengths = []
        counter = 0
        for node in source_nodes:
            contig = self.traverse_iteratitve(node)
            if (len(contig) >= min):
                contigs.append(contig)
                contig_lengths.append(len(contig))
            else:
                counter += 1
                # used for counting the amount of contigs less than min value
                #print("chaff: " + str(counter))
        return contigs, contig_lengths
        
    # given a starting node, returns the contig after reaching a sink node or an
    # already visted node
    def traverse_iteratitve(self, starting_node):
        result = ""
        if (len(self.nodes[starting_node].outcoming) == 0): #base case
            return starting_node #entire string
        else:
            result = starting_node
            next_node = self.nodes[starting_node].outcoming[0]
            while (next_node != "" and (not self.nodes[next_node].is_visited)):
                result = result + next_node[-1]
                self.nodes[next_node].is_visited = True
                if (len(self.nodes[next_node].outcoming) == 0):
                    next_node = ""
                else:
                    next_node = self.nodes[next_node].outcoming[0]
            return result


#returns a list of sequence reads from a file given through input
def get_sequence_reads(input):
    reads = []
    for line in input:
        reads.append(line.strip())
    return reads
    

# writes all elements in a given list to a given filename, each element on
# on a different line
def print_to_file(list, filename):
    with open(filename, 'w') as file:
        for i in range(len(list)):
            file.write(str(list[i]))
            file.write("\n")


if __name__ == "__main__":
    reads = get_sequence_reads(sys.stdin) 
    k = 31
    graph = DeBruijnGraph(reads, k)
    print("******************** BUILDING ORIGINAL GRAPH ********************")
    graph.build_graph()
    #graph.printGraph()

    print("******************** GETTING GOOD READS ********************")
    good_reads = graph.get_good_reads()

    print("******************** PRINTING GOOD READS ********************")
    print_to_file(good_reads, "good_reads")

    print("******************** REBUILDING GRAPH USING GOOD READS ********************")
    new_graph = DeBruijnGraph(good_reads, k)
    new_graph.build_graph()
    #new_graph.printGraph()
    
    print("******************** GETTING CONTIGS ********************")
    contigs, contig_lengths = new_graph.get_contigs()
    contig_lengths.sort()

    print("******************** PRINTING CONTIGS ********************")
    print_to_file(contigs, "output_contigs")
    print_to_file(contig_lengths, "contig_lengths")
    
    # Used for determining N50 value
    #print(sum(contig_lengths))




    
    
    