# ==========
# Name: Teh Jia Xuan
# Student ID: 32844700
# ==========
from collections import deque 

# ==========
# Q1
"""
    To find substring of a string
    given that starting and ending of a string
"""
class OrfFinder:
    
    def __init__(self, word):
        """ 
        a magic method to initialise attribute and create 
        prefix of a suffix in the trie

        precondition: a string is given
        postcondition: prefix of a suffix has been created in the trie

        Input:
            word: a string

        Return:
            None

        Time complexity:
            Best case analysis: O(N^2) where N is the number of character in the word as inserting 
                                N prefix for N times so total it would be N^2
            Worst case analysis: same as best case

        Space complexity:
            Input space analysis: O(N) where N is the length of word
            Aux space analysis: O(N^2) as inserting all the prefix into trie it will be N^2 number of node
                                in the trie and each node has a list size of 5 to link them all together so it would be 5N^2
                                = O(N^2)
        """
        #time: O(1)
        #space: O(1)
        self.trie = Trie()
        current = self.trie.root
        #put the empty string into the trie before we start the rest of substring 
        self.trie.link_current(current, 0, True, 0)
        self.prefix = ""

        #time: O(N^2) where N is the number of character in the word as inserting N prefix and N times
        #space: O(N^2) where N is the number of character in the word as there will be N^2 nodes in the trie and for each node
        # we are creating a list to link them
        for i in range(len(word) - 1, -1, -1):
            #adding each char to prefix every loop
            self.prefix = word[i] + self.prefix
            #insert each prefix to the trie
            #time: O(M) where M is the length of prefix as inserting char into the trie and traverse to nx char so is O(M)
            #space: O(M) where M is the length of prefix as creating M number of list in trie
            self.trie.insert(self.prefix, i)


    def find(self, starting, ending):
        """ 
        find the substring within starting and ending
        return all of the substring

        precondition: a valid starting and ending string 
        postcondition: return all of the substring within starting and ending of the string

        Input:
            starting: a string representing starting of the string
            ending: a string representing ending of the string

        Return:
            None
        
        Approach description:
        Intially i have created a suffix, which is inserting all the prefix into a trie and while inserting them
        i assigned the position that the char appear in the word into the node attribute called position. For example lets say 
        for a word aabbbc i have inserted bc so in my node of b i have position of 5 now in the position list. Then i insert
        bbc. The position list of first b will be [5,4] and when i insert bbbc it would be [5,4,3]. So i save the position
        of character. After inserting all the prefix. If i want to find substring between ab and bc for starting i loop through
        a and get the position of the maximum index lets say ab i get position b. For ending i loop through ending to get
        the position of the maximum index lets say bc i get c. After getting all the index, i check if the starting and ending
        overlap by checking their index, if is not overlap then enter the slicing loop and get all the substring between the index
        of starting and ending. If it overlap lets say starting index > ending index, it doesnt enter the loop and return an empty list.

        Time complexity:
            Best case analysis: O(1) where the program cant even find the first word of starting or ending in the trie
                                while looping through starting and ending and they will just return empty list

            Worst case analysis: O(T + U + V) where T is the length of starting string, U is the length of ending string
                                 V is the number of character in the return list. 
                                 Looping through ending to get the position is O(U) where U is the length of ending string
                                 Once we got the ending index, we need to loop through the starting string to get the position,
                                 so is O(T) where T is the length of starting String.
                                 Once we got both we loop through both of them and slice the substring accordingly so is O(X*Y*N) = O(V)
                                 where X is the number of starting index in the list and Y is the number of ending index in the list
                                 and V is the number of character in the result list. So overall is O(T + U + V). 
                                 
                                 But there is a scenario where it could be O(T + U) is when the starting and the ending index is overlapping 
                                 so it is impossible to generate substring that is overlap. So if it doesnt meet the condition it would be O(T + U)
                                 so it wont go through the looping and generate result. But worst case is still O(T + U + V)
        Space complexity:
            Input space analysis: O(T + U) where T is length of starting and U is length of ending string
            Aux space analysis: O(V) where V is the number of character in the result list. as appending V number of character it is O(V)
        """
        current = self.trie.root
        result = []
        ending_index = []
        starting_index = []
    
        #time O(U) where U is length of ending
        #space O(1) in place as it is just traversing and getting position
        for i in ending:
            #getting ascii index
            index = ord(i) - 65 + 1
            # if ending is not there -> word doesnt exist return []
            if (current.links[index] is None):
                return []
            #getting position and link
            ending_index = current.links[index].position
            current = current.links[index]

        current = self.trie.root
        #time O(T) where T is the length of starting
        #space O(1) in place, as it just traversing and getting index
        for i in starting:
            index = ord(i) - 65 + 1
            #if no starting -> doesnt exist return []
            if (current.links[index] is None):
                return []
            #getting index
            starting_index = current.links[index].position
            current = current.links[index]

        #if it overlap then dont go in this loop so it would be O(1) if it overlap
        if (ending_index[0] - len(ending) >= starting_index[len(starting_index) - 1]):
            #time:  O(V) where V is the number of character in result list
            #Space: O(V) where V is the number of character in the result list
            #for every starting index and ending index slice them accordingly to get all the substring
            #within starting and ending
            for i in starting_index:
                for j in ending_index:
                    #if it overlap then dont slice it
                    if(j - len(ending) >= i):
                        #O(N) where N is the length of word
                        result.append(self.prefix[i-len(starting) : j])
     
        return result
            

class Node:
    """ 
        magic method that initialise all the attribute 

        precondition: a valid size with integer is given, data is optional
        postcondition: necessary attributes are created

        Input:
            data(optional): the data to put in the node
            size: size of the range of alphabet which is A-D and terminal so is 5 

        Return:
            None

        Time complexity:
            Best case analysis: O(1) same as worst case
            Worst case analysis: O(1) as creating list of size "size" is constant

        Space complexity:
            Input space analysis: O(1) as size is constant 
            Aux: O(1) as creating attributes
    """
    def __init__(self, data = None, size = 5):
        self.data = data
        self.links = [None] * size
        self.end = False
        self.position = []
        self.index = 0

"""
    to insert string as a trie 
    have all the necessary method of trie
"""
class Trie:
    
    def __init__(self):
        """ 
        magic method that initialise the root attribute

        precondition: None
        postcondition: the root is created

        Input:
            None

        Return:
            None

        Time complexity:
            Best case analysis: O(1) as creating a node
            Worst case analysis: O(1) same as best case

        Space complexity:
            Input space analysis: O(1)  
            Aux: O(1) in place as node has a constant size of 5
        """
        self.root = Node()

    
    def insert(self, key, word_position):
        """ 
        insert string into trie

        precondition: a valid string is needed to insert, and the position of prefix of a suffix in the string 
        postcondition: the string is inserted into the trie

        Input:
            key: a string representing the substring string to insert to the trie
            word_position: the position of the substring in the string 

        Return:
            None

        Time complexity:
            Best case analysis: O(N) where N is the length of key, the rest is O(1) as those are assignment and link_current is O(1)
            Worst case analysis: O(N) where N is the length of key, same as best case

        Space complexity:
            Input space analysis: O(1) as key is a string and word position is an integer
            Aux: O(N) where N is the length of key as creating N lists for the character in key
    """
        current = self.root
        #time: O(N) where N is the length of key
        #space: O(N) where N is the length of key as creating N lists for the character in key
        for i in range(len(key)):
            #getting the index in the node list
            #position is the position of char in the word
            index = ord(key[i]) - 65 + 1
            position = word_position + i + 1

            #link it to the next node of index O(1) time and space aux
            current = self.link_current(current, index, False, position)
            current.index = position  
        #terminate after loop to the end of the key
        #linking the node list to 0 -> indicating terminates
        index = 0
        self.link_current(current, index, True, 0)
    

    
    def link_current(self, current, index, end_indicator, word_position):
        """ 
        link current node to the next char node using list

        precondition: a current representing current node, index representing which char index in the list
                      a boolean representing the node is the end of the substring
                      word position is the position of the char in the string

        postcondition: the node is linked to the next char node and the current is next char node

        Input:
            current: a node representing current node
            index: integer representing char index in the list
            end_indicator: a boolean indicates this is the end of the substring
            word_position: an integer representing the position of the char in the string

        Return:
            current: a node representing current node

        Time complexity:
            Best case analysis: O(1) as everything is indexing and appending word_position
            Worst case analysis: O(1) same as best case

        Space complexity:
            Input space analysis: current is a node so O(1), index is an integer O(1), end indicator is boolean and word
                                  position is an integer so is O(1). Overall is O(1) space 
            Aux: O(1) as everything is indexing and appending word position
    """
        #time:O(1) as indexing and appending word position
        #space:O(1) as indexing and appending word position
        #if there already existing word traverse to that node and append the word position
        if current.links[index] is not None:
            current = current.links[index]
            current.position.append(word_position)
        else:

            #time:O(1) as indexing and appending word position
            #space:O(1) as indexing and appending word position
            #if is none then create a new node in that index and append position
            #if end indicator is true then it is the end of the substring 
            current.links[index] = Node()
            current = current.links[index]
            current.position.append(word_position)
            current.end = end_indicator

        return current

# ==========
# Q2

class Graph:
    """
    To create a graph with vertex and edges
    """
    def __init__(self, preferences, officers_per_org, min_shift, max_shift, day):
        """ 
        magic function to initialise variable

        precondition: a valid list of list preferences of officer shift, a valid list of list shift needed by a 
                      company, an integer representing min shift, an integer representing max shift and an integer 
                      representing day

        postcondition: a graph object representing the graph of flow network is created

        Input:
            preferences: list of list representing officer preferred shift
            officers_per_org: list of list representing shift that needed by a company
            min_shift: an integer representing minimum shift
            max_shift: an integer representing maximum shift
            day: an integer representing day

        Return:
            None

        Time complexity:
            Best case analysis: O(2N) = O(N) where N is the total number of node in the graph which is
                                total_officer + total_shift + total day + total company + sink + source + min max node
                                and loop through all the node to create vertex so is O(2N) = O(N)
            Worst case analysis: same as best case

        Space complexity:
            Input space analysis: O(M + U) where M is the number of preference list in preferences and U
                                  is the number of company in officers per org

            Aux space analysis: O(N) where N is the total number of node
        """
    
        #creating index to access vertex
        total_day = day
        self.source_index = 0
        self.min_max_node = 1
        self.officer_count = len(preferences) + 1
        self.days_count = self.officer_count + total_day #CHANGE DAY
        self.shift_count = self.days_count + (total_day * 3) #as each day has 3 shifts
        self.companies_count = self.shift_count + len(officers_per_org)
        self.sink_index = self.companies_count + 1
        self.vertices = [None] * (self.sink_index + 1) #space O(N) where N is the number of node
        #creating vertex
        #time O(N) where creating N vertex into edge list where N is the number of node
        #space O(1)
        for i in range(len(self.vertices)):
            self.vertices[i] = Vertex(i)
        
        
    
class Vertex:
    """
    To create a vertex with this class
    """
    def __init__(self, id) -> None:
        """ 
        magic function to initialise value of vertex

        precondition: id is an integer
        postcondition: a vertex is created with the id

        Input:
            id representing id of vertex

        Return:
            None

        Time complexity:
            Best case analysis: O(1) as it is initialisation
            Worst case analysis: O(1) as it is initialisation

        Space complexity:
            Input space analysis: O(1)as it is initialisation
            Aux space analysis: O(1) as it is initialisation

        """
        self.id = str(id)
        self.edges = []
        self.demand = 0
        self.shift_needed = [None]
        self.visited = False
        self.discovered = False
        self.parent = None
        




    def __str__(self):
        """ 
        magic function for vertex to print

        precondition: a vertex object is created
        postcondition: print a string representing vertex 

        Input:
            None
        Return:
            return_string: a string represent vertex

        Time complexity:
            Best case analysis: O(1) as it is just initialisation
            Worst case analysis: O(1) as it is just initialisation

        Space complexity:
            Input space analysis: O(1)as it is just initialisation
            Aux space analysis: O(1)as it is just initialisation
        """
        return_string = str(self.id)
        return return_string

class Edge:
    """
    class to create edge
    """
    def __init__(self, u, v, lower_bound = None, capacity = None):
        """ 
        magic function to initialise edge

        precondition: u and v is a valid integer representing vertex, lower bound representing the lower bound
                      of the edge, capacity represent capacity of the edge
        postcondition: an edge object is created to link the vertex

        Input:
            u: vertex representing vertex starting 
            v: vertex representing vertex that ending
            lower_bound: an integer representing lower bound of the edge
            capacity: an integer representing capacity of the ddge
        Return:
            None

        Time complexity:
            Best case analysis: O(1) as it is just initialisation
            Worst case analysis: O(1)as it is just initialisation

        Space complexity:
            Input space analysis: O(1)as it is just initialisation
            Aux space analysis: O(1)as it is just initialisation
        """
        self.u = u
        self.v = v
        self.flow = 0
        self.capacity = capacity
        self.lower_bound = lower_bound
        self.balance = 0
        self.forward_residue_flow = 0
        self.backward_residue_flow = 0
        self.forward = None
        self.backward = None
        self.company_balance = 0
    

def allocate(preferences, officers_per_org, min_shifts, max_shifts):
    """ 
        Allocate security officer to the company that needed security officer
        and allocate to their shift and days accordingly

        precondition: a valid list of list preferences that indicates officer preferences
                      a valid list of list officers per org that indicates the amount of security needed by the company
                      an integer indicating minimum shift of each officer
                      an integer indicating maximum shift of each officer

        postcondition: return a list that allocate officer to the company if is feasible
                       return None if is not feasible

        Input:
            a list of list preferences that indicates officer preferences
            a list of list officers per org that indicates the amount of security needed by the company
            an integer indicating minimum shift of each officer
            an integer indicating maximum shift of each officer

        Return:
            None

        Approach description:
        A flow network with demands and lower bound is created, since every officer has minimum shift and maximum shift, so i 
        have a vertex at the start of the flow network which is min max node. Min max node will provide a lower bound and capacity of 
        min_shift and max_shift. So this min max node will connect to every single officer providing the lowerbound and capacity.
        and each officer will have an edge to every day with a lower bound of 0 and capacity of 1. So each officer can only work for one day
        After that, each day have 3 shifts. So each day will have 3 edges to 3 shifts. The lower bound will be 0 and the capacity will be the 
        number of officers that signed up for the shift. After that each shift will have an edge to each company with lower bound of company 
        needed for that particular shift and capacity of company needed for that particular shift. So for it to be possible every shift 
        needs to have a lower bound of the capacity to make sure the company get the amount of security they requested. After that
        the company node will have positive demand of total security that requested by all company. Min max node have negative demand with
        the total security that that requested by all company. After assigning all that, i resolved the lower bound of every edge using -incoming
        + outgoing method. After that i remove all the demands by connecting source to all the negative demand with the balance of the absolute value of negative demand
        positive demand will be connecting to the sink where the balance is the positive demand. Then i run ford fulkerson to ensure that the assignment is feasible
        If is feasible ill return an allocate list of security officer to the company. If is not feasible ill return None. 

        Time complexity:
            Best case analysis: O(M + N) where M is the number of company and N is the number of officer
                                inserting shift needed by company is O(M) where M is the number of company 
                                as i need to loop through every company to insert
                                connecting officer with days is O(30N) where N is the number of officer
                                as i need to loop through every officer and create an edge with all the days
                                so is O(30N) = O(N). Connecting day with shift is O(90) as looping through
                                30 days and each day connect with 3 shift which is 90 in total so is O(1) 
                                since it is constant. Next is connecting shift with the company which is 
                                O(M*30*3) = O(M) where M is the number of company as connecting 90 shifts 
                                with M company so it is O(M). Next connecting min max node to the officer
                                is O(N) where N is the number of officer as i need to loop through N and 
                                create N edges between min max node and officer. next is connecting source and sink
                                is O(N + M) where N is the number of Officer and M is the number of company as i looping through all the vertex and
                                the rest are constant so N and M is the only variable that dominates the complexity
                                Running ford fulkerson is O(N + M) where N is the number of officer in the graph
                                and M is the number of company. As going through all the path is same as going through 
                                all the nodes so is O(N+M). If the flow network is not feasible it returns none so the 
                                overall complexity is O(M + N) where M is the number of company and N is the number
                                of officer
                                

            Worst case analysis: O(M * N * N) where M is the number of company and N is the number of officer.
                                 the process same as best case but instead of return None. It returns a list where
                                 the security will be allocate to the company and shift accordingly. The complexity is 
                                 O(M * N * N) due to we need to traverse all path that can lead to that all the company
                                 to assign to the security officer in the list. So is O(M * N * N)

        Space complexity:
            Input space analysis: O(N + M) where N is the length of preferences and M is the length of officer_per_orgs
            Aux space analysis: connecting officer with day is O(30N) where N is the number of officer, we need to connect
                                officer with days so we need 30N edges which is O(N). O(1) for connecting day with shift as it is
                                constant. Connecting shift with company is O(M*30*3) where M is the number of company.
                                as we are creating M*30*3 of edges while connecting them. Connecting min max to officer is 
                                O(N) where N is the number of officers as we are creating N number of edges.
                                Connecting source and sink is O(N+M) where N is the number of officers and M is the number of company
                                as we are creating edges to all the nodes. So overall best case aux space complexity is O(N + M). If 
                                it returns none
            
            worst case aux space: O(M*N*N) where M is the number of company and N is the number of officer. The process is same as
                                  best case but if it is feasible it will loop through all the path to the company and append each 
                                  shift to the list so it will be O(M*N*N) in worst case
    """
    day = 30
    flow_reduce = Graph(preferences, officers_per_org, min_shifts, max_shifts, day)
    source = flow_reduce.vertices[0]
    sink = flow_reduce.vertices[len(flow_reduce.vertices) - 1]
    
    #inserting the shift needed by company
    #time:O(M) where M is the number of company
    #space: O(1) in place
    for i in range(len(officers_per_org)):
        flow_reduce.vertices[flow_reduce.shift_count + 1 + i].id = "Company" + str(i)
        flow_reduce.vertices[flow_reduce.shift_count + 1 + i].shift_needed = officers_per_org[i]

    #time: O(30N) where N is the number of officer
    #space: O(30N) where N is the number of officer 
    connect_officer_with_day(flow_reduce, preferences)
    

     #time O(30*3) = O(90) = O(1) as days and shift are constant
    #space O(30*3) = O(90) = O(1) as creating edge from day to shift and it is constant
    connect_day_with_shift(flow_reduce, preferences)

    #time O(M*30*3) = O(M) where M is the number of company 
    #space = O(M*30*3) = O(M) creating an edge and append to list for each shift to company so is O(M) 
    #where M is the number of company
    connect_shift_with_company(flow_reduce)

     #Time O(N) where N is the number of officer in the graph
    #space: O(N) where N is the number of officer, creating edge to all the officer
    connect_min_max_to_officer(flow_reduce, min_shifts, max_shifts)

    #time O(N + M) where N is number of officer and M is company number
    #space O(N + M)
    connect_source_sink(flow_reduce, source, sink)
    
    
    #Time: O(N + M) where N is the number of officer in the graph and M is the number of company in the graph
    #Space: O(N+M) if all the node need to have an edge with the source and sink
    feasible, flow = ford_fulkerson(flow_reduce)

    if(feasible):
        allocation = []
        #allocation[securityO1][company01][day][shift1] = 1
        day_list = []
        officer_list = []
        company = []
        shift_list = [0,0,0] 
        counter = 0  #counter of officer
        #Time O(M*N*N) where M is the number of company, N is the number of officer
        #space O(M*N*N) where M is the number of company and N is the number of officer
        #as going through all the edges and the vertex
        for i in range(flow_reduce.min_max_node + 1, flow_reduce.officer_count + 1):
            for officer in flow_reduce.vertices[i].edges:
                if officer.flow == 1:
                    for day in officer.v.edges:
                        if day.flow >= 1:
                            day.flow -= 1
                            
                            for shift in range(len(day.v.edges)):
                                if day.v.edges[shift].flow >= 1 and preferences[counter][shift]:
                                    day.v.edges[shift].flow -= 1
                                    shift_list[shift] = 1
                                    
                            day_list.append(shift_list)
                            shift_list = [0,0,0] #reset the shift list
                        company.append(day_list)
                else:
                    day_list.append(shift_list)
                officer_list.append(day_list)
                day_list = []
            allocation.append(officer_list)
            officer_list = []
            counter += 1

        

        return allocation

    return None

def connect_officer_with_day(flow_reduce, preferences):
    """ 
        connect officer with day 

        precondition: a valid graph with vertex and connected edges

        postcondition: graph where officer has connected edge to the day 

        Input:
            flow_reduce: a graph that connect with vertex and edges
            
        Return:
            None

        Time complexity:
            Best case analysis: time: O(30N) where N is the number of officer
                                connecting each officer to each day which is 30 a constant
                                so is O(N)

            Worst case analysis: same as best case

        Space complexity:
            Input space analysis: O(1) a graph object 
            Aux space analysis: space O(30N) where N is the number of officer 
                                creating 30 edges for every N so is O(N)


    """
    #connect officer with days
     #time: O(30N) where N is the number of officer
     #space: O(30N) where N is the number of officer 
    for i in range (len(preferences)):
        officer_index = flow_reduce.min_max_node + 1 + i
        flow_reduce.vertices[officer_index].id = "Officer" + str(i)
        #loop through the number of days which is O(30) = O(1) time
        #space: O(1) in place
        for day in range(flow_reduce.officer_count + 1, flow_reduce.days_count + 1):
            #looping through all the days and create edge from officer to day with lowerbound of 0 and 1 cap
            #edge balance is cap - lower bound
            #assign a forward and backward edge for residual graph later
            edge = Edge(flow_reduce.vertices[officer_index], flow_reduce.vertices[day], 0, 1)
            flow_reduce.vertices[officer_index].edges.append(edge)
            edge.balance = edge.capacity - edge.lower_bound
            edge.backward = Edge(flow_reduce.vertices[day], flow_reduce.vertices[officer_index])
            edge.forward = Edge(flow_reduce.vertices[officer_index], flow_reduce.vertices[day])
            edge.forward.forward_residue_flow = edge.balance

def connect_day_with_shift(flow_reduce, preferences):
    """ 
        connect day with shift each day can have 3 shift

        precondition: a valid graph with vertex and connected edges

        postcondition: graph where day has connected shift to the company

        Input:
            flow_reduce: a graph that connect with vertex and edges
            
        Return:
            None

        Time complexity:
            Best case analysis: time O(30*3 + N) = O(90N) = O(N) where N is the number of officer as days and shift are constant so
                                connect edges between them are constant

            Worst case analysis: same as best case

        Space complexity:
            Input space analysis: O(1) a graph object 
            Aux space analysis: space O(30*3) = O(90) = O(1) as creating edge from day to shift and it is constant


    """
    shift_preference = [0,0,0]
    #count preference
    #getting count for each shift preference
    #time O(N*3) = O(N) where N is the number of officer as each preference has 3 shift so is N*3
    #space O(1)
    for i in range(len(preferences)):
        #time: O(3) = O(1)
        for j in range(len(preferences[i])):
            shift_preference[j] += preferences[i][j] 

    #connect days with shift
     #getting starting index and ending index for shift as each day has its own 3 shift
    start_shift_index = flow_reduce.days_count + 1
    ending_shift_index = start_shift_index + 3 #exclusive
    #time O(30*3) = O(90) = O(1) as days and shift are constant
    #space O(30*3) = O(90) = O(1) as creating edge from day to shift and it is constant
    for day_index in range(flow_reduce.officer_count + 1, flow_reduce.days_count + 1):
        counter = 0
        #time O(3) = O(1) as shift is constant
        for j in range(start_shift_index, ending_shift_index):
            #creating an edge between day and shift 
            #balance =  capacity - lower bound and have a forward for residual graph later
            #using counter as a shift preference counter it will reset to 0 every time it enters this loop
            #forward has the same flow with balance
            edge = Edge(flow_reduce.vertices[day_index], flow_reduce.vertices[j], 0, shift_preference[counter])
            flow_reduce.vertices[day_index].edges.append(edge)
            edge.balance = edge.capacity - edge.lower_bound
            edge.forward = Edge(flow_reduce.vertices[day_index], flow_reduce.vertices[j])
            edge.backward = Edge(flow_reduce.vertices[j], flow_reduce.vertices[day_index])
            edge.forward.forward_residue_flow = edge.balance
            counter += 1
        #getting the next shift set for the next day
        start_shift_index = ending_shift_index
        ending_shift_index += 3

def connect_shift_with_company(flow_reduce):
    """ 
        connect shift to each company

        precondition: a valid graph with vertex and connected edges

        postcondition: shift is connected to each company with approperiate lowerbound and capacity

        Input:
            flow_reduce: a graph that connect with vertex and edges
            
        Return:
            None

        Time complexity:
            Best case analysis: time O(M*30*3) = O(M) where M is the number of company 
                                as creating an edge from shift to each company and shift is constant which is 90
                                so is O(M) 

            Worst case analysis: same as best case

        Space complexity:
            Input space analysis: O(1) a graph object 
            Aux space analysis: space:= O(M*30*3) = O(M) creating an edge and append to list for each shift to company so is O(M) space 
                                where M is the number of company

    """
    #connect shift with company
    #time connect shift with company each company will be connected with all shift
    #time O(M*30*3) = O(M) where M is the number of company 
    #space = O(M*30*3) = O(M) creating an edge and append to list for each shift to company so is O(M) 
    #where M is the number of company
    for i in range(flow_reduce.shift_count + 1, flow_reduce.companies_count + 1):
        company = flow_reduce.vertices[i]
        shift_count = 0
        #O(30 * 3) = O(1) as days and shift are constant so is O(1)
        for j in range(flow_reduce.days_count + 1, flow_reduce.shift_count + 1):
            #creating each edge from shift to company 
            #use counter to find the corresponding needed
            #assign demand to min max node which is the second node in the list
            #create forward and backward for residual later
            shift = flow_reduce.vertices[j]
            shift_count = shift_count % 3
            edge = Edge(shift, company, company.shift_needed[shift_count], company.shift_needed[shift_count])
            edge.company_balance = company.shift_needed[shift_count]
            shift.edges.append(edge)
            shift.demand += company.shift_needed[shift_count]

            #assign demand to min max node
            flow_reduce.vertices[1].demand -= company.shift_needed[shift_count]

            edge.balance = edge.capacity - edge.lower_bound
            edge.forward = Edge(shift, company)
            edge.backward = Edge(company, shift)
            edge.forward.forward_residue_flow = edge.balance
            shift_count += 1

def connect_min_max_to_officer(flow_reduce, min_shifts, max_shifts):
    """ 
        connect min max vertex to officer to state the lower bound and capacity

        precondition: a valid graph with vertex and connected edges

        postcondition: the officer vertex are connected by min max node to indicate the lower bound flow and capacity

        Input:
            flow_reduce: a graph that connect with vertex and edges
            
        Return:
            None

        Time complexity:
            Best case analysis: O(N) where N is the number of officer in the graph as looping each and every
                                officer to create an edge from min max node to the officer of is O(N)

            Worst case analysis: same as best case

        Space complexity:
            Input space analysis: O(1) a graph object 
            Aux space analysis: space: O(N) where N is the number of officer, creating edge to all the officer so is O(N)

    """
    #connect min max node to every officer
    #Time O(N) where N is the number of officer in the graph
    #space: O(N) where N is the number of officer, creating edge to all the officer
    for i in range(flow_reduce.min_max_node + 1, flow_reduce.officer_count + 1):
        #connect min max node to each officer of the graph to provide the lower bound
        #create forward and backward for residual 
        min_max_node = flow_reduce.vertices[1]
        officer = flow_reduce.vertices[i]
        edge = Edge(min_max_node, officer, min_shifts, max_shifts)
        min_max_node.edges.append(edge)
        officer.demand -= min_shifts
        min_max_node.demand += min_shifts
        edge.balance = edge.capacity - edge.lower_bound
        edge.forward = Edge(min_max_node, officer)
        edge.backward = Edge(officer, min_max_node)
        edge.forward.forward_residue_flow = edge.balance

            
def connect_source_sink(flow_reduce, source, sink):
    """ 
        connect source and sink to negative and positive demand node

        precondition: a valid graph with vertex and connected edges, a vertex representing source
                      a vertex representing sink

        postcondition: the necessary vertex are connected to source and sink

        Input:
            flow_reduce: a graph that connect with vertex and edges
            source: a vertex representing source
            sink: a vertex representing sink
        Return:
            None

        Time complexity:
            Best case analysis: O(N + M) where N is the number of officer and M is the number of company
                                looping through all the vertex to identify which needs to connnect is O(N+M+(30*3)+30+3) = O(N + M)

            Worst case analysis: same as best case

        Space complexity:
            Input space analysis: O(1) a graph object and vertex
            Aux space analysis: O(N + M) as adding edges of all vertex to source and sink is O(N + M)
    """
    #connect source and sink
    #time: O(N+M+(30*3)+30+3) = O(N + M) where N is the number of officer and M is company number
    #Space: O(N+M) if all the node need to have an edge with the source and sink
    for i in range(len(flow_reduce.vertices) - 1):
        node = flow_reduce.vertices[i]
        #creating edge from source to node
        if(node.demand < 0):
            edge = Edge(source, node)
            source.edges.append(edge)
            edge.balance = abs(node.demand)
            edge.forward_residue_flow = edge.balance
            edge.forward = Edge(source, node)
            edge.backward = Edge(node, source)
            edge.forward.forward_residue_flow = edge.balance
            node.demand = 0
        
        #creating edge from node to sink
        #have forward balance and everything
        elif (node.demand > 0):
            edge = Edge(node, sink)
            edge.balance = abs(node.demand)
            edge.forward_residue_flow = edge.balance
            node.edges.append(edge)
            edge.forward = Edge(node, sink)
            edge.backward = Edge(sink, node)
            edge.forward.forward_residue_flow = edge.balance
            node.demand = 0

def ford_fulkerson(graph):
    """ 
        obtain the flow and check the graph is it feasible

        precondition: a valid graph that has edge and vertex connected with each other
        postcondition: able to determine it is feasible or no and return flow

        Input:
            graph: a valid graph that has edge and vertex connected with each other
        Return:
            boolean representing feasible and flow of it

        Time complexity:
            Best case analysis: O(V + E) where V is the number of vertex in the graph and E is the total edges in the graph
                                loop through all the path will be O(V + E) but we need to loop 90 times as the maximum it can be
                                is the number of shift that connected to sink. The shift is constant which is 90. So it is O(90(V+E)).
                                traversing from node to node in the path maximum we can traverse is V to get the min flow of the path where V is the number of 
                                vertex in the graph. So now is O(90(V + E + V)), now we need to update the residue flow according to the min flow
                                So we need to traverse maximum V so now is O(90(V + E + V + V)). Next we need to check is the graph feasible, checking 
                                as it loop through all the vertex that source connected and all the vertex that connected to the sink. So maximum is V 
                                where V is the number of vertex. So now the complexity is O(90(V + E + V) + V) = O(V + E) after removing all the constant

            Worst case analysis: Same as best case

        Space complexity:
            Input space analysis: O(1) a graph object
            Aux space analysis: O(V) as in bfs we need O(V) to store the vertex
    """
    flow = 0
    #while there is a path
    path_exist = True
    min_flow = 0
    
    #maximum path is the number of the node that connected to sink
    #O(90(2V + E + V + V)) where N is the total vertex in graph and V is the number of vertex in graph
    while path_exist:
        flow += min_flow
        min_flow = 30

        #time O(2V + E)
        #Space O(V)
        path = bfs(graph)
        path_reference = path
        path_exist = path

        #traverse from sink to source 
        #time O(V) where V is the number of vertex in the graph
        while path is not None:
            if path is not None:
                min_flow = min(min_flow, path.forward.forward_residue_flow)
                path = path.u.parent
        
        #loop until no path
        #time O(V) where V is the number of vertex in the graph
        while path_reference is not None:
            if path_reference is not None:
                path_reference.forward.forward_residue_flow -= min_flow
                path_reference.forward_residue_flow -= min_flow
                path_reference.backward.backward_residue_flow += min_flow
                path_reference.backward_residue_flow += min_flow
                path_reference.flow += min_flow
                path_reference = path_reference.u.parent
    #time: O(V)
    #space: O(1)
    feasible = check_feasible(graph)
  
    return feasible, flow

def check_feasible(graph):
    """ 
        check it is feasible for the given graph

        precondition: a valid graph that has edge and vertex connected with each other
        postcondition: able to determine it is feasible or no

        Input:
            graph: a valid graph that has edge and vertex connected with each other
        Return:
            boolean representing feasible

        Time complexity:
            Best case analysis: O(V + 90) = O(V) where V is the number of vertex in the graph. as we need to loop through all the of vertex
                                we need to loop through all of the shift so is O(V + 90). Since 90 is constant so O(V)

            Worst case analysis: Same as best case

        Space complexity:
            Input space analysis: O(1) a graph object
            Aux space analysis: O(1) in place
    """
    feasible = True
    #flow through each edges in source
    #if is not equal then not feasible
    #O(V) where V is the total number of vertex
    for i in graph.vertices[0].edges:
        if i.balance != i.flow:
            feasible = False

    #flow through each edges of shift if is not max then not feasible
    #O(30 * 3) = O(1) it has 90 shift so is constant
    for i in range(graph.days_count + 1, graph.shift_count + 1):
        shift = graph.vertices[i]
        if (shift.edges[-1].balance != shift.edges[-1].flow):
            feasible = False

    #check min max node if there is edge to sink then check 
    if graph.vertices[1].edges[-1].v.id == str(len(graph.vertices) - 1):
        if (i.v.balance != i.v.flow):
               feasible = False
    
    return feasible


def bfs(graph):
    """ 
        breath first search that traverse all of the edges and vertex

        precondition: a valid graph that has edge and vertex connected with each other
        postcondition: able to return a path or None if no path

        Input:
            graph: a valid graph that has edge and vertex connected with each other
        Return:
            an edge between the node and the sink

        Time complexity:
            Best case analysis: O(N) where N is the number of vertex as traversing each of them to 
                                assign false to discovered and visited. When the source is connected to sink at the 
                                first edge. So when it encountered sink after source then it return so is O(N)

            Worst case analysis: O(2V + E) = O(V + E)where N is the number of vertex as traversing each of them to assign, is O(V)
                                 V is the number of vertex in the graph and E is the number of edge in the graph. 
                                 it is worst case when it traverse every single node before meeting the sink. So is O(N+V+E)

        Space complexity:
            Input space analysis: O(1) a graph object
            Aux space analysis: O(V) where V is the number of vertex in the graph for the queue to store the node, if a vertex connected to V number of node it
                                will store V number of node into the queue before it started pop
    """
    #assign a queue
    queue = deque()
    #setting all the graph vertices discovered and visited to false
    #time: O(V) where V is the number of vertex
    #space : O(1) in place
    for i in graph.vertices:
        i.discovered = False
        i.visited = False
        
    #queue append source first
    queue.append(graph.vertices[0])

    while len(queue) > 0:
        #get the node to traverse
        u = queue.popleft() 
        u.visited = True
        
        #traverse all the edges of the node
        #time: O(E) where E is the number of edge of the node
        #space :O(1)
        for edge in u.edges:
            #if the forward is 0 then don go to the edge
            if edge.forward.forward_residue_flow > 0:
                v = edge.v
                if v.discovered == False:
                    queue.append(v)
                    v.discovered = True
                    v.parent = edge
                    #if traverse until a sink then return 
                    if(v.id == str(len(graph.vertices) - 1)):
                        return v.parent
        
    return None

