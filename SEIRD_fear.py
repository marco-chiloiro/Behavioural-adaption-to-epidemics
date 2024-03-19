import numpy as np   
import math
import matplotlib
import matplotlib.pyplot as plt
import networkx as nx
import random
import copy
import scipy
from tqdm import tqdm
import pandas as pd
import seaborn as sns

#initialization of the infection in a graph
def initialize_SEIRD(G, I0,visualize=True):
    """
    Initializes the SIR model on a given graph.

    Parameters:
    G (networkx.Graph): The graph on which the SIR model will be initialized.
    I0 (int): The number of initially infected nodes.

    Returns:
    G (networkx.Graph): The graph with the SIR model initialized.
    S (list): The list of susceptible nodes.
    E (list): The list of exposed nodes.
    I (list): The list of infected nodes.
    R (list): The list of recovered nodes.
    D (list): The list of dead nodes.
    """
    # state 0 = susceptible
    # state 1 = infected
    
    #set attribute state to all nodes and set it to 0 (susceptible)
    nx.set_node_attributes(G, 0, 'state')

    #set the state of I0 randomly chosen nodes to 1 (infected)
    for i in range(I0):
        #extract a random node (not really random because of the node degree)
        node = random.choice(list(nx.nodes(G)))

        #if it is already infected extract another one
        while G.nodes[node]['state'] == 1:
            node = random.choice(list(nx.nodes(G)))

        #set the state to infected
        G.nodes[node]['state'] = 1

    N = G.number_of_nodes() 

    S = [N-I0]
    E = [0] 
    I = [I0]
    R = [0]
    D = [0]

    #initialize next state
    nx.set_node_attributes(G, 0, 'new_state')

    if visualize:
        #visualize the initial conditions
        pos = nx.spring_layout(G)
        nx.draw(G, pos, node_size=10, node_color='blue')
        nx.draw_networkx_nodes(G, pos, node_size=10, nodelist=[i for i in nx.nodes(G) if G.nodes[i]['state']==1], node_color='red')
        plt.title('Initial conditions with infected nodes in red', fontsize=8)

        
    return G, S, E, I, R, D



def SEIRD(G, S, E, I, R, D, sigma ,beta,mu, theta, simulation_time=1, info=True):

    # 0 = susceptible
    # 1 = exposed
    # 2 = infected
    # 3 = recovered
    # 4 = dead

    #beta is the infection rate
    #mu is the recovery rate
    #sigma is the incubation period E->I
    #theta is the probability of death
    

    #set the simulation

    t=0
    
    #perform the simulation
    while True:

        #loop over all nodes
        for i in nx.nodes(G):
            #if the node is susceptible
            if G.nodes[i]['state'] == 0:
                #loop over all neighbors
                for j in nx.all_neighbors(G, i):
                    #if the neighbor is infected
                    if G.nodes[j]['state'] == 2:
                        #infect with probability beta (S->E)
                        if np.random.random() < beta:
                            G.nodes[i]['new_state'] = 1
                            break #break the loop over neighbors when a neighbor is infected
                            
                        else:
                            G.nodes[i]['new_state'] = 0
                
            #if the node is exposed
            elif G.nodes[i]['state'] == 1:
                #become ill with probability sigma (E->I)
                if np.random.random() < sigma:
                    G.nodes[i]['new_state'] = 2
                else:
                    G.nodes[i]['new_state'] = 1
                    
            #if the node is infected recover with probability mu (I->R) or die with probability theta (I->D)
            elif G.nodes[i]['state'] == 2:
                if np.random.random() < mu:
                    G.nodes[i]['new_state'] = 3
                elif np.random.random() < theta:
                    G.nodes[i]['new_state'] = 4
                else:
                    G.nodes[i]['new_state'] = 2
       
                
        #update the state of all nodes
        for i in nx.nodes(G):
            G.nodes[i]['state'] = G.nodes[i]['new_state']

        #update the time
        t += 1

        s,e,infected,r,d = 0,0,0,0,0
        #count the number of nodes in each state
        for i in nx.nodes(G):
            
            if G.nodes[i]['state'] == 0:
                s += 1
            elif G.nodes[i]['state'] == 1:
                e += 1
            elif G.nodes[i]['state'] == 2:
                infected += 1
            elif G.nodes[i]['state'] == 3:
                r += 1
            elif G.nodes[i]['state'] == 4:
                d += 1

        #add the new state to the list
        S.append(s)
        E.append(e)
        I.append(infected)
        R.append(r)
        D.append(d)
        
        #print the number of infections every 100 time steps (gives an overview of the simulation)
        if t % 100 == 0:
            print(t, infected)

        #if there are no more infected nodes or reach the simulation time stop the simulation
        if infected == 0 or t == simulation_time:
            break

    if info:
        print("converged in: {} steps".format(t))
        print("maximal number of infected nodes: {}".format(max(I)), "at time step: {}".format(I.index(max(I))))

    return G,t,S, E, I, R, D #,n_rewired

def plot_seird(G_final, S, E, I, R, D, graph=True):
    

    

    if graph:

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 3))
        # Plot 2: Network graph
        pos = nx.spring_layout(G_final)
        nx.draw(G_final, pos, node_size=10, node_color='blue', ax=ax1)
        nx.draw_networkx_nodes(G_final, pos, node_size=10, nodelist=[i for i in nx.nodes(G_final) if G_final.nodes[i]['state']==1], node_color='red', ax=ax1)
        ax1.set_title('final graph', fontsize=8)


        # Plot 1: SIR curves
        ax2.plot(range(len(S)), S, label='susceptible', color='blue')
        ax2.plot(range(len(E)), E, label='exposed', color='orange')
        ax2.plot(range(len(I)), I, label='infected', color='red')
        ax2.plot(range(len(R)), R, label='recovered', color='green')
        ax2.plot(range(len(D)), D, label='dead', color='black')
        ax2.legend()
        ax2.set_xlabel('time')

    else:

        # Plot 1: SIR curves
        plt.plot(range(len(S)), S, label='susceptible', color='blue')
        plt.plot(range(len(E)), E, label='exposed', color='orange')
        plt.plot(range(len(I)), I, label='infected', color='red')
        plt.plot(range(len(R)), R, label='recovered', color='green')
        plt.plot(range(len(D)), D, label='dead', color='black')
        plt.legend()
        plt.xlabel('time')
    
    plt.show()

    #idea is to link the fear level to the probability of infection and the probability of death
#to start with we can link the fear level to the probability of death
#the fear level is a normal distribution with mean and standard deviation as parameters


def starting_fear(F, mean,std,visualize=True):

    #idea is that each node has a fear value that is sampled from a normal distribution. 

    
    #set attribute state to all nodes and set it to 0 (no fear)
    nx.set_node_attributes(F, 0, 'state')

    for i in nx.nodes(F):
        #sample the fear value for each node from a normal distribution
        F.nodes[i]['state'] = np.random.normal(mean, std)


    if visualize:
        #visualize the initial conditions
        pos = nx.spring_layout(F)
        nx.draw(F, pos, node_size=10, node_color='blue')
        nodelist = []
        
        for i in nx.nodes(F):
            if F.nodes[i]['state'] > 0.5:
                nodelist.append(i)

        nx.draw_networkx_nodes(F, pos, node_size=10, nodelist=nodelist, node_color='red')
        plt.title('initial fear level', fontsize=8)
        plt.show()
        

        
    return F



# fear based on the number of infected, recovered and dead nodes. each node has a fear value that is updated at each time step
# we link fear to the betas

#higher deaths, higher infections -> lower beta


def fear_level(neighbors, G ,a=1, b=1):
    """
    Calculate the fear value which is based on the level of infected or dead neighbors.

    Parameters:
    - neighbors: list
        The list of neighbors of a node.
    - a: float
        The parameter to weight the infected.
    - b: float
        The parameter to weight the deaths.

    Returns:
    - fear: float
        1 is no fear, closer to zero is higher fear.
    """
    #calculate the fear level based on the number of infected and dead neighbors in the epidemic layer G
    #ik are the fraction of infected neighbors

    ik= sum([1 for n in neighbors if G.nodes[n]['state'] == 2])/len(neighbors)


    #dk are the fraction of dead neighbors
    dk= sum([1 for n in neighbors if G.nodes[n]['state'] == 4])/len(neighbors)
    
    #calculate the fear level as a exponential of the fraction of infected and dead neighbors
    fear = np.exp(-(a*ik + b*dk))

    
    return fear


#select the node with the highest degree in a network
def highest_degree_node(G):
    """
    Select the node with the highest degree in the network.

    Parameters:
    - G: NetworkX graph object
        The graph representing the population and connections between individuals.

    Returns:
    - node: int
        The node with the highest degree in the network.
    """
    #get the node with the highest degree
    node = max(G.degree(), key=lambda x: x[1])[0]
    return node


#link seird to fear level

def SEIRD_fear(G,F,S,E,I,R,D,sigma ,beta,mu, theta, simulation_time=1, info=True,a=1,b=1):
    
        # 0 = susceptible
        # 1 = exposed
        # 2 = infected
        # 3 = recovered
        # 4 = dead
    
        #beta is the infection rate
        #mu is the recovery rate
        #sigma is the incubation period E->I
        #theta is the probability of death
        
        #F is the fear net
        #G is the epidemic net
    
        #set the simulation
    
        t=0

        # we want to follow the fear level of one node for each degree value

        degree_sequence = [d for n, d in F.degree()]

        binwidth=1
        histogram=np.histogram(degree_sequence, bins=np.arange(min(degree_sequence), max(degree_sequence) + binwidth, binwidth))

        list_degree = histogram[1] #bins
        weights = histogram[0] #frequency
        weights = weights[weights != 0] #remove the zeros

        #shuffle the list of nodes
        nodes = list(F.nodes())
        random.shuffle(nodes)

        #loop over the nodes and save the node with the degree in the list
        node_to_follow = []
        for i in list_degree:
            for j in nodes:
                if F.degree(j) == i:
                    node_to_follow.append(j)
                    break

        #create an array of size t x size of the list of nodes to follow the fear level of each node
        fear_t = np.zeros((simulation_time,len(node_to_follow)))


        #perform the simulation
        while True:\
        
            #----------update the fear level--------------------------------------------------------------------------------
            
            #loop over the nodes in the social network and update the fear level of each node
            
            for i in nx.nodes(F):

                #calculate the fear level based on the number of infected and dead neighbors in the social net
                neighbors = list(nx.neighbors(F, i)) # take the social neibourghs of the node and calculate the fear level with their state in the epidemic net
                fear = fear_level(neighbors, G, a=a, b=b) 
                    
                #not considering the previous fear level 
                F.nodes[i]['state'] = fear    

                if i in node_to_follow:
                    fear_t[t][node_to_follow.index(i)] = fear

        
            #----------update the state of the epidemic net-----------------------------------------------------------------  

            #loop over all nodes in the epidemic net
            for i in nx.nodes(G):
                #if the node is susceptible
                if G.nodes[i]['state'] == 0:
                    #loop over all neighbors
                    
                    for j in nx.all_neighbors(G, i):
                        #if the neighbor is infected
                        if G.nodes[j]['state'] == 2:
                            #infect with probability beta (S->E) and the probability of infection is linked to the fear level
                            #calculate the fear level based on the number of infected and dead neighbors in the social net
                            #neighbors=list(nx.all_neighbors(F, i))
                            fear = F.nodes[i]['state']
                            #infect with probability beta*fear
                            if np.random.random() < beta*fear: 

                                G.nodes[i]['new_state'] = 1
                                
                                break #break the loop over neighbors when a neighbor is infected
                                
                            else:
                                G.nodes[i]['new_state'] = 0
                    
                #if the node is exposed
                elif G.nodes[i]['state'] == 1:
                    #become ill with probability sigma (E->I)
                    if np.random.random() < sigma:
                        G.nodes[i]['new_state'] = 2
                    else:
                        G.nodes[i]['new_state'] = 1
                        
                #if the node is infected recover with probability mu (I->R) or die with probability theta (I->D)
                elif G.nodes[i]['state'] == 2:

                    r=np.random.random()
                    r2=np.random.random()+1
                    
                    while r < mu and r2 < theta+1:
                        r=np.random.random()
                        r2=np.random.random()+1
                    
                    if r < mu:
                        G.nodes[i]['new_state'] = 3
                    elif r2 < theta+1:
                        G.nodes[i]['new_state'] = 4
                    else:
                        G.nodes[i]['new_state'] = 2
        
                    
            #update the state of all nodes
            for i in nx.nodes(G):
                G.nodes[i]['state'] = G.nodes[i]['new_state']


            #update the time
            t += 1
    
            s,e,infected,r,d = 0,0,0,0,0
            #count the number of nodes in each state
            for i in nx.nodes(G):
                
                if G.nodes[i]['state'] == 0:
                    s += 1
                elif G.nodes[i]['state'] == 1:
                    e += 1
                elif G.nodes[i]['state'] == 2:
                    infected += 1
                elif G.nodes[i]['state'] == 3:
                    r += 1
                elif G.nodes[i]['state'] == 4:
                    d += 1

            #add the new state to the list
            S.append(s)
            E.append(e)
            I.append(infected)
            R.append(r)
            D.append(d)
            
            #print the number of infections every 100 time steps (gives an overview of the simulation)
            if t % 1000 == 0:
                print(t, infected)

            #if there are no more infected nodes or reach the simulation time stop the simulation
            if infected == 0 or t == simulation_time:
                break

        if info:
            print("converged in: {} steps".format(t))
            print("maximal number of infected nodes: {}".format(max(I)), "at time step: {}".format(I.index(max(I))))

        fear_t = pd.DataFrame(fear_t, columns=node_to_follow)
        #cut when the epidemic is over
        fear_t = fear_t.iloc[:t]    
        weighted_fear = np.average(fear_t, axis=1, weights=weights)
        #append the weighted fear to the last column of the fear_t dataframe
        fear_t['weighted_fear'] = weighted_fear

        
        

        return G,t,S, E, I, R, D,fear_t #,n_rewired

def plot_fear(fear_t, F):
    
    #plot also the average
    plt.figure(figsize=(10,5))


    plt.plot(fear_t["weighted_fear"], label='weighted average fear level', color='blue', linestyle='dashed',zorder=10)
    #plot the fear level for every node
    plt.plot(fear_t, color='grey', alpha=0.3)

    #plot the fear level for the node with the highest degree
    #get the node with the highest degree
    node = highest_degree_node(F)
    #take the column of the fear level of the node with the highest degree
    fear_node = fear_t[node]
    plt.plot(fear_node, label='fear level of the node with the highest degree', color='red')


    plt.title('fear level of the nodes followed')
    plt.xlabel('time')
    plt.ylabel('fear level')
    plt.legend()
    plt.show()
    

def compare_fear(A,BA, I0, sigma,beta, mu,theta,sim_time, visualize=True,a=1,b=1):
    

    BAf=BA.copy()

    S,E,I,R,D=[],[],[],[],[]


    G_start,S_start,E_start,I_start,R_start,D_start=initialize_SEIRD(BA,I0=I0,visualize=False)


    S.append(S_start[0]) 
    I.append(I_start[0])
    R.append(R_start[0])
    E.append(E_start[0])
    D.append(D_start[0])

                                                    #E->I        S->E    I->R      I->D
    G_start,t,S,E,I,R,D = SEIRD(G_start,S,E,I,R,D,sigma=sigma,beta=beta,mu=mu, theta=theta,simulation_time=sim_time,info=False)

    if visualize:

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 3))

        # Plot 1: SEIRD curves
        ax1.plot(range(len(S)), S, label='susceptible', color='blue')
        ax1.plot(range(len(E)), E, label='exposed', color='orange')
        ax1.plot(range(len(I)), I, label='infected', color='red')
        ax1.plot(range(len(R)), R, label='recovered', color='green')
        ax1.plot(range(len(D)), D, label='dead', color='black')
        ax1.legend()
        ax1.set_xlabel('time')
        ax1.set_title('SEIRD curves')


    #initialize to 1 = no fear
    F = starting_fear(A,mean=1,std=0,visualize=False)

    Sf,Ef,If,Rf,Df=[],[],[],[],[]


    Gf_start,Sf_start,Ef_start,If_start,Rf_start,Df_start=initialize_SEIRD(BAf,I0=I0,visualize=False)



    Sf.append(Sf_start[0])
    If.append(If_start[0])
    Rf.append(Rf_start[0])
    Ef.append(Ef_start[0])
    Df.append(Df_start[0])
                                                                    #E->I   S->E    I->R   I->D
    Gf_start,tf,Sf,Ef,If,Rf,Df,fear_t = SEIRD_fear(Gf_start,F,Sf,Ef,If,Rf,Df,sigma=sigma,beta=beta,mu=mu, theta=theta,simulation_time=sim_time,info=False,a=a,b=b)

    if visualize:
        # Plot 2: SEIRD curves
        ax2.plot(range(len(Sf)), Sf, label='susceptible', color='blue')
        ax2.plot(range(len(Ef)), Ef, label='exposed', color='orange')
        ax2.plot(range(len(If)), If, label='infected', color='red')
        ax2.plot(range(len(Rf)), Rf, label='recovered', color='green')
        ax2.plot(range(len(Df)), Df, label='dead', color='black')
        ax2.legend()
        ax2.set_xlabel('time')
        ax2.set_title('SEIRD curves with fear')


        plt.show()
        plot_fear(fear_t, F)
    return G_start,Gf_start,I,If,t,tf




