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
def initialize_SIR(G, I0,visualize=True):
    """
    Initializes the SIR model on a given graph.

    Parameters:
    G (networkx.Graph): The graph on which the SIR model will be initialized.
    I0 (int): The number of initially infected nodes.

    Returns:
    G (networkx.Graph): The graph with the SIR model initialized.
    S (list): The list of susceptible nodes.
    I (list): The list of infected nodes.
    R (list): The list of recovered nodes.
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
    I = [I0]
    R = [0]

    #initialize next state
    nx.set_node_attributes(G, 0, 'new_state')

    if visualize:
        #visualize the initial conditions
        pos = nx.spring_layout(G)
        nx.draw(G, pos, node_size=10, node_color='blue')
        nx.draw_networkx_nodes(G, pos, node_size=10, nodelist=[i for i in nx.nodes(G) if G.nodes[i]['state']==1], node_color='red')
        plt.title('Initial conditions with infected nodes in red', fontsize=8)

        
    return G, S, I, R


def SIR_adaptive(G, S, I, R, beta, mu, adaptive=False, theta=0.5, simulation_time=1, info=True):
    """
    Simulates the SIR (Susceptible-Infected-Recovered) model on a given graph.

    Parameters:
    - G (networkx.Graph): The initial graph representing the network.
    - S (list): List to store the number of susceptible nodes at each time step.
    - I (list): List to store the number of infected nodes at each time step.
    - R (list): List to store the number of recovered nodes at each time step.
    - beta (float): Infection probability.
    - mu (float): Recovery probability.
    - adaptive (bool): Flag indicating whether to use adaptive rewiring.
    - theta (float): Probability of rewiring an edge during adaptive rewiring.
    - simulation_time (int): Maximum number of time steps to simulate.
    - info (bool): Flag indicating whether to print simulation information.

    Returns:
    - G_start (networkx.Graph): The initial graph before simulation.
    - G (networkx.Graph): The final graph after simulation.
    - t (int): Number of time steps taken to converge.
    - S (list): List of the number of susceptible nodes at each time step.
    - I (list): List of the number of infected nodes at each time step.
    - R (list): List of the number of recovered nodes at each time step.
    """
    

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
                    if G.nodes[j]['state'] == 1:
                        #infect with probability beta
                        if np.random.random() < beta:
                            G.nodes[i]['new_state'] = 1
                            break #break the loop over neighbors when a neighbor is infected
                            
                        else:
                            G.nodes[i]['new_state'] = 0
                
            #if the node is infected
            elif G.nodes[i]['state'] == 1:
                #recover with probability mu
                if np.random.random() < mu:
                    G.nodes[i]['new_state'] = 2
                else:
                    G.nodes[i]['new_state'] = 1
            #if the node is recovered
            elif G.nodes[i]['state'] == 2:
                G.nodes[i]['new_state'] = 2

        if adaptive:

            rewire= [] # link to rewire

            #loop over edges of the graph exluding self-loops 
            for i,j in nx.edges(G):

                #avoid self loop
                if i == j:
                     #print("self loop")
                     #skip the self loop
                    continue
                     
                #avoid double counting
                if i > j:
                    #print("double count")
                    continue
                
                #if the edge is between a susceptible and an infected node rewire the link to a 
                #susceptible node that is the neighbor of the infected node with probability theta

                if G.nodes[i]['state'] == 0 and G.nodes[j]['state'] == 1:

                    #TODO: probability of rewiring dependent on the degree of the infected node and the number of infected neighbors, how to weight them?
                    
                    #rewire with probability theta
                    if np.random.random() < theta:

                        #extract a random not infected node
                        node = random.choice(list(nx.nodes(G)))

                        #check that the node is not a neighor of the infected node and not the infected node itself
                        while G.nodes[node]['state'] == 1 or node in nx.neighbors(G,i) or node==i or node==j: # or node in nx.neighbors(G,j)?

                            node = random.choice(list(nx.nodes(G)))


                        #print("rewire {}-{} to {}-{}".format(i,j,i,node))
                        #save the link and rewire it after the loop
                        rewire.append((i,j,node))
                        #for debugging
                        #n_rewired.append((i,j,node))


            #rewire the links without self loops
            for i,j,node in rewire:
                if i!=j!=node:
                    G.remove_edge(i,j)
                    G.add_edge(i,node)
                    
            #remove self loops if any left
            #G.remove_edges_from(nx.selfloop_edges(G))

                
        #update the state of all nodes
        for i in nx.nodes(G):
            G.nodes[i]['state'] = G.nodes[i]['new_state']

        #update the time
        t += 1

        s,infected,r = 0,0,0
        #count the number of nodes in each state
        for i in nx.nodes(G):
            
            if G.nodes[i]['state'] == 0:
                s += 1
            elif G.nodes[i]['state'] == 1:
                infected += 1
            elif G.nodes[i]['state'] == 2:
                r += 1

        #add the new state to the list
        S.append(s)
        I.append(infected)
        R.append(r)
        
        #print the number of infections every 100 time steps (gives an overview of the simulation)
        if t % 1000 == 0:
            print(t, infected)

        #if there are no more infected nodes or reach the simulation time stop the simulation
        if infected == 0 or t == simulation_time:
            break

    if info:
        print("converged in: {} steps".format(t))
        print("maximal number of infected nodes: {}".format(max(I)), "at time step: {}".format(I.index(max(I))))

    return G,t,S,I,R #,n_rewired

def plot_sim(G_final, S, I, R):
    

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 3))

    # Plot 2: Network graph
    pos = nx.spring_layout(G_final)
    nx.draw(G_final, pos, node_size=10, node_color='blue', ax=ax1)
    nx.draw_networkx_nodes(G_final, pos, node_size=10, nodelist=[i for i in nx.nodes(G_final) if G_final.nodes[i]['state']==1], node_color='red', ax=ax1)
    ax1.set_title('final graph', fontsize=8)


    # Plot 1: SIR curves
    ax2.plot(range(len(I)), I, label='infected', color='red')
    ax2.plot(range(len(S)), S, label='susceptible', color='blue')
    ax2.plot(range(len(R)), R, label='recovered', color='green')
    ax2.legend()
    ax2.set_xlabel('time')
    
    plt.show()


def plot_net(G_final):


    #fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 3))

    # Plot 2: Network graph
    pos = nx.spring_layout(G_final)
    nx.draw(G_final, pos, node_size=10, node_color='blue')
    nx.draw_networkx_nodes(G_final, pos, node_size=10, nodelist=[i for i in nx.nodes(G_final) if G_final.nodes[i]['state']==1], node_color='red')

    nx.draw_networkx_nodes(G_final, pos, node_size=10, nodelist=[i for i in nx.nodes(G_final) if G_final.nodes[i]['state']==1], node_color='red')
    plt.title('final graph', fontsize=8)
    #plt.savefig('gif/final_graph_{}.png')
    plt.show()

def fitR0( y, start = 0, mu = 1/7, log = True, info=True):

    # x are timesteps
    # y is the number of infected nodes

    x =range(len(y))

    # fit only the data before the maximum number of infected nodes, i.e. the exponential growth phase
    n_points = int(y.index(max(y))-0.05*len(y)/2) #number of points to fit. if it does not work change 0.05*len(y)/2 (usually when the epidemic is slow to grow)

    xfit = x[start:start+n_points]
    yfit = y[start:start+n_points]
    
    def exp(X, I0, G):
        return I0*np.exp(G*X)
    
    # fit curve
    try:
        popt, _ = scipy.optimize.curve_fit(exp, xfit, yfit)
    except:
        #print("fit failed")
        return 0
        
    I0, G = popt
    
    if info:
        x_line = np.linspace(xfit[0], xfit[-1], 20)
        y_line = exp(x_line, I0, G)
        plt.scatter(x, y, c='navy', label='data')
        plt.plot(x_line, y_line, '--', color='red', label='exponential fit')
        if log:
            plt.yscale("log")
        plt.xlabel('time')
        plt.ylabel('$I$')
        plt.legend()
        plt.show()
    
        print("G =", G)
        print("\nR0 =", G/mu+1)

    return G/mu+1

## Static vs adaptive case on er and ba graphs

def compare(G, I0=10, beta=0.3, mu=0.1, theta=0.5, sim_time=1000, info=True):
    """
    Compare the SIR simulation results between a graph with adaptive behavior and a graph without adaptive behavior.

    Parameters:
    - G: NetworkX graph object
        The graph representing the population and connections between individuals.
    - I0: int, optional
        The initial number of infected individuals. Default is 10.
    - beta: float, optional
        The infection rate. Default is 0.3.
    - mu: float, optional
        The recovery rate. Default is 0.1.
    - theta: float, optional
        The threshold for adaptive behavior. Default is 0.5.
    - sim_time: int, optional
        The simulation time in days. Default is 1000.
    - info: bool, optional
        Whether to print additional information during the simulation. Default is True.

    Returns:
    - G: NetworkX graph object
        The graph representing the population and connections between individuals without adaptive behavior.
    - H: NetworkX graph object
        The graph representing the population and connections between individuals with adaptive behavior.
    """
    
    #make a copy of the graph (not sure if this is the right way to do it)
    H = G.__class__()
    H.add_nodes_from(G)
    H.add_edges_from(G.edges())

    S,I,R=[],[],[]
    
    G_start,S_start,I_start,R_start=initialize_SIR(G,I0=I0,visualize=info)

    S.append(S_start[0]) # N-I0
    I.append(I_start[0]) # I0
    R.append(R_start[0]) # 0

    

    G_start,t_stat,S,I,R = SIR_adaptive(G_start,S,I,R, beta=beta, mu=mu,adaptive=False,simulation_time=sim_time,theta=theta,info=info)
    
    if info:
        plot_sim(G_start, S, I, R)

    # adaptive case
        
    S_dyn,I_dyn,R_dyn=[],[],[]
    
    H_dyn,S_dyn_start,I_dyn_start,R_dyn_start=initialize_SIR(H,I0=I0,visualize=info)

    S_dyn.append(S_dyn_start[0])
    I_dyn.append(I_dyn_start[0])
    R_dyn.append(R_dyn_start[0])

    

    H_dyn,t_dyn,S_dyn,I_dyn,R_dyn = SIR_adaptive(H_dyn,S_dyn,I_dyn,R_dyn, beta=beta, mu=mu,adaptive=True,simulation_time=sim_time,theta=theta,info=info)
    
    if info:
        plot_sim(H_dyn, S_dyn, I_dyn, R_dyn)



    return G,H,I,I_dyn,t_stat,t_dyn

