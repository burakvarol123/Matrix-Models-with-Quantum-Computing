#This script is used to collect the useful ansatzes that worked on the way

from qiskit import *
from qiskit.circuit import Parameter

def ansatz_review_exact(N, depth):  # i only use ry gates now, exactly like in the paper
    """
    Ansatz designed in the way of the paper Variational QBMs. Close but just an approximate.
    """
    circuit = QuantumCircuit(N)
    thetas = []
    for l in range(N * (depth + 1)):
        thetas.append(Parameter('θ_' + str(chr(l))))
    counter = 0
    for i in range(depth):
        for l in range(N):
            if l + counter < len(thetas):
                circuit.ry(thetas[l + counter + N], l )
        counter = counter + N
        for p in range(N):
            circuit.cx( (p+N-1)%N,(p + N) % N)
        circuit.barrier()
    for j in range(N):
        circuit.ry(thetas[j], j)
    for k in range(N):
        if k + int(N / 2) < N:
            circuit.cx(k, (k + int(N / 2)))
    return circuit
    
def ansatz_review_exact_non(N, depth):
    """
    Ansatz with Hadamart gates in the beginning and end of the depth layers. Was usefull to see that the H gates "activates" the layers
    
    """
    circuit = QuantumCircuit(N)
    thetas = []
    for l in range(N * (depth + 1)):
        thetas.append(Parameter('θ_' + str(chr(l))))
    counter = 0
    for i in range(depth):
        for l in range(N):
            if l + counter < len(thetas):
                circuit.h(l)
                circuit.ry(thetas[l + counter + N], l )
                
        counter = counter + N
        for p in range(N):
            circuit.cx( (p+N-1)%N,(p + N) % N)
            circuit.h(p)
        circuit.barrier()
    for j in range(N):
        circuit.ry(thetas[j], j)
    for k in range(N):
        if k + int(N / 2) < N:
            circuit.cx(k, (k + int(N / 2)))
    return circuit

def ansatz_review_exact_full(N, depth): 
    """
    FUll Entanglement Ansatz with CNOTs. Does not work that well
    """
    circuit = QuantumCircuit(N)
    thetas = []
    for l in range(N * (depth + 1)):
        thetas.append(Parameter('θ_' + str(chr(l))))
    counter = 0
    for i in range(depth):
        for l in range(N):
            if l + counter < len(thetas):
                circuit.h(l)
                circuit.ry(thetas[l + counter + N], l )       
        counter = counter + N
        for p in range(N):
            for j in range(N):
                if p != p+j and p+j <= N-1:
                    circuit.cx( p,p +j)
            if p == N-1:
                circuit.cx(p,0)
            circuit.h(p) 
        circuit.barrier()
    for j in range(N):
        circuit.ry(thetas[j], j)
    for k in range(N):
        if k + int(N / 2) < N:
            circuit.cx(k, (k + int(N / 2)))
    return circuit

def ansatz_review_exact_full_cry(N, depth): 
    """
    Full Entanglement Ansatz with Controlled Ry Gates. The best by far, but takes too much time. Can be massively optimised.
    """
    circuit = QuantumCircuit(N)
    thetas = []
    for l in range(N * (depth + 1)+ int(depth*N*(N-1)/2)+100):
        thetas.append(Parameter('θ_' + str(chr(l))))
    counter = 0
    for i in range(depth):
        for l in range(N):
            if l + counter < len(thetas):
                circuit.h(l)
                circuit.ry(thetas[l + counter + N], l ) 
        counter_3 = 0
        for p in range(N):
            counter_2 = 0
            for j in range(N):
                if p != p+j and p+j <= N-1:
                    circuit.cry(thetas[j+p+counter+2*N+counter_3-1],p,p +j)
                    counter_2 += 1
            counter_3 += counter_2-1
            if p == N-1:
                circuit.cx(p,0)
            circuit.h(p)
        counter = counter + 2*N +counter_3
        circuit.barrier()
    for j in range(N):
        circuit.ry(thetas[j], j)
    for k in range(N):
        if k + int(N / 2) < N:
            circuit.cx(k, (k + int(N / 2)))
    print(counter_3)
    return circuit