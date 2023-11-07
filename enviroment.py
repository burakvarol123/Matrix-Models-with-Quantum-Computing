import gymnasium as gym
import numpy as np
from gymnasium import spaces
from qiskit import QuantumCircuit
from qiskit.circuit import Parameter
from qiskit.circuit.library import RXGate, RYGate, RZGate
from qiskit.primitives import Estimator
from scipy.optimize import minimize


class Circuit(gym.Env):

    def __init__(self, numqubits, numgates, hamiltonian):
        self.numqubits = numqubits
        self.numgates = numgates
        self.num_rot_gate = 0
        self.cur_numgates = 0
        self.hamiltonian = hamiltonian
        self.fake_min_energy = -5
        self.accuracy = 7
        self.action_space = spaces.Discrete(numqubits * (numqubits + 2))
        self.observation_space = spaces.Box(low=0, high=numqubits, shape=(4*numgates,))

        self.action_to_gate = []

        for i in range(self.numqubits):
            for j in range(1, 4):
                self.action_to_gate.append([i, j])

        for i in range(self.numqubits):
            for k in range(self.numqubits):
                if i != k:
                    self.action_to_gate.append([i, k])

    def __get_obs__(self):
        return {"circuit": np.array(self.circuit_state).flatten()}

    def reset(self, seed=None, options=None):
        super().reset(seed=seed)

        self.circuit_state = [
            [4] * self.numgates,
            [0] * self.numgates,
            [4] * self.numgates,
            [0] * self.numgates
        ]
        self.circuit_state = np.array(self.circuit_state)
        self.cur_numgates = 0
        self.num_rot_gate = 0

        self.circuit = QuantumCircuit(self.numqubits)

        self.backward_energy = 2

        observation = self.__get_obs__()

        return observation, {}

    def R_gate(self, axis):
        if axis == 1:
            self.num_rot_gate += 1
            return RXGate(Parameter('θ_' + str(self.num_rot_gate)))
        elif axis == 2:
            self.num_rot_gate += 1
            return RYGate(Parameter('θ_' + str(self.num_rot_gate)))
        elif axis == 3:
            self.num_rot_gate += 1
            return RZGate(Parameter('θ_' + str(self.num_rot_gate)))
        else:
            print("Wrong gate")
            return 1

    def step(self, action):
        gate = self.action_to_gate[action]
        rot_added = False

        if self.numqubits * 3 > action:
            rot_added = not rot_added
            self.circuit_state[2:, self.cur_numgates] = gate
            self.circuit.append(self.R_gate(gate[1]), [gate[0]])
        else:
            self.circuit_state[:2, self.cur_numgates] = gate
            self.circuit.cx(gate[0], gate[1])

        self.circuit.barrier()
        self.cur_numgates += 1

        if rot_added:
            vqe = VQE(self.circuit, self.hamiltonian, np.random.uniform(low=0, high=np.pi, size=self.num_rot_gate))
            optimized = vqe.optimization()

            if optimized.fun < self.accuracy:
                reward = 5
            elif self.cur_numgates >= self.numgates and optimized.fun >= self.accuracy:
                reward = -5
            else:
                reward = max((self.backward_energy - optimized.fun) / (self.backward_energy - self.fake_min_energy), -1)

            self.backward_energy = optimized.fun
        else:
            reward = 0

        terminated = np.equal(self.cur_numgates, self.numgates)
        observation = self.__get_obs__()

        return observation, reward, terminated, False, {}


class VQE:

    def __init__(self, ansatz, hamiltonian, param):
        self.ansatz = ansatz
        self.hamiltonian = hamiltonian
        self.param = param

    def cost_function(self, param):
        bound = self.ansatz.bind_parameters({self.ansatz.parameters[i]: param[i] for i in range(len(param))})
        estimator = Estimator()
        expectation_value = estimator.run(bound, self.hamiltonian).result().values[0]
        return expectation_value

    def optimization(self):
        optimized = minimize(fun=self.cost_function, x0=self.param, method="L-BFGS-B")
        return optimized
