import stim
import numpy as np

def genCliffordEigenStates(c,N_QUBITS):
  paulis = ["I"] * N_QUBITS
  initial = ["I"] * N_QUBITS
  s_flippers = ["I"] * N_QUBITS
  measurement = ["I"] * N_QUBITS
  # this looks dreadful but basically provides operators for each clifford to:
  #   prepare each eigenstate
  #   flip between eigenstates
  #   correct for measurement in the computational basis
  pauli_map = {
      "H":"I",
      "S":"Z",
      "I":"I",
      "CX":"ZX"
  }
  initial_map = {
      "H":"H",
      "S":"Z",
      "I":"I",
      "CX":"IH"
  }
  s_map = {
      "H":"Z",
      "S":"X",
      "I":"X",
      "CX":"XZ"
  }
  m_map= {
      "H":"I",
      "S":"I",
      "I":"I",
      "CX":"IH"
  }

  maps = [pauli_map,initial_map,s_map,m_map]
  lists = [paulis,initial,s_flippers,measurement]

  for gate in c:

    for i in range(len(maps)):

      split_gate = str(gate).split(" ")
      current_map = maps[i]
      current_list = lists[i]
      unitary = current_map[split_gate[0]]

      if(len(unitary) > 1):
        current_list[int(split_gate[1])] = unitary[0]
        current_list[int(split_gate[2])] = unitary[1]
      else:
        for i in range(len(split_gate)-1):
          current_list[int(split_gate[i + 1])] = unitary

  return paulis,initial,s_flippers,measurement

def getPCAPhase(circuit, pauli_string):
  c = stim.Tableau.from_circuit(circuit)
  pct = c(pauli_string)
  return pct


def generate_circuit(N_QUBITS,s,clifford_circuit):
  prep_circuit = stim.Circuit()
  paulis,initial,s_flippers,m_prep = genCliffordEigenStates(clifford_circuit,N_QUBITS)

  # Prepare Pauli Eigenstates
  for i in range(N_QUBITS):
    prep_circuit.append(initial[i],i)

  # Choose Specific Eigenstate from s
  for i in range(N_QUBITS):
    if s[i]:
      prep_circuit.append(s_flippers[i],i)

  pauli_string = stim.PauliString(paulis)
  phase = getPCAPhase(clifford_circuit,pauli_string).sign

  # sandwich state prep circuit with clifford circuit
  c = prep_circuit.__iadd__(clifford_circuit)

  for i in range(N_QUBITS):
    c.append(m_prep[i],i)

  c.append("M", range(N_QUBITS))

  return c,phase

def createBitString(i,N_QUBITS):
  s_bitstring = bin(i)

  # create an s bitstring
  s = []
  for j in range(len(s_bitstring)-1,1,-1):
    s.append(bool(int(s_bitstring[j])))
  for j in range(N_QUBITS - len(s_bitstring)+2):
    s.append(False)
  return s

def fullSampler(clifford_circuit):
  N_SAMPLES = 2045
  N_QUBITS = clifford_circuit.num_qubits
  batch = []

  for i in range( 2 ** N_QUBITS): # for every bit string s
    s = createBitString(i,N_QUBITS)
    c, phase = generate_circuit(N_QUBITS,s,clifford_circuit)

    samples = c.compile_sampler().sample(N_SAMPLES)
    sample = np.mean(samples, axis=0).round(3)
    # flip every sample
    for j in range(N_QUBITS):
      if(s[j]):
        sample[j] = 1- sample[j]

    batch.append(sample)

  batch = np.array(batch)
  totals = [0] * N_QUBITS

  for j in range(len(batch[0])):
    for i in range(len(batch)):
        totals[j] += 1-batch[i][j]

  for i in range(len(totals)):
    totals[i] = totals[i] / (2**N_QUBITS)

  # print(batch)
  # print(totals)
  circuit_eigenvalue = (np.sum(totals) * phase / N_QUBITS).round(3)

  return circuit_eigenvalue
