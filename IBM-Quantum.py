#!/bin/python3
import numpy as np
from qiskit.circuit.library import IQP
from qiskit import QuantumCircuit, transpile
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
from qiskit.quantum_info import random_hermitian
from qiskit_ibm_runtime import QiskitRuntimeService, Sampler, Options, Session

nbQubit = 16
repeats = 8
max_time_seconds = 599

service = QiskitRuntimeService(channel="ibm_quantum", token="YOUR-API-KEY-HERE")
backendProperties = service.get_backend("ibm_brisbane").properties() #ibmq_qasm_simulator
backend = service.get_backend("ibm_brisbane") #"//least_busy(simulator=True, operational=True) ibm_brisbane
options = Options(resilience_level=1, optimization_level=3)
options.transpilation.skip_transpilation=True

rng = np.random.default_rng()
mats = [np.real(random_hermitian(nbQubit, seed=rng)) for _ in range(repeats)]
circuits = [IQP(mat) for mat in mats]
for circuit in circuits:
    circuit.measure_all()
  
pm = generate_preset_pass_manager(0,backend=backend)
isa_circuits = pm.run(circuits)

compiled_circuits = transpile(isa_circuits, backend=backend, backend_properties=backendProperties)

with Session(service, backend=backend) as session:
    sampler = Sampler(backend=backend, options=options, session=session)
    while True:
        f = open("Qubits.txt", 'a')
        job = sampler.run(compiled_circuits, shots=1024)
        print(f"Job ID: {job.job_id()}")
        result = job.result()
        for job in result.quasi_dists:
            for result in job:
                f.write(str(job[result]) + "\n")
        f.close()
