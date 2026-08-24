from aiida.engine import run, submit 
from aiida.orm import load_code, SinglefileData, Dict
from aiida import load_profile 

#load_profile("rajany")  # This is not required if running in a verdi shell environment 

builder = load_code("nwchemchemsh@rajanylaptop").get_builder() 
builder.structure = SinglefileData(file="/home/rajany/chemsh-nwchem-esp/tests/charges/nwchem/h2o.pun")
#builder.qm_parameters = Dict({"theory": "CP2K", "basis": "DZVP-MOLOPT-GTH", "functional": "lda"})
builder.qm_parameters = Dict({"theory": "NWChem", "method": "HF", "basis": "sto-3g"})
builder.chargefitting_parameters = Dict({'method':'resp', 'npoints':20, 'type':'shell', 'vdw_scale':1.5, 'nlayers':1, 'tolerance':1e-12})

computer = builder.code.computer
builder.metadata.options.withmpi = True
builder.metadata.options.resources = {
    'num_machines': 1,         # Number of compute nodes/machines
#    'num_mpiprocs_per_machine': 4     # Number of MPI processes per node
    'num_mpiprocs_per_machine': computer.get_default_mpiprocs_per_machine()
}
#print(run(builder))
results, node = run.get_node(builder)

print("Final Energy = ", results.get("energy"))
