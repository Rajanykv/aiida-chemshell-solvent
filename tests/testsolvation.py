
from aiida.engine import run, submit , run_get_node
from aiida.orm import load_code, SinglefileData, Dict
from aiida import load_profile 
from aiida.plugins import WorkflowFactory
"""Submit the Solvation WorkChain."""
builder = WorkflowFactory("chemshell.solvation").get_builder()  # pyright: ignore[reportFunctionMemberAccess]
builder.chemsh.code = load_code("nwchemchemsh")
#builder.chemsh.structure = SinglefileData(file="/home/jovyan/work/h2o_bq.pun")
builder.chemsh.structure = SinglefileData(file="/home/rajany/aiida-chemshell-solv22/tests/h2o_bq.pun")
builder.chemsh.qm_parameters = Dict(
    {
        "theory": "NWChem",
        "method": "dft",
        #"functional": 
        "basis": "sto-3g",
    }
)

builder.chemsh.metadata.options.resources = {
    "num_mpiprocs_per_machine": 4,
    "num_cores_per_machine": 1,
    "num_machines": 1,
    "tot_num_mpiprocs": 4
}
builder.chemsh.chargefitting_parameters = {
    'method':'resp', 
    'npoints': 50,
    'type':'shell',
    'vdw_scale':1.5, 
    'nlayers':1, 
    'nlayers':1, 
    'tolerance':1e-12
    }
     # Only set ``withmpi`` when the code itself does not declare it.
if builder.chemsh.code.with_mpi is None:
    builder.chemsh.metadata.options.withmpi = True
node = submit(builder)
print(f"Submitted WorkChain PK: {node.pk}")
