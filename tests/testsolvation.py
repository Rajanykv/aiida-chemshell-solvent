
from aiida.engine import run, submit , run_get_node
from aiida.orm import load_code, SinglefileData, Dict
from aiida import load_profile 
from aiida.plugins import WorkflowFactory
"""Submit the Solvation WorkChain."""
builder = WorkflowFactory("chemshell.solvation").get_builder()  # pyright: ignore[reportFunctionMemberAccess]
builder.code = load_code("chemshdlpolysol@rajanylaptop")
#builder.structure = SinglefileData(file="/home/jovyan/work/h2o_bq.pun")
builder.structure = SinglefileData(file="/home/rajany/solventwork/aiida-chemshell-solvent/tests/h2o_bq.pun")
builder.qm_parameters = Dict(
    {
        "theory": "NWChem",
        "method": "dft",
        #"functional": 
        "basis": "sto-3g",
    }
)

#Note metadata is in a namespace chemsh. other inputs are not.
builder.chemsh.metadata.options.resources = {
                                              "tot_num_mpiprocs": 4,}

# "num_mpiprocs_per_machine": 2,"num_machines": 1, }
#"num_cores_per_machine": 1,
builder.chargefitting_parameters = Dict({
    'method':'resp', 
    'npoints': 50,
    'type':'shell',
    'vdw_scale':1.5, 
    'nlayers':1, 
    'nlayers':1, 
    'tolerance':1e-12
    })
     # Only set ``withmpi`` when the code itself does not declare it.
#if builder.code.with_mpi is None:
#    builder.metadata.options.withmpi = True
#node = submit(builder)
#builder.metadata.dry_run = True
results, node = run.get_node(builder)
#print("Final Energy = ", results.get("energy"))
print(f"Submitted WorkChain PK: {node.pk}")
