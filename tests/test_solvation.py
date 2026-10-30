
from aiida.engine import run, submit , run_get_node
from aiida.orm import load_code, SinglefileData, Dict
from aiida import load_profile 
from aiida.plugins import WorkflowFactory
"""Submit the Solvation WorkChain."""
builder = WorkflowFactory("chemshell.solvation").get_builder()
builder.code = load_code("chemsh-mamba@rajanylaptop")
#builder.structure = SinglefileData(file="/home/jovyan/work/h2o_bq.pun")
builder.structure = SinglefileData(file="/home/rajany/solventwork/chemsh-dial-solvent/repeat/fl-nh2.xyz")
builder.qm_parameters = Dict(
    {
        "theory": "NWChem",
        "method": "dft",
        "functional": 'blyp',
        "basis": "3-21g",
    }
)

#Note metadata is in a namespace chemsh. other inputs are not.
builder.chemsh.metadata.options.resources = {
"num_cores_per_machine": 1,
 "num_mpiprocs_per_machine": 4,"num_machines": 1, }
#"tot_num_mpiprocs": 4,}
# Only set ``withmpi`` when the code itself does not declare it.
if builder.code.with_mpi is None:
    builder.metadata.options.withmpi = True

#builder.chargefitting_parameters = Dict({
#    'method':'resp', 
#    'npoints': 50,
#    'type':'shell',
#    'vdw_scale':1.5, 
#    'nlayers':1, 
#    'nlayers':1, 
#    'tolerance':1e-12
#    })
builder.dryrun_opt = True
builder.dryrun_charge = True
#rajany note: validator at base needs mm_parameters to be provided explicitly before it can be set any default value by process,
builder.mm_parameters = Dict({"theory": "DL_POLY", 'ff' : 'charmm'})
builder.md_parameters = Dict({
                'nvt'            : {'length' : 40, 'minimisation': 20},
                'npt'            : {'length' : 40, 'minimisation': 20},
                'max_ncycles'           : 10,
                'production'     : {'length' : 20},
                'padding'               : 40.0,
                'nsnapshots'            : 2,
        })

builder.solvent_parameters = Dict( {
                'solvent_name'         :'chloroform' ,
                'box_halflength'       : 40,
        })

#node = submit(builder)
#builder.metadata.dry_run = True
results, node = run.get_node(builder)
#print("Final Energy = ", results.get("energy"))
print(f"Submitted WorkChain PK: {node.pk}")
