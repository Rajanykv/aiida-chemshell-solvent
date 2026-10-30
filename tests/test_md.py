from aiida_shell import launch_shell_job
from aiida.engine import run, submit , run_get_node
from aiida.orm import load_code, SinglefileData, Dict, Bool
from aiida import load_profile 
from aiida.plugins import WorkflowFactory
from aiida_chemshell.calculations.solvation import SolventCalculation
"""Submit an MD calculation directly not using the WorkChain."""
load_profile("rajany")
builder = SolventCalculation.get_builder()
builder.code = load_code("chemsh-mamba@rajanylaptop")
#builder.structure = SinglefileData(file="h2o_bq.pun")
builder.structure = SinglefileData(file="/home/rajany/solventwork/chemsh-dial-solvent/repeat/fl-nh2.xyz")

#Note metadata is in a namespace chemsh. other inputs are not.
#builder.metadata.options.resources = {"tot_num_mpiprocs": 8,
#}

builder.metadata.options.resources = {
 "num_mpiprocs_per_machine": 4,"num_machines": 1, }
#"num_cores_per_machine": 1,
# Only set ``withmpi`` when the code itself does not declare it.
if builder.code.with_mpi is None:
    builder.metadata.options.withmpi = True

builder.mm_parameters = Dict({"theory": "DL_POLY", 'ff' : 'charmm'})

builder.md_parameters = Dict({
                'nvt'            : {'length' : 40, 'minimisation': 20},
                'npt'            : {'length' : 40, 'minimisation': 20},
                'max_ncycles'           : 10,
                'production'     : {'length' : 20},
                'padding'               : 40.0,
                'nsnapshots'            : 2,
        })

#builder.dryrunmd = True
builder.solvent_parameters = Dict( {
                'solvent_name'         :'chloroform' ,
                'box_halflength'       : 40,
        })

#node = submit(builder)
#builder.metadata.dry_run = True
results, node = run.get_node(builder)
#print("Final Energy = ", results.get("energy"))
print(f"Submitted MD Calculation PK: {node.pk}")
