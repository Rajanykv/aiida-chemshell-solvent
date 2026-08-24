from aiida.engine import run, submit , run_get_node
from aiida.orm import load_code, SinglefileData, Dict
from aiida import load_profile 
from aiida.plugins import WorkflowFactory

#load_profile("rajany")  # This is not required if running in a verdi shell environment 
#builder = load_code("chemsh").get_builder() 
#print(builder)
#builder.structure = SinglefileData(file="/home/jovyan/work/h2o_bq.pun")
#builder.qm_parameters = Dict({"theory": "NWChem", "method": "HF", "basis": "sto-3g"})
#builder.chargefitting_parameters = Dict({'method':'resp', 'npoints':20, 'type':'shell', 'vdw_scale':1.5, 'nlayers':1, 'tolerance':1e-12})
##computer = builder.code.computer
#builder.metadata.options.withmpi = True
#builder.metadata.options.resources = {
#    'num_machines': 1,         # Number of compute nodes/machines
#    'num_mpiprocs_per_machine': 4     # Number of MPI processes per node
#    'num_mpiprocs_per_machine': computer.get_default_mpiprocs_per_machine()
#}
#print(run(builder))
#results, node = run.get_node(builder)

myworkchain = WorkflowFactory("chemshell.solvation")
builder = myworkchain.get_builder()
print(builder)
#print(builder.chemsh)
code = load_code("chemsh")
f= SinglefileData(file="/home/jovyan/work/h2o_bq.pun")

chemsh_inp =  { 'structure' : f, 
               'code' : code,
               'metadata':{
               'options':{ 'resources' : { 'num_machines' :1 , 'num_mpiprocs_per_machine' : 4 } ,
                           'withmpi' : True ,
                         },
               },
              'qm_parameters' : {"theory": "NWChem", "method": "HF", "basis": "sto-3g"},
              'chargefitting_parameters' : {'method':'resp', 'npoints':20, 'type':'shell', 
                                       'vdw_scale':1.5, 'nlayers':1, 'tolerance':1e-12}
             }
inputs = { 'chemsh' : chemsh_inp }

results, node = run_get_node(myworkchain, **inputs)
print("Final Energy = ", results.get("Final_energy"))
print("Charges fitted = ", results.get("Fitted_charges").get_list())
