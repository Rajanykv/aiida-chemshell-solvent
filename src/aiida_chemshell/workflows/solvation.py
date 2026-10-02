"""Workflows for geometry optimisation based taks."""

from aiida.common.exceptions import MissingEntryPointError
from aiida.engine import ToContext, WorkChain, if_
from aiida.orm import ArrayData, Bool, Code, Dict, Float, SinglefileData, List, FolderData
from aiida.orm import StructureData
from aiida.plugins.factories import CalculationFactory

from aiida_chemshell.calculations.solvation import SolventCalculation
from aiida_chemshell.calculations.base import ChemShellCalculation
from aiida_chemshell.workflows.isolated_atoms import IsolatedAtomicEnergiesWorkChain
from aiida.engine import ExitCode

class SolvationWorkChain(WorkChain):
    """Steps in the Solvation Work Flow."""

    @classmethod
    def define(cls, spec) -> None:
        """Define the AiiDA process specification for the WorkChain."""
        super().define(spec)

        ## Inputs ##
        spec.expose_inputs(SolventCalculation,exclude=['metadata'])
        spec.expose_inputs(SolventCalculation,include=['metadata'], namespace="chemsh")

        ## Workflow ##
        #rajany todo add full/correct workflow
        #if_(cls.not_md_dryrun)(
        #   )
        spec.outline(
            cls.validate_inputs_1,
            cls.qm_optimise,
            cls.result_opt,
            cls.charge_fit,
            cls.result_charge,
            cls.validate_inputs_2,
            cls.solvate_md,
            cls.result_md,
            #cls.setup_qmmm,
            #cls.qmmm_opt,
            #cls.result,
        )

        return

    def validate_inputs_1(self):
        """Validate the inputs provided to the WorkChain."""
        has_file = "structure" in self.inputs
        if not has_file:
            return SolventCalculation.exit_codes.ERROR_NO_INPUTS

        return None

    def validate_inputs_2(self):
        """Validate the inputs provided to the WorkChain."""
        has_box = "solvent_box" in self.inputs
        #has_ff = "ff" in self.inputs.mm_parameters.get_dict()
        #has_ff = "force_field_file" in self.inputs
        #if has_ff and "mm_parameters" not in self.inputs:
        if not has_box:
            return SolventCalculation.exit_codes.ERROR_NO_INPUTS
        return None


    def qm_optimise(self):
        """Perform the geometry optimisation."""

        inputs = self.exposed_inputs(SolventCalculation)
        inputs.update({
                  "do_init_optimise": Bool(True),
        })

        if "qm_parameters" not in self.inputs:
            inputs["qm_parameters"] = Dict(
                {
                "theory": "NWChem",
                "method": "dft",
                "functional": "B3LYP",
                "basis": "cc-pvdz",
                }
            )

        if "optimisation_parameters" not in self.inputs:
            inputs["optimisation_parameters"] = Dict({})

        if 'metadata' in self.inputs.chemsh:
            inputs["metadata"] = self.inputs.chemsh["metadata"]

        #to avoid parsing output
        inputs.update({"chargefitting_parameters": {},})

        future = self.submit(SolventCalculation, **inputs)
        future.label = SolventCalculation.default_process_label(future)
        future.description = (
            f"Geometry optimisation step from WorkChainNode pk: {self.node.pk}"
        )
        if inputs.dryrun:
            self.ctx.optimise = future
        else:
            return ToContext(optimise=future)

        del inputs
        return

    def charge_fit(self):
        """Perform the charge fitting."""
        inputs = self.exposed_inputs(SolventCalculation)

        if "qm_esp_parameters" not in self.inputs:
                inputs["qm_parameters"] = Dict(
                {
                "theory": "NWChem",
                "method": "hf",
                "functional": "B3LYP",
                "basis": "cc-pvdz",
                }
            )
        else:
                inputs["qm_parameters"] = self.inputs["qm_esp_parameters"]

        if not inputs.dryrun:
            if 'optimise' in self.ctx and self.ctx.optimise.is_finished:

                if not self.ctx.optimise.is_finished_ok:
                    return ( "Optimisation has not finished successfully")

                structure = self.ctx.optimise.outputs.optimised_structure

                inputs.update({
                    "structure"       : structure,
                })

            else:
                return("Optimisation has not finished successully")

        inputs.update({
                    "do_charge_fit" :  Bool(True),
        })

        if "chargefitting_parameters" not in self.inputs:

            inputs["chargefitting_parameters"] = Dict({
              "method" : "resp",
              "npoints" : 50,
              "type" : "shell",
              "vdw_scale" : 1.5,
              "nlayers" : 1,
              "tolerance" : 1e-12,
              })

        if 'metadata' in self.inputs.chemsh:
            inputs["metadata"] = self.inputs.chemsh["metadata"]

        #avoid parsing output at each step
        inputs["optimisation_parameters"] = {}

        #do not submit if a dryrunmd; for testing:rajany#todo remove later
        if not self.inputs.dryrunmd.value:
            future = self.submit(SolventCalculation, **inputs)
            future.label = SolventCalculation.default_process_label(future)
            future.description = (
                f"Charge Fitting Calculation Node pk: {self.node.pk}"
            )
            if inputs.dryrun:
                self.ctx.chargefit = future
            else:
                return ToContext(chargefit=future)

        del inputs
        return

    def solvate_md(self):
        """Perform the classical MD equillibration step."""

        inputs = self.exposed_inputs(SolventCalculation)
        if not inputs.dryrun:
            inputs.update({
                     "structure" : self.ctx.optimise.outputs.optimised_structure
        })
                #if qmmm_chk:
                #"qm_parameters": self.ctx.energy.inputs.qm_parameters,

        inputs.update({
                     "do_md_equillibrate" : Bool(True),
        })

        if "mm_parameters" not in self.inputs:
            mm_parameters = {"theory": "DL_POLY", 
        }
        elif "mm_parameters" in self.inputs:
            mm_parameters = self.inputs["mm_parameters"].get_dict()

        if "force_field_file" not in self.inputs and "ff" not in self.inputs.mm_parameters.get_dict():
            mm_parameters.update({ "ff" : "charmm"})

        #else:
        #rajany todo
        #generate/access prepared force field

        inputs["mm_parameters"] = Dict(mm_parameters)

        if "md_parameters" not in inputs:
            md_parameters = {
                'length_npt'            : 50,         # in fs (timestep: 2 fs)
                'length_nvt'            : 20,         # in fs (timestep: 2 fs)
                'length_production'     : 20,        # in fs (timestep: 2 fs)
                'max_ncycles'           : 20,
                'minimisation_npt'      : 5,
                'minimisation_nvt'      : 5,
                'solutes_dist'          : 3.0,
                'padding'               : 50.0,
                'nsnapshots'            : 10,
                'fixed_npt'             : '',

        }
        else:
            md_parameters = self.inputs["md_parameters"].get_dict()

        md_parameters.update({
                'driver'                :'mm_theory',
                'neutralise'            : True,
        })
        inputs["md_parameters"] = Dict(md_parameters)

        #if "qmmm_parameters" not in inputs:
        #    inputs["qmmm_parameters"] = Dict({"qm_region": []})

        if 'metadata' in self.inputs.chemsh:
            inputs["metadata"] = self.inputs.chemsh["metadata"]

        #to avoid parsing output
        inputs.update({"chargefitting_parameters": {},})

        future = self.submit(SolventCalculation, **inputs)
        future.label = SolventCalculation.default_process_label(future)
        future.description = (
            f"Solvation Calculation Node pk: {self.node.pk}"
        )
        return ToContext(md=future)

        del inputs
        return

    def setup_qmmm(self):
        """Prepare QMMM and do a dryrun."""

        inputs = self.exposed_inputs(SolventCalculation)
        if not inputs.dryrun:
            inputs.update({
                     "structure" : self.ctx.optimise.outputs.optimised_structure
        })
        if 'qmmm_parameters' in inputs:
            if 'qm_parameters' in inputs.qmmm_parameters.get_dict():        #if qmmm_chk:
                "qm_parameters" = Dict(inputs.qmmm_parameters.get_dict().qm_parameters)
            elif 'qm_parameters' in inputs:
                "qm_parameters" = inputs.qm_parameters,
            else:
                "qm_parameters" = self.ctx.optimise.inputs.qm_parameters

        inputs.update({
                     "do_qmmm" : Bool(True),
        })

        if "mm_parameters" not in self.inputs:
            mm_parameters = {"theory": "DL_POLY", 
        }
        elif "mm_parameters" in self.inputs:
            mm_parameters = self.inputs["mm_parameters"].get_dict()

        if "force_field_file" not in self.inputs and "ff" not in self.inputs.mm_parameters.get_dict():
            mm_parameters.update({ "ff" : "charmm"})

        #else:
        #rajany todo
        #generate/access prepared force field

        inputs["mm_parameters"] = Dict(mm_parameters)

        if "md_parameters" not in inputs:
            md_parameters = {
                'length_npt'            : 50,         # in fs (timestep: 2 fs)
                'length_nvt'            : 20,         # in fs (timestep: 2 fs)
                'length_production'     : 20,        # in fs (timestep: 2 fs)
                'max_ncycles'           : 20,
                'minimisation_npt'      : 5,
                'minimisation_nvt'      : 5,
                'solutes_dist'          : 3.0,
                'padding'               : 50.0,
                'nsnapshots'            : 10,
                'fixed_npt'             : '',

        }
        else:
            md_parameters = self.inputs["md_parameters"].get_dict()

        md_parameters.update({
                'driver'                :'mm_theory',
                'neutralise'            : True,
        })
        inputs["md_parameters"] = Dict(md_parameters)

        #if "qmmm_parameters" not in inputs:
        #    inputs["qmmm_parameters"] = Dict({"qm_region": []})

        if 'metadata' in self.inputs.chemsh:
            inputs["metadata"] = self.inputs.chemsh["metadata"]

        #to avoid parsing output
        inputs.update({"chargefitting_parameters": {},})

        future = self.submit(SolventCalculation, **inputs)
        future.label = SolventCalculation.default_process_label(future)
        future.description = (
            f"Solvation Calculation Node pk: {self.node.pk}"
        )
        return ToContext(md=future)

        del inputs
        return


#template for outputs:Add/remove if extra outputs to be parsed.
    def result_opt(self):
        """Extract the final workflow results."""
        if "optimised_structure" not in self.ctx.optimise.outputs:
            return(ChemShellCalculation.exit_codes.ERROR_MISSING_OPTIMISED_STRUCTURE_FILE)
        return

    def result_charge(self):
        """Extract the final workflow results."""
        if "chargefit" not in self.ctx:
            return
        if not "fitted_charges" in self.ctx.chargefit.outputs:
            if not "charges_file" in self.ctx.chargefit.outputs:
                return( ChemShellCalculation.exit_codes.ERROR_CHARGES_NOT_FOUND)
        return

    def result_md(self):
        if "do_md_equillibrate" in self.ctx.md.inputs and self.ctx.md.is_finished:

            #if SolventCalculation.FILE_SOLVATED_STRUCT in self.ctx.md.outputs.retrieved.list_object_names():
            if "solvated_structure" not in self.ctx.md.outputs:
                return SolventCalculation.exit_codes.ERROR_SOLVATED_STRUCTURE_NOT_FOUND
            if "snapshots" not in self.ctx.md.outputs:
                if self.ctx.md.inputs.dryrunmd.value:
                    pass
                else:
                    return SolventCalculation.exit_codes.ERROR_MD_SNAPSHOTS_NOT_FOUND

        else:
                return SolventCalculation.exit_codes.ERROR_MD_NOT_FINISHED

        status = self.finished_status(self.ctx.md)
        if status:
            return ExitCode(988, status)
        return

    def not_md_dryrun(self):
        return not self.inputs.dryrunmd.value

    def finished_status(self, calc_node):
        if calc_node.is_killed:
            return(f"Step {calc_node.process_label} was manually killed by user or daemon.")

        elif calc_node.is_excepted:
            return(f"Step {calc_node.process_label} crashed due to an unhandled Python exception.")

        elif calc_node.is_failed:
            return(f"Step {calc_node.process_label} finished with non-zero exit code: {calc_node.exit_status}.")

        elif "_scheduler-stderr.txt" in calc_node.outputs.retrieved.list_object_names():
            stderr = calc_node.outputs.retrieved.get_object_content("_scheduler-stderr.txt")
            if stderr:
               return(f"Step {calc_node.process_label} exited with errors; inspect _scheduler-stderr.txt/output.log.")

        elif not calc_node.is_finished_ok:
            return(f"Step {calc_node.process_label} is either still running or stopped abnormally.")
        return
