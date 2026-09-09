"""Workflows for geometry optimisation based taks."""

from aiida.common.exceptions import MissingEntryPointError
from aiida.engine import ToContext, WorkChain
from aiida.orm import ArrayData, Bool, Code, Dict, Float, SinglefileData, List
from aiida.plugins.factories import CalculationFactory

from aiida_chemshell.calculations.solvation import SolventCalculation
from aiida_chemshell.calculations.base import ChemShellCalculation
from aiida_chemshell.workflows.isolated_atoms import IsolatedAtomicEnergiesWorkChain

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
        #rajany todo add full workflow
        spec.outline(
            cls.validate_inputs_1,
            cls.qm_optimise,
            cls.energy,
            #cls.charge_fit,
            cls.validate_inputs_2,
            cls.solvate_md,
            #cls.setup_qmmm,
            #cls.qmmm_opt,
            #cls.result,
        )

        return

    def validate_inputs_1(self):
        """Validate the inputs provided to the WorkChain."""
        has_file = "structure" in self.inputs
        if not has_file:
            return self.exit_codes.ERROR_NO_INPUTS

        return None

    def validate_inputs_2(self):
        """Validate the inputs provided to the WorkChain."""
        has_file = "structure" in self.inputs
        has_box = "solvent_box" in self.inputs
        #has_ff = "force_field_file" in self.inputs
        #if has_ff and "mm_parameters" not in self.inputs:
        if not has_file and not has_box:
            return self.exit_codes.ERROR_NO_INPUTS
        return None


    def qm_optimise(self):
        """Perform the geometry optimisation."""

        inputs = self.exposed_inputs(SolventCalculation)
        inputs.update({
                  "do_charge_fit"   : Bool(False),
                  "do_init_optimise": Bool(True),
                  "do_opt_equillibrate" : Bool(False),
                  "do_md_equillibrate"  : Bool(False),
        })

        if "qm_parameters" not in self.inputs:
            self.inputs["qm_parameters"] = Dict(
                {
                "theory": "NWChem",
                "method": "dft",
                "functional": "B3LYP",
                "basis": "cc-pvdz",
                }
            )
        inputs["qm_parameters"] = self.inputs["qm_parameters"]

        if "optimisation_parameters" not in self.inputs:
            inputs["optimisation_parameters"] = Dict({})

        if 'metadata' in self.inputs.chemsh:
            inputs["metadata"] = self.inputs.chemsh["metadata"]

        future = self.submit(SolventCalculation, **inputs)
        future.label = SolventCalculation.default_process_label(future)
        future.description = (
            f"Geometry optimisation step from WorkChainNode pk: {self.node.pk}"
        )
        if inputs.dryrun:
            self.ctx.optimise = future
        else:
            return ToContext(optimise=future)

    def energy(self):
        """Perform a single point energy calculation on the optimised structure."""

        inputs = self.exposed_inputs(SolventCalculation)

        if inputs.dryrun:
            if "qm_parameters" not in self.inputs:
                self.inputs["qm_parameters"] = Dict(
                {
                "theory": "NWChem",
                "method": "dft",
                "functional": "B3LYP",
                "basis": "cc-pvdz",
                }
            )
            qm_parameters = self.inputs["qm_parameters"]
            structure = self.inputs.structure

        elif 'optimise' in self.ctx and self.ctx.optimise.is_finished:
 
            if not self.ctx.optimise.is_finished_ok:
                return ( "Optimisation is not finished successfully")

            structure = self.ctx.optimise.outputs.optimised_structure
            qm_parameters = self.ctx.optimise.inputs.qm_parameters.get_dict()

        inputs.update({
                    "structure"       : structure,
                    "qm_parameters"   : qm_parameters,
                    "do_charge_fit"   : Bool(False),
                    "do_init_optimise": Bool(False),
                    "do_opt_equillibrate" : Bool(False),
                    "do_md_equillibrate"  : Bool(False),
        })

        if 'metadata' in self.inputs.chemsh:
            inputs["metadata"] = self.inputs.chemsh["metadata"]

        future = self.submit(SolventCalculation, **inputs)
        future.label = SolventCalculation.default_process_label(future)
        future.description = (
                f"Energy calculation on optimised structure step from WorkChainNode "
                f"pk: {self.node.pk}"
        )
        if inputs.dryrun:
            self.ctx.energy = future
        else:
            return ToContext(energy=future)
        return None

    def charge_fit(self):
        """Perform the charge fitting."""
        inputs = self.exposed_inputs(SolventCalculation)

        if self.dryrun:
            inputs.update({
                    "structure": self.inputs.structure,
                    "qm_parameters": self.inputs["qm_parameters"],
        })
        else:
            inputs.update({
                    "structure": self.ctx.energy.inputs.structure,
                    "qm_parameters": self.ctx.energy.inputs.qm_parameters,
        })
        inputs.update({
                    "do_charge_fit" :  Bool(True),
                    "do_init_optimise" : Bool(False),
                    "do_opt_equillibrate" : Bool(False),
                    "do_md_equillibrate" : Bool(False),
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

        future = self.submit(SolventCalculation, **inputs)
        future.label = SolventCalculation.default_process_label(future)
        future.description = (
            f"Charge Fitting Calculation Node pk: {self.node.pk}"
        )
        if inputs.dryrun:
            self.ctx.chargefit = future
        else:
            return ToContext(chargefit=future)

    def solvate_md(self):
        """Perform the classical MD equillibration step."""

        inputs = self.exposed_inputs(SolventCalculation)
        if self.inputs.dryrun:
            inputs.update({
                     "structure": self.inputs.structure,
        })
        else:
            inputs.update({
                     "structure": self.ctx.energy.inputs.structure,
        })
                #if qmmm_chk:
                #"qm_parameters": self.ctx.energy.inputs.qm_parameters,

        inputs.update({
                     "do_charge_fit" :  Bool(False),
                     "do_init_optimise" : Bool(False),
                     "do_opt_equillibrate" : Bool(False),
                     "do_md_equillibrate" : Bool(True),
        })

        if "mm_parameters" not in self.inputs:
            mm_parameters = {"theory": "DL_POLY", 
                              "ff"   : "charmm",
        }
        elif "mm_parameters" in self.inputs:
            mm_parameters = self.inputs["mm_parameters"].get_dict()

        #rajany note. ff is not directly passed like this
        #if "force_field_file" in self.inputs:
        #    mm_parameters.update({
        #                         'ff' : self.inputs.force_field_file.filename,
        #})
        #else:
        #rajany todo
        #generate/access force field

        inputs["mm_parameters"] = Dict(mm_parameters)

        if "md_parameters" not in inputs:
            md_parameters = {
                'driver'                :'mm_theory',
                'length_npt'            : 50,         # in fs (timestep: 2 fs)
                'length_nvt'            : 20,         # in fs (timestep: 2 fs)
                'length_production'     : 20,        # in fs (timestep: 2 fs)
                'max_ncycles'           : 20,
                'minimisation_npt'      : 5,
                'minimisation_nvt'      : 5,
                'neutralise'            : True,
                'solutes_dist'          : 3.0,
                'padding'               : 50.0,
                'nsnapshots'            : 10,
                'fixed_npt'             : '',

        }
        else:
            md_parameters = self.inputs["md_parameters"]
        #md_parameters.update()

        inputs["md_parameters"] = Dict(md_parameters)

        #if "qmmm_parameters" not in inputs:
        #    inputs["qmmm_parameters"] = Dict({"qm_region": []})

        if 'metadata' in self.inputs.chemsh:
            inputs["metadata"] = self.inputs.chemsh["metadata"]

        future = self.submit(SolventCalculation, **inputs)
        future.label = SolventCalculation.default_process_label(future)
        future.description = (
            f"Solvation Calculation Node pk: {self.node.pk}"
        )
        if inputs.dryrunmd:
            self.ctx.md = future
        else:
            return ToContext(md=future)


    def result(self):
        """Extract the final workflow results."""
        self.out(
            "Charges_file", self.ctx.chargefit.outputs.charges_file
        )
        self.out("Fitted_charges", self.ctx.chargefit.outputs.fitted_charges)
        self.out("Final_energy", self.ctx.chargefit.outputs.energy)
        return
