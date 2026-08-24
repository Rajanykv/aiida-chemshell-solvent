"""Workflows for geometry optimisation based taks."""

from aiida.common.exceptions import MissingEntryPointError
from aiida.engine import ToContext, WorkChain
from aiida.orm import ArrayData, Bool, Code, Dict, Float, SinglefileData, List
from aiida.plugins.factories import CalculationFactory

from aiida_chemshell.calculations.solvation import SolventCalculation
from aiida_chemshell.workflows.isolated_atoms import IsolatedAtomicEnergiesWorkChain

class SolvationWorkChain(WorkChain):
    """Steps in the Solvation Work Flow."""

    @classmethod
    def define(cls, spec) -> None:
        """Define the AiiDA process specification for the WorkChain."""
        super().define(spec)

        ## Inputs ##
        spec.expose_inputs(SolventCalculation, namespace="chemsh")

        spec.output(
            "Final_energy",
            valid_type=Float,
            required=True,
            help="The final energy for the structure.",
        )
        spec.output(
            "Charges_file",
            valid_type=SinglefileData,
            required=True,
            help="The file containing the fitted  charges of atoms.",
        )
        spec.output(
            "Fitted_charges",
            valid_type=List,
            required=False,
            help="The calculated fitted charges for the structure",
           )
        #rajany todo examine later correct inputs
        spec.exit_code(
            308,
            "ERROR_NO_INPUTS",
            message=(
                "Required Inputs are not provided."
            ),
        )

        ## Workflow ##
        #rajany todo add full workflow
        spec.outline(
            cls.validate_inputs_1,
            cls.qm_optimise,
            cls.energy,
            cls.charge_fit,
            cls.validate_inputs_2,
            #cls.setup_qmmm,
            #cls.classical_md,
            cls.result,
        )

        return

    def validate_inputs_1(self):
        """Validate the inputs provided to the WorkChain."""
        has_file = "structure" in self.inputs.chemsh
        if not has_file:
            return self.exit_codes.ERROR_NO_INPUTS
        return None

    def validate_inputs_2(self):
        """Validate the inputs provided to the WorkChain."""
        has_file = "structure" in self.inputs.chemsh
        has_box = "solvent_box" in self.inputs.chemsh
        #if not has_files and not has_box:
        #if not has_box:
        #    return self.exit_codes.ERROR_NO_INPUTS
        return None

    def qm_optimise(self):
        """Perform the geometry optimisation."""

        inputs = self.exposed_inputs(SolventCalculation, namespace="chemsh")

        inputs["do_charge_fit"] =  Bool(False)
        inputs["do_init_optimise"] = Bool(True)
        inputs["do_opt_equillibrate"] = Bool(False)
        inputs["do_md_equillibrate"] = Bool(False)
        if "qm_parameters" not in inputs:
            inputs["qm_parameters"] = Dict(
                {
                "theory": "NWChem",
                "method": "dft",
                "functional": "B3LYP",
                "basis": "cc-pvdz",
                }
            )
        if "optimisation_parameters" not in inputs:
            inputs["optimisation_parameters"] = Dict({})

        future = self.submit(SolventCalculation, **inputs)
        future.label = SolventCalculation.default_process_label(future)
        future.description = (
            f"Geometry optimisation step from WorkChainNode pk: {self.node.pk}"
        )
        return ToContext(optimise=future)

    def energy(self):
        """Perform a single point energy calculation on the optimised structure."""

        if 'optimise' in self.ctx and self.ctx.optimise.is_finished_ok:
            structure=self.ctx.optimise.outputs.optimised_structure
            qm_parameters=self.ctx.optimise.inputs.qm_parameters,
        else:
            structure=self.inputs.chemsh.structure
            qm_parameters=self.inputs.chemsh.qm_parameters
        inputs = {
                "code": self.exposed_inputs(SolventCalculation, namespace="chemsh")[
                    "code"
                ],
                "metadata": self.exposed_inputs(
                    SolventCalculation, namespace="chemsh"
                )["metadata"],
                "structure": structure,
                "qm_parameters": qm_parameters,
                "do_charge_fit" :  Bool(False),
                "do_init_optimise" : Bool(False),
                "do_opt_equillibrate" : Bool(False),
                "do_md_equillibrate" : Bool(False),
        }

        future = self.submit(SolventCalculation, **inputs)
        future.label = SolventCalculation.default_process_label(future)
        future.description = (
                f"Energy calculation on optimised structure step from WorkChainNode "
                f"pk: {self.node.pk}"
            )
        return ToContext(energy=future)
        return None

    def charge_fit(self):
        """Perform the charge fitting."""

        inp = self.exposed_inputs(SolventCalculation, namespace="chemsh")
        inputs = {
                "code": inp["code"],
                "metadata": inp["metadata"],
                "structure": self.ctx.energy.inputs.structure,
                "qm_parameters": self.ctx.energy.inputs.qm_parameters,
                "do_charge_fit" :  Bool(True),
                "do_init_optimise" : Bool(False),
                "do_opt_equillibrate" : Bool(False),
                "do_md_equillibrate" : Bool(False),
        }

        if "chargefitting_parameters" not in inp:
            inputs["chargefitting_parameters"] = Dict({
              "method" : "resp",
              "npoints" : 50,
              "type" : "shell",
              "vdw_scale" : 1.5,
              "nlayers" : 1,
              "tolerance" : 1e-12,
              })
        else:

            inputs["chargefitting_parameters"] = inp["chargefitting_parameters"]

        future = self.submit(SolventCalculation, **inputs)
        future.label = SolventCalculation.default_process_label(future)
        future.description = (
            f"Charge Fitting Calculation Node pk: {self.node.pk}"
        )
        return ToContext(chargefit=future)

    def classical_md(self):
        """Perform the classical MD equillibration step."""

        inp = self.exposed_inputs(SolventCalculation, namespace="chemsh")
        inputs = {
                "code": inp["code"],
                "metadata": inp["metadata"],
                "structure": inp["solvent_box"],
                "do_charge_fit" :  Bool(False),
                "do_init_optimise" : Bool(False),
                "do_opt_equillibrate" : Bool(False),
                "do_md_equillibrate" : Bool(True),
        }

        if "md_parameters" not in inp:
            inputs["md_parameters"] = Dict({
              "driver" : "dl_poly",
              "ensemble" : "nvt",
              "nsteps" : 2000,
              "nsnapshots" : 10,
              "temperature" : 300.0,
              "timestep" : 0.5,
              })
        else:

            inputs["md_parameters"] = inp["md_parameters"]

        future = self.submit(SolventCalculation, **inputs)
        future.label = SolventCalculation.default_process_label(future)
        future.description = (
            f"Charge Fitting Calculation Node pk: {self.node.pk}"
        )
        return ToContext(md=future)


    def result(self):
        """Extract the final workflow results."""
        self.out(
            "Charges_file", self.ctx.chargefit.outputs.charges_file
        )
        self.out("Fitted_charges", self.ctx.chargefit.outputs.fitted_charges)
        self.out("Final_energy", self.ctx.chargefit.outputs.energy)
        return
