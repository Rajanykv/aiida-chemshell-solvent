"""Workflows for geometry optimisation based taks."""

from aiida.engine import ToContext, WorkChain
from aiida.orm import Dict, Float, SinglefileData, List

from aiida_chemshell.calculations.base import ChemShellCalculation

class ChargeFitWorkChain(WorkChain):
    """Charge fitting calculation."""

    @classmethod
    def define(cls, spec) -> None:
        """Define the AiiDA process specification for the WorkChain."""
        super().define(spec)

        ## Inputs ##
        spec.expose_inputs(ChemShellCalculation, namespace="chemsh")

        ## Outputs ##
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

        ## Workflow ##
        spec.outline(
            #cls.validate_inputs,
            cls.charge_fit,
            cls.result,
        )
        return

    def charge_fit(self):
        """Perform the charge fitting."""
        inputs = self.exposed_inputs(ChemShellCalculation, namespace="chemsh")
        if "qm_parameters" not in inputs:
            inputs["qm_parameters"] = Dict(
                {
                    "theory": "NWChem",
                    "method": "dft",
                    "functional": "B3LYP",
                    "basis": "cc-pvdz",
                }
            )
        if "chargefitting_parameters" not in inputs:
            inputs["chargefitting_parameters"] = Dict({
              "method" : "resp",
              "npoints" : 50,
              "type" : "shell",
              "vdw_scale" : 1.5,
              "nlayers" : 1,
              "tolerance" : 1e-12,
              })

        future = self.submit(ChemShellCalculation, **inputs)
        future.label = ChemShellCalculation.default_process_label(future)
        future.description = (
            f"Charge Fitting Calculation Node pk: {self.node.pk}"
        )
        return ToContext(chargefit=future)

    def result(self):
        """Extract the final workflow results."""
        self.out(
            "Charges_file", self.ctx.chargefit.outputs.charges_file
        )
        self.out("Fitted_charges", self.ctx.chargefit.outputs.fitted_charges)
        self.out("Final_energy", self.ctx.chargefit.outputs.energy)
        return
