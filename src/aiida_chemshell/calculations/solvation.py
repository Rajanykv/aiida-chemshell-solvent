"""Core ChemShell calculations module."""

from aiida.common import CalcInfo, CodeInfo
from aiida.common.folders import Folder
from aiida.engine import CalcJob, CalcJobProcessSpec, PortNamespace
from aiida.orm import (
    ArrayData,
    Dict,
    Float,
    Int,
    SinglefileData,
    StructureData,
    TrajectoryData,
    List,
    Bool,
    FolderData,
)

from aiida_chemshell.units import UnitsConverter
from aiida_chemshell.utils import ChemShellMMTheory, ChemShellQMTheory
from aiida_chemshell.calculations.base import ChemShellCalculation


class SolventCalculation(ChemShellCalculation):
    """
    AiiDA calculation plugin wrapper for Solvation calculations to avoid clutter in base.

    Currently supports the following tasks:
      - Single point energy
      - Geometry optimisation

    """
    FOLDER_SNAPSHOTS = "_snapshots"

    @classmethod
    def define(cls, spec: CalcJobProcessSpec) -> None:
        """
        Define the inputs, outputs and metadata of the ChemShell calculation.

        Parameters
        ----------
        spec : CalcJobProcessSpec
            The AiiDA Process specification object for the job.
        """
        super().define(spec)

        spec.input(
            "do_sp",
             valid_type = Bool,
             default=lambda: Bool(False),
             required = True,
             help = "Whether to do a single point energy calculation. (default False)"
        )
        spec.input(
            "do_init_optimise",
             valid_type = Bool,
             default=lambda: Bool(False),
             required = True,
             help = "Whether to do optimisation on the initial structure. (default False)"
        )
        spec.input(
            "do_charge_fit",
             valid_type = Bool,
             default=lambda: Bool(False),
             required = True,
             help = "Whether to do a QM RESP charge fitting(default False)"
        )
        spec.input(
            "do_md_equillibrate",
             valid_type = Bool,
             default=lambda: Bool(False),
             required = True,
             help = "Whether to do an MD equillibration step (default False)"
        )
        spec.input(
            "do_opt_equillibrate",
             valid_type = Bool,
             default=lambda: Bool(False),
             required = True,
             help = "Whether to do an optimisation instead of MD equillibration(default False)"
        )
        spec.input(
            "dryrunmd",
             valid_type = Bool,
             default=lambda: Bool(False),
             required = True,
             help = "Whether to do a dry run with the MD parameters(default False)"
        )
        spec.input(
            "dryrun",
             valid_type = Bool,
             default=lambda: Bool(False),
             required = True,
             help = "Whether to do a dry run of all the steps in solvation(default False)"
        )
        spec.input(
            "solvent_box",
            valid_type=SinglefileData,
            validator=cls.validate_inputs_solvent,
            required=False,
            help=(
                "The solvent boxes input structure to be used for the ChemShell solvation"
                "Choose from the list in format '.pqr' or '.pdb'"
            ),
        )
        spec.input(
            "md_parameters",
            valid_type=Dict,
            required=False,
            validator=cls.validate_md_parameters,
            help="A dictionary of parameters for the ChemShell MD Solvation.",
        )

        #rajany check metadata is not inherited
        spec.inputs["metadata"]["options"]["resources"].default = {
            "num_machines": 1,
            "num_mpiprocs_per_machine": 4,
        }
        spec.inputs["metadata"]["options"]["parser_name"].default = "chemshell"

        #rajany todo examine later correct specs
        spec.output(
            "snapshots",
            valid_type=FolderData,
            required=False,
            help=(
                "Snapshots from an MD simulation."
            ),
        )

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

        spec.exit_code(
            307,
            "ERROR_NO_INPUTS",
            message=(
                "Required Inputs are not provided."
            ),
        )
        spec.exit_code(
            308,
            "ERROR_MD_NOT_FINISHED",
            message=(
                "MD snapshots not found. MD run has not completed."
            ),
        )
        spec.exit_code(
            306,
            "ERROR_MD_NOT_FINISHED",
            message=(
                "Failed to complete the MD equillibration."
            ),
        )


        return

    @classmethod
    def validate_inputs_solvent(cls, value: SinglefileData | None, _) -> str | None:
        """Validate the inputs provided to the WorkChain."""
        if isinstance(value, SinglefileData):
            if value.filename[-4:] not in [".xyz", ".pun", ".pqr", ".pdb"]:
                if value.filename[-6:] != ".cjson":
                    return (
                        "Structure file must be either an '.xyz', '.pun', '.pqr', '.pdb' or "
                        "'.cjson' formatted structure file."
                    )

        return None


    #rajany todo-retain only necessary
    @classmethod
    def get_valid_md_parameters(cls) -> dict[str:type]:
        """
        Return a tuple of valid parameter keys for the Chemshell MD calculation.

        Returns
        -------
        validKeys : dict[str: type]
            A tuple of valid parameter keys for the ChemShell Charge fitting calculation.
        """
        #rajany todo,
        return {
                'driver'                : str,
                'length_npt'            : int,         # in fs (timestep: 2 fs)
                'length_nvt'            : int,         # in fs (timestep: 2 fs)
                'length_production'     : int,        # in fs (timestep: 2 fs)
                'max_ncycles'           : int,
                'minimisation_npt'      : int,
                'minimisation_nvt'      : int,
                'neutralise'            : bool,
                'solute'                : str,
                'solutes_dist'          : float,
                'solvent'               : str,
                'padding'               : float,
                'nsnapshots'            : int,
                'fixed_npt'             : str,

        }
    @classmethod
    def validate_md_parameters(cls, value: Dict | None, _) -> str | None:
        """
        Validate the MD parameters.

        Parameters
        ----------
        value : Dict | None
            A dictionary of parameters for the ChemShell Charge fitting calculation.
            If None, no validation is performed.

        Returns
        -------
        str | None
            Returns None if the parameters are valid, otherwise returns an error
            message string.
        """
        valid_keys = cls.get_valid_md_parameters()

        # Check for valid parameter keys
        invalid_keys = set(value.keys()).difference(set(valid_keys.keys()))
        if invalid_keys:
            return (
                "The following parameter keys are invalid: "
                f"{', '.join(invalid_keys):s}. Valid keys are: "
                f"{', '.join(valid_keys.keys()):s}"
            )

        # Check for valid parameter types
        for key, val in value.items():
            if not isinstance(val, valid_keys[key]):
                return (
                    f"The parameter '{key:s}' must be of type "
                    f"{valid_keys[key].__name__:s}."
                )

        return None


    def _build_process_label(self) -> str:
        """
        AiiDA Process label definition.

        Defines the process label to be associated with the created ProcessNode
        stored in the AiiDA database.

        Returns
        -------   
        str
            The process label based on what inputs have been provided.
        """
        return SolventCalculation.default_process_label(self)

    @classmethod
    def default_process_label(cls, node) -> str:
        """
        AiiDA Process label definition.

        Defines the process label to be associated with the created ProcessNode
        stored in the AiiDA database.

        Returns
        -------
        str
            The process label based on what inputs have been provided.
        """
        if node.inputs.do_opt_equillibrate.value:
            job_str = "_OptStep"
        elif node.inputs.do_init_optimise.value:
            job_str = "_OptStep"
        elif node.inputs.do_charge_fit.value:
            job_str = "_ESPChargeStep"
        elif node.inputs.do_md_equillibrate.value:
            job_str = "_MDStep"
        else:
            job_str = "_SPStep"
        return "Chemshell_Solvation" + job_str

    def chemsh_script_generator(self) -> str:
        """
        Generate the input script for a ChemShell Solvation Workflow.

        Returns
        -------
        script : str
            A string containing the ChemShell input script for the calculation.
        """
        qm_theory = None
        mm_theory = None

        script = "from chemsh import Fragment\n"
        if isinstance(self.inputs.structure, StructureData):
            if len(self.inputs.structure.sites) < 2:
                atom_names = [site.kind_name for site in self.inputs.structure.sites]
                coords = [
                    [UnitsConverter.angstrom_to_bohr(r) for r in site.position]
                    for site in self.inputs.structure.sites
                ]
                script += (
                    f"structure = Fragment(coords={str(coords):s}, names="
                    f"{str(atom_names):s})\n"
                )
            else:
                script += "structure = Fragment(coords="
                script += f"'{SolventCalculation.FILE_TMP_STRUCTURE:s}')\n"
        elif isinstance(self.inputs.structure, TrajectoryData):
            script += "structure = Fragment(coords="
            script += f"'{SolventCalculation.FILE_TMP_STRUCTURE:s}')\n"
        else:  # SinglefileData
            script += (
                f"structure = Fragment(coords='{self.inputs.structure.filename:s}')\n"
            )

        if "qm_parameters" in self.inputs:
            # Creates a quantum mechanics Theory object
            qm_theory = ChemShellQMTheory[
                self.inputs.qm_parameters.get("theory").upper()
            ]

            if qm_theory != ChemShellQMTheory.NONE:
                qm_theory_key = SolventCalculation.get_qm_theory_key(qm_theory)

                script += f"from chemsh import {qm_theory_key:s}\n"
                script += f"qmtheory = {qm_theory_key:s}(frag=structure"
                param_str = ""
                if "qm_parameters" in self.inputs:
                    for key in self.inputs.qm_parameters.keys():
                        if key == "theory":
                            continue
                        val = self.inputs.qm_parameters.get(key)
                        if isinstance(val, str):
                            param_str += ", " + key + "='" + val + "'"
                        else:
                            param_str += ", " + key + "=" + str(val)
                script += param_str + ")\n"

        script_opt = ""
        if self.inputs.do_init_optimise.value or self.inputs.do_opt_equillibrate.value:
            # Run a geometry optimisation using DL_FIND
            #rajany todo- determine if mm opt or qmmmopt option is necessary in the workflow

            script_opt += "from chemsh import Opt\n"
            theory_str = "qmtheory"
            #rajany todo
            #if qmmm or mm
            #theory_str = "mmtheory"
            #theory_str = "qmmmtheory"
            opt_str = f"job=Opt(theory={theory_str}"
            if "optimisation_parameters" in self.inputs:
                for key in self.inputs.optimisation_parameters.keys():
                    if isinstance(self.inputs.optimisation_parameters.get(key), str):
                        opt_str += ", " + key + "='"
                        opt_str += self.inputs.optimisation_parameters.get(key) + "'"
                    else:
                        opt_str += ", " + key + "="
                        opt_str += str(self.inputs.optimisation_parameters.get(key))
            script_opt += opt_str + ")\n"
            if self.inputs.dryrun:
                    script_opt += "job.run(dryrun=True)\njob.result.save()\n"
            else:
                    script_opt += "job.run(dryrun=False)\njob.result.save()\n"

            script += script_opt
            if (not self.inputs.optimisation_parameters.get("thermal", False) 
                   and self.inputs.optimisation_parameters.get("neb", "no") not in [
                   "free", "frozen", "perpendicular",
                    ]):
                script += f'structure.save("{SolventCalculation.FILE_DLFIND}")\n'
            return script

        script_ch = ""
        if self.inputs.do_charge_fit.value:
            # Run a Charge fitting task
            script_ch += "from chemsh import ChargeFitting\n"
            fit_str = f"job = ChargeFitting(theory = qmtheory"
            if "chargefitting_parameters" in self.inputs:
                for key in self.inputs.chargefitting_parameters.keys():
                    if isinstance(self.inputs.chargefitting_parameters.get(key), str):
                        fit_str += ", " + key + "='"
                        fit_str += self.inputs.chargefitting_parameters.get(key) + "'"
                    else:
                        fit_str += ", " + key + "="
                        fit_str += str(self.inputs.chargefitting_parameters.get(key))
            script_ch += fit_str + ")\n"
            if self.inputs.dryrun:
                    script_ch += "job.run(dryrun=True)\njob.result.save()\n"
            else:
                    script_ch += "job.run(dryrun=False)\njob.result.save()\n"

            script_ch += f"from numpy import column_stack, savetxt\n"
            script_ch += f"charges = column_stack([structure.names.astype(str), structure.charges])\n"
            script_ch += f"savetxt('{SolventCalculation.FILE_CHARGES}', charges, delimiter=' ', fmt='%s')\n"
            script    += script_ch
            return script

        # Perform an MD
        if self.inputs.do_md_equillibrate.value:
            if "mm_parameters" not in self.inputs:
                return("mm_parameters not provided for md step")
            else:
               # Creates a molecular mechanics Theory object
               mm_theory = ChemShellMMTheory[
                    self.inputs.mm_parameters.get("theory").upper()
               ]
               if mm_theory != ChemShellMMTheory.NONE:
                    mm_theory_key = SolventCalculation.get_mm_theory_key(mm_theory)

                    script += f"from chemsh import {mm_theory_key:s}\n"
                    param_str = ""
                    #if qmmm_chk:
                    #else:
                    script += f"mmtheory = {mm_theory_key:s}"
                    script += f"(ff='{self.inputs.force_field_file.filename:s}'"

                    for key in self.inputs.mm_parameters.keys():
                        if key == "theory":
                            continue
                        val = self.inputs.mm_parameters.get(key)
                        if isinstance(val, str):
                            param_str += ", " + key + "='" + val + "'"
                        else:
                            param_str += ", " + key + "=" + str(val)
                    script += f"{param_str:s})\n"

        # Create Solvent box structure object if requested
            if "solvent_box" in self.inputs:
                if isinstance(self.inputs.solvent_box, SinglefileData):
                    fname = self.inputs.solvent_box.filename
                else:
                    raise Exception("Solvent box type not recognized")
                script += f"solvent_structure = Fragment(coords='{fname:s}')\n"

            script += f"solute_structure = Fragment(coords='{SolventCalculation.FILE_TMP_STRUCTURE:s}')\n"

            theory_str = "mmtheory"
            if not self.inputs.solvent_box:
                raise Exception("Solvent box not provided")
                #qmmm_chk = "qm_parameters" in self.inputs and "mm_parameters" in self.inputs

                # If both QM and MM are specified, create a QM/MM interface object
                #if qmmm_chk:
                #    theory_str = "qmmm"
                #    script += "from chemsh import QMMM\n"
                #    script += "qmmm = QMMM(frag=structure, qm=qmtheory, mm=mmtheory, "
                #    qm_region_str = str(self.inputs.qmmm_parameters.get("qm_region", []))
                #    script += f"qm_region={qm_region_str:s})\n"
                #elif mm_theory:
                #    theory_str = "mmtheory"
                #else:
                #    theory_str = "qmtheory"

            script += "from chemsh import Solvation\n"
            script_md = ""
            script_md += f"job = Solvation(driver={theory_str:s}, solute=solute_structure, solvent=solvent_structure"
            for key in self.inputs.md_parameters.keys():

                    if key == "driver" or key == "solute" or key == "solvent":
                        continue
                    if isinstance(self.inputs.md_parameters.get(key), str):
                        script_md += ", " + key + "='"
                        script_md += self.inputs.md_parameters.get(key) + "'"
                    else:
                        script_md += ", " + key + "="
                        script_md += str(self.inputs.md_parameters.get(key))
            script_md += ")\n"

            if self.inputs.dryrunmd:
                    script_md += "job.run(dryrun=True)\njob.result.save()\n"
            else:
                    script_md += "job.run(dryrun=False)\njob.result.save()\n"

            script += script_md
            return script

 # Perform a single point energy calculation (default calculation type)
        if self.inputs.do_sp:
            script += "from chemsh import SP\n"
            if "calculation_parameters" not in self.inputs:
                # Assign default values if none are given
                self.inputs.calculation_parameters = Dict(dict={})


            theory_str = "qmtheory"
            script_qm = ""
            # Runs a QM single point energy calculation
            script_qm += f"job = SP(theory={theory_str:s}, "
            grad_str = str(self.inputs.calculation_parameters.get("gradients", False))
            script_qm += f"gradients={grad_str:s}, "
            hess_str = str(self.inputs.calculation_parameters.get("hessian", False))
            script_qm += f"hessian={hess_str:s})\n"

            if self.inputs.dryrun:
                    script_qm += "job.run(dryrun=True)\njob.result.save()\n"
            else:
                    script_qm += "job.run(dryrun=False)\njob.result.save()\n"

            script += script_qm

        return script

    def prepare_for_submission(self, folder: Folder) -> CalcInfo:
        """
        Prepare the ChemShell calculation for submission.

        Params
        ------
        folder : Folder
            An `aiida.common.folders.Folder` specifying the temporary working
            directory for the calculation.

        Returns
        -------
        calcInfo : CalcInfo
            An `aiida.common.CalcInfo` instance.
        """

        #rajany diag
        inputs_list = [self.inputs.do_opt_equillibrate.value, 
                       self.inputs.do_init_optimise.value,
                       self.inputs.do_charge_fit.value,
                       self.inputs.do_md_equillibrate.value]
        with folder.open('inputs.dat', 'w') as f:
            f.write(f"INPUTS = \n{inputs_list}\n")

        #end of diag

        # Create the ChemShell input script
        input_script = self.chemsh_script_generator()
        with folder.open(SolventCalculation.FILE_SCRIPT, "w") as f:
            f.write(input_script)

        # Define the AiiDA code parameters
        code_info = CodeInfo()
        code_info.code_uuid = self.inputs.code.uuid
        if "chemsh.x" in str(self.inputs.code.filepath_executable):
            code_info.cmdline_params = [
                SolventCalculation.FILE_SCRIPT,
            ]
        else:
            n_machines = self.inputs.metadata.options.resources.get("num_machines")
            n_mpi_pm = self.inputs.metadata.options.resources.get(
                "num_mpiprocs_per_machine"
            )
            tot_mpi = n_machines * n_mpi_pm
            code_info.cmdline_params = [
                "-np",
                self.inputs.metadata.options.resources.get("tot_num_mpiprocs", tot_mpi),
                SolventCalculation.FILE_SCRIPT,
            ]
        code_info.stdout_name = SolventCalculation.FILE_STDOUT

        # Setup the calculation information object
        calc_info = CalcInfo()
        calc_info.codes_info = [code_info]
        calc_info.retrieve_temporary_list = [
            SolventCalculation.FILE_RESULTS,
        ]
        calc_info.provenance_exclude_list = []
        calc_info.retrieve_list = [
            SolventCalculation.FILE_STDOUT,
        ]
        calc_info.local_copy_list = []

        if isinstance(self.inputs.structure, StructureData):
            with folder.open(SolventCalculation.FILE_TMP_STRUCTURE, "wb") as f:
                f.write(self.inputs.structure._prepare_xyz()[0])
        elif isinstance(self.inputs.structure, TrajectoryData):
            with folder.open(SolventCalculation.FILE_TMP_STRUCTURE, "wb") as f:
                index = self.inputs.structure_index.value
                structure = self.inputs.structure.get_step_structure(index=index)
                f.write(structure._prepare_xyz()[0])
        else:
            calc_info.local_copy_list.append(
                (
                    self.inputs.structure.uuid,
                    self.inputs.structure.filename,
                    self.inputs.structure.filename,
                ),
            )

        # Copy data for second input fragment if required
        if "structure2" in self.inputs:
            if isinstance(self.inputs.structure2, StructureData):
                with folder.open(SolventCalculation.FILE_TMP_STRUCTURE, "wb") as f:
                    f.write(self.inputs.structure2._prepare_xyz()[0])
            else:
                calc_info.local_copy_list.append(
                    (
                        self.inputs.structure2.uuid,
                        self.inputs.structure2.filename,
                        self.inputs.structure2.filename,
                    ),
                )

        # If running with an MM theory a force field file is required and copied
        if "force_field_file" in self.inputs:
            calc_info.local_copy_list.append(
                (
                    self.inputs.force_field_file.uuid,
                    self.inputs.force_field_file.filename,
                    self.inputs.force_field_file.filename,
                )
            )
        if "solvent_box" in self.inputs:
            calc_info.local_copy_list.append(
                (
                    self.inputs.solvent_box.uuid,
                    self.inputs.solvent_box.filename,
                    self.inputs.solvent_box.filename,
                )
            )

        # If performing a geometry optimisation retrieve the generated _dl_find.pun
        # file containing the optimised structure
        if "optimisation_parameters" in self.inputs:
            if not self.inputs.optimisation_parameters.get(
                "thermal", False
            ) and self.inputs.optimisation_parameters.get("neb", "no") not in [
                "free",
                "frozen",
                "perpendicular",
            ]:
                calc_info.retrieve_temporary_list.append(
                    SolventCalculation.FILE_DLFIND
                )
            if self.inputs.optimisation_parameters.get("save_path", False):
                calc_info.retrieve_temporary_list.append(
                    "_dl_find/" + SolventCalculation.FILE_TRJPTH
                )
                calc_info.retrieve_temporary_list.append(
                    "_dl_find/" + SolventCalculation.FILE_TRJFRC
                )
            if self.inputs.optimisation_parameters.get("neb", "no") in [
                "free",
                "frozen",
                "perpendicular",
            ]:
                calc_info.retrieve_temporary_list.append("nebinfo")
                calc_info.retrieve_temporary_list.append("nebpath.xyz")

        if "chargefitting_parameters" in self.inputs:
            calc_info.retrieve_list.append(f"{SolventCalculation.FILE_CHARGES}")
        if "md_parameters" in self.inputs:
            calc_info.retrieve_list.append(f"{SolventCalculation.FOLDER_SNAPSHOTS}")

        return calc_info
