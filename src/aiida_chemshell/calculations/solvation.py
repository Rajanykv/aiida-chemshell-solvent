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
)

from aiida_chemshell.units import UnitsConverter
from aiida_chemshell.utils import ChemShellMMTheory, ChemShellQMTheory


class SolventCalculation(ChemShellCalculation):
    """
    AiiDA calculation plugin wrapper for Solvation calculations to avoid clutter in base.

    Currently supports the following tasks:
      - Single point energy
      - Geometry optimisation

    """
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
            "do_init_optimise",
             valid_type = Bool,
             default=lambda: Bool(True),
             required = False,
             help = "Whether to do optimisation on the initial structure. (default True)"
        )
        spec.input(
            "do_charge_fit",
             valid_type = Bool,
             default=lambda: Bool(True),
             required = False,
             help = "Whether to do a QM RESP charge fitting(default True)"
        )
        spec.input(
            "do_md_equillibrate",
             valid_type = Bool,
             default=lambda: Bool(True),
             required = False,
             help = "Whether to do an MD equillibration step (default True)"
        )
        spec.input(
            "do_opt_equillibrate",
             valid_type = Bool,
             default=lambda: Bool(True),
             required = False,
             help = "Whether to do an optimisation instead of MD equillibration(default False)"
        )

        spec.input(
            "solvent_box",
            valid_type=(SinglefileData),
            validator=cls.validate_inputs_2,
            required=True,
            help=(
                "The solvent boxes input structure to be used for the ChemShell solvation"
                "Choose from the list in format '.pqr' or '.pdb'"
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
    def validate_inputs_2(self):
        """Validate the inputs provided to the WorkChain."""
        has_file = "structure_file" in self.inputs
        has_box = "solvent_box" in self.inputs
        #if not has_files and not has_box:
        if not has_box:
            return self.exit_codes.ERROR_NO_INPUTS
        return None

    @classmethod
    def get_valid_md_parameters(cls) -> dict[str:type]:
        """
        Return a tuple of valid parameter keys for the Chemshell MD calculation.

        Returns
        -------
        validKeys : dict[str: type]
            A tuple of valid parameter keys for the ChemShell Charge fitting calculation.
        """
        #rajany todo, copied all now. retain only required ones.
        return {
                'active'                :[],
                'boundary'              :'periodic',
                'constraints'           :[],
                'dcd'                   :'',
                'density_variance'      :   0.2,
                'driver'                :'dl_poly',
                'ensemble'              :'NPT',
                'ensemble_method'       :'langevin',
                'ensemble_barostat_coupling'  : 0.0, # in fs
                'ensemble_barostat_friction'  : 0.0, # in 1/fs
                'ensemble_thermostat_coupling': 0.0, # in fs
                'ensemble_thermostat_friction': 0.0, # in 1/fs
                'equilibrate'           : 0,
                'ff'                    :'charmm',
                'fix'                   :[],
                'freq_out_energy'       :-1,
                'freq_out_pressure'     :-1,
                'freq_rescale'          :-1,
                'frozen'                : full((1), -1, dtype=int64),
                'langevin'              : True,
                'langevin_temperature'  : None,
                'langevin_damping'      : 1,
                'langevin_H'            : False,
                'langevin_piston'       : False,
                'langevin_piston_target':   1.01325,        # bar
                'langevin_piston_period': 100.0,
                'langevin_piston_decay' :  50.0,
                'langevin_piston_temp'  : 293.15,
                'minimise'              : 100,
                'nsnapshots'            : 10,
                'nsteps'                : 100,
                'pbc_wrap'              :[],                # list of atom indices to perform PBC wrapping
                'plumed'                : plumed.PLUMED(),
                'pressure'              : 0.001,               # katm
                'result'                : resultmd.ResultMD(),
                'result_theory'         : None,
                'rigid'                 :'all',
                'save_trajectory'       : True,
                'seed'                  : 2020,
                'shake'                 : True,
                'shake_max_iter'        : 250,
                'shake_tolerance'       : 1e-05,               # Ang
                'spheric_r0'            :-1.0,
                'spheric_k'             : 10.0,
                'spheric_exponent'      : 2,
                'temperature'           : 293.15,
                'theory'                : None,
                'timestep'              : 1.0,               # femtosecond
                'velocities'            : zeros(shape=(1,3), dtype=float64),
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

        # Check for valid parameter values if value options are restricted
        #if "method" in value.keys():
        #    method = value.get("method").upper()
        #   if method not in ["ESP", "RESP"]:
        #      return f"The specified method key ('{method:s}') is not valid."

        return None

    def default_process_label(self) -> str:
        """
        AiiDA Process label definition.

        Defines the process label to be associated with the created ProcessNode
        stored in the AiiDA database.

        Returns
        -------
        str
            The process label based on what inputs have been provided.
        """
        return "Chemshell_Solvation_ Workflow"

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
        elif "structure_index" in self.inputs:
            print("Not yet supported.")
            raise Exception("SinglefileData trajectories not yet supported.")
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
                script += f"qmtheory = {qm_theory_key:s}(frag=structure"
                script += param_str + ")\n"

        script_opt = ""
        if "optimisation_parameters" in self.inputs:
            # Run a geometry optimisation using DL_FIND
            if self.inputs.do_init_optimise or self.inputs.do_opt_equillibrate:

                script_opt += "from chemsh import Opt\n"
                opt_str = f"job = Opt(theory=qmtheory"
                for key in self.inputs.optimisation_parameters.keys():
                    if isinstance(self.inputs.optimisation_parameters.get(key), str):
                        opt_str += ", " + key + "='"
                        opt_str += self.inputs.optimisation_parameters.get(key) + "'"
                    else:
                        opt_str += ", " + key + "="
                        opt_str += str(self.inputs.optimisation_parameters.get(key))
                script_opt += opt_str + ")\n"
                script_opt += "job.run()\njob.result.save()\n"
                script += script_opt

        # Perform a single point energy calculation (default calculation type)
        script += "from chemsh import SP\n"
        if "calculation_parameters" not in self.inputs:
            # Assign default values if none are given
            self.inputs.calculation_parameters = Dict(dict={})

        script_job = ""
        # Runs a QM single point energy calculation
        script_job += f"job = SP(theory={theory_str:s}, "
        grad_str = str(self.inputs.calculation_parameters.get("gradients", False))
        script_job += f"gradients={grad_str:s}, "
        hess_str = str(self.inputs.calculation_parameters.get("hessian", False))
        script_job += f"hessian={hess_str:s})\n"

        script_job += "job.run()\njob.result.save()\n"
        script += script_job

        script_ch = ""
        if "chargefitting_parameters" in self.inputs and self.inputs.do_charge_fit:
            # Run a Charge fitting task
            script = script.replace(script_opt, "")
            script = script.replace(script_job, "")
            script_ch += "from chemsh import ChargeFitting\n"
            fit_str = f"job = ChargeFitting(theory = qmtheory"
            for key in self.inputs.chargefitting_parameters.keys():
                if isinstance(self.inputs.chargefitting_parameters.get(key), str):
                    fit_str += ", " + key + "='"
                    fit_str += self.inputs.chargefitting_parameters.get(key) + "'"
                else:
                    fit_str += ", " + key + "="
                    fit_str += str(self.inputs.chargefitting_parameters.get(key))
            script_ch += fit_str + ")\n"
            script_ch += "job.run()\njob.result.save()\n"
            script_ch += f"from numpy import column_stack, savetxt\n"
            script_ch += f"charges = column_stack([structure.names.astype(str), structure.charges])\n"
            script_ch += f"savetxt('{SolventCalculation.FILE_CHARGES}', charges, delimiter=' ', fmt='%s')"
            script    += script_ch

        # Create Solvent box structure object if requested
        if "solvent_box" in self.inputs:
            if isinstance(self.inputs.solvent_box, SinglefileData):
                fname = self.inputs.solvent_box.filename
            else:
                raise Exception("Solvent box type not recognized")
            print(fname)
            script += f"box = Fragment(coords='{fname:s}')\n"

        ## Setup Theory objects

        # Perform an MD
        if "mm_parameters" in self.inputs and self.inputs.do_md_equillibrate:
            # Creates a molecular mechanics Theory object
            script = script.replace(script_opt, "")
            script = script.replace(script_job, "")
            script = script.replace(script_ch, "")
            mm_theory = ChemShellMMTheory[
                self.inputs.mm_parameters.get("theory").upper()
            ]
            if mm_theory != ChemShellMMTheory.NONE:
                mm_theory_key = SolventCalculation.get_mm_theory_key(mm_theory)

                script += f"from chemsh import {mm_theory_key:s}\n"
                param_str = ""
                for key in self.inputs.mm_parameters.keys():
                    if key == "theory":
                        continue
                    val = self.inputs.mm_parameters.get(key)
                    if isinstance(val, str):
                        param_str += ", " + key + "='" + val + "'"
                    else:
                        param_str += ", " + key + "=" + str(val)
                #if qmmm_chk:
                #    script += f"mmtheory = {mm_theory_key:s}"
                #    script += f"(ff='{self.inputs.force_field_file.filename:s}'"
                #    script += f"{param_str:s})\n"
                #else:
                #    script += f"mmtheory = {mm_theory_key:s}(frag=structure, "
                #    script += f"ff='{self.inputs.force_field_file.filename:s}'"
                #    script += f"{param_str:s})\n"

                script += f"mmtheory = {mm_theory_key:s}(frag=box, "
                script += f"ff='{self.inputs.force_field_file.filename:s}'"
                script += f"{param_str:s})\n"
            theory_str = "mmtheory"

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

        if "md_parameters" in self.inputs:
            script = script.replace(script_opt, "")
            script = script.replace(script_job, "")
            script = script.replace(script_ch, "")
            script += "from chemsh import MD\n"

            script_job = ""
            # Runs a QM single point energy calculation
            script_job += f"job = MD(theory={theory_str:s} "
            for key in self.inputs.md_parameters.keys():
                if isinstance(self.inputs.md_parameters.get(key), str):
                        script_job += ", " + key + "='"
                        script_job += self.inputs.optimisation_parameters.get(key) + "'"
                else:
                        script_job += ", " + key + "="
                        script_job += str(self.inputs.optimisation_parameters.get(key))
                script_job += script_job + ")\n"

    
            script_job += "job.run()\njob.result.save()\n"
            script += script_job

        if "optimisation_parameters" in self.inputs:
            if not self.inputs.optimisation_parameters.get(
                "thermal", False
            ) and self.inputs.optimisation_parameters.get("neb", "no") not in [
                "free",
                "frozen",
                "perpendicular",
            ]:
                script += f'structure.save("{SolventCalculation.FILE_DLFIND}")\n'

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

        return calc_info
