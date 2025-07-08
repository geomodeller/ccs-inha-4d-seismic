import numpy as np
import os
import shutil
from datetime import datetime
import pyvista as pv
class Sim_CCS():
    """    
    # author: Honggeun Jo
    # Date: 2024-01-20
    __version__: Returns the version of the class.
    
    __init__: Initializes the object with the current directory, template directory, and total injection amount.
    reset_injection_amount: Resets the total injection amount.
    reset_parameters_to_play: Resets the simulation parameters.
    reset_bhp_constraint: Resets the bottom hole pressure constraint.
    reset_grid_size: Resets the grid size.
    reset_Ptime: Resets the simulation time.
    reset_CMG_exe: Resets the path to the CMG executable.
    run_simulation: Runs a simulation with given input parameters and an optional bottom hole pressure constraint.
    read_grid_results: Reads grid results from the specified run directory.
    read_well_results: Reads well results from a specified run directory.
    remove_out_files: Removes the output files associated with the specified run directory.
    """
    version = '0.1.0'
    CMG =  '"C:\\Program Files\\CMG\\GEM\\2022.10\\Win_x64\\EXE\\gm202210.exe"'

    scf_to_kg = 0.05189
    scf_to_ton = scf_to_kg / 1E3
    total_injection_per_year = 300_000 # tons/year
    total_injection_per_day = total_injection_per_year/365
    total_injection_per_day_in_scf = total_injection_per_day/scf_to_ton

    NX, NY, NZ = 32, 32, 16
    Ptime = [2025,2026,2027]    
    list_of_property = ['TotalGasinPlace,SCcuft', 'GasSaturation', 'WaterSaturation', 'Pressure(psia)']

    currert_directory = '.'
    template_dir   = os.path.join(currert_directory,'template')
    CMG_GEM_data_file = os.path.join(template_dir, 'cmg_ccs_run_file_v2.dat') 
    CMG_GEM_inc_files = []
    inc_files_list =  ['fluid_data.inc','permeability.inc','porosity.inc',
                        'rock_fluid_data.inc','rocktype.inc',
                        'xcorn.inc','ycorn.inc', 'zcorn.inc']
    for inc_file in inc_files_list:
        CMG_GEM_inc_files.append(os.path.join(template_dir, inc_file))

    params = ['well_1_perf_top','well_1_perf_bottom',
              'well_2_perf_top','well_2_perf_bottom',
              'well_3_perf_top','well_3_perf_bottom',
              'well_4_perf_top','well_4_perf_bottom',
              'sch_ratio_year_1',
              'sch_ratio_year_2',
              'sch_ratio_year_3']

    bhp_constraint = None

    def __version__(self):
        print(f'The current version is {self.__version__}')

    def __init__(self, 
                 currert_directory = None, 
                 template_dir = None,
                 total_injection_ammount = None):
        """
        Initializes the object with the current directory, template directory, and total injection amount.

        Parameters:
            currert_directory (str): The current directory path.
            template_dir (str): The template directory path.
            total_injection_ammount (int or float): The total injection amount.

        Returns:
            None
        """
        
        if currert_directory is not None:
            self.currert_directory = currert_directory
        if template_dir is not None:
            self.template = os.path.join(self.currert_directory,template_dir)

        if type(total_injection_ammount) in [int, float]:
            self.total_injection_ammount = total_injection_ammount
        elif total_injection_ammount is None:
            print(f'Default injection amount is {self.total_injection_per_day_in_scf:,.0f} scf/day')
            self.total_injection_ammount = self.total_injection_per_day_in_scf
        else:
             assert False, "only int and float are allowed for total_injection_ammount"

        for file in self.CMG_GEM_inc_files:
            assert os.path.isfile(file), f'{file} does not exist'
        assert os.path.isfile(self.CMG_GEM_data_file), f'{self.CMG_GEM_data_file} does not exist'

    def reset_injection_amout(self, total_injection_ammount = None):
        if type(total_injection_ammount) in [int, float]:
            self.total_injection_ammount = total_injection_ammount
        else:
            assert False, "only int and float are allowed"
    def reset_parameters_to_play(self, keys):
        self.params = keys 
    def reset_bhp_constraint(self, bhp_constraint):
        self.bhp_constraint = bhp_constraint 
    def reset_grid_size(self, nx, ny, nz):
        self.NX, self.NY, self.NZ = int(nx), int(ny), int(nz)
    def reset_Ptime(self, Ptime):
        self.Ptime = Ptime
    def reset_CMG_exe(self, CMG):
        self.CMG = CMG
    def run_simulation(self, run_dir = 'test_run', input_parameters = {}, bhp_constraint = None):
        """
        Run a simulation with the given input parameters and optional BHP constraint.

        Args:
            run_dir (str): The directory for the simulation run.
            input_parameters (dict): Dictionary of input parameters for the simulation.
            bhp_constraint (float or int, optional): BHP constraint for the simulation.

        Returns:
            None
        """
        out_file =  self.CMG_GEM_data_file.split('\\')[-1].replace('.dat','.out')
        data_file_full_dir =  os.path.join(run_dir, self.CMG_GEM_data_file.split('\\')[-1])
        data_file = self.CMG_GEM_data_file.split('\\')[-1]

        if bhp_constraint is not None:
            assert type(bhp_constraint) in [float, int], f'bhp_constraint should be either float or int'
            self.bhp_constraint = bhp_constraint

        if os.path.isdir(run_dir) is False: 
            os.mkdir(run_dir)
        if os.path.isdir(run_dir) and os.path.isfile(os.path.join(run_dir, out_file)):
            print('simulation result is already there ... please check your run_directory')
            return
        for keyword in input_parameters.keys():
            assert keyword in self.params, f'{keyword} is not defined'
        for keyword in self.params:
            assert keyword in input_parameters.keys(), f'{keyword} is missing in the input parameter'

        # copy inc files:
        for file in self.CMG_GEM_inc_files:
            shutil.copyfile(file, os.path.join(run_dir, file.split('\\')[-1]))

        well_params = []
        well_params_values = []
        sch_params = []
        sch_params_values = []

        for key, val in input_parameters.items():
            if key.startswith('well'):
                well_params.append('$'+key)
                well_params_values.append(val)
            if key.startswith('sch'):
                sch_params.append('$'+key)
                sch_params_values.append(val)
        
        # copy data file with replacing well contraints:
        shutil.copyfile(self.CMG_GEM_data_file, data_file_full_dir)
        self._search_and_replace(data_file_full_dir, well_params, well_params_values)
        
        # create schedules:
        for i, ratio in enumerate(sch_params_values):
            self._make_new_schedule(ratio = ratio, step = i+1, run_dir = run_dir)

        # run sim:
        os.chdir(run_dir)
        os.system(f"{self.CMG} -f {data_file}")
        os.chdir('..')

    def read_grid_results(self, run_dir = 'test_run'):
        """
        Reads grid results from the specified run directory.

        Args:
            run_dir (str): The directory where the run results are located. Defaults to 'test_run'.

        Returns:
            The grid results read from the output file.
        """
        output_file = os.path.join(run_dir, self.CMG_GEM_data_file.split('\\')[-1].replace('.dat', '.out'))
        assert os.path.isfile(output_file), 'no output file located - check if the simulation is done'
        grid_out = self._read_out_grid(output_file)
        return grid_out
    def read_well_results(self, run_dir = 'test_run'):
        """
        Reads well results from a specified run directory.

        Args:
            run_dir (str): The directory from which to read the well results. Defaults to 'test_run'.

        Returns:
            dict: A dictionary containing the well results.
        """
        output_file = os.path.join(run_dir, self.CMG_GEM_data_file.split('\\')[-1].replace('.dat', '.out'))
        assert os.path.isfile(output_file), 'no output file located - check if the simulation is done'
        well_out = self._read_out_well(output_file)
        print('outputs are: time, date, WBHP, WGBP, WWIR, WCWI ')
        return well_out
    
    def remove_out_files(self, run_dir = 'test_run'):
        """
        Removes the output files associated with the specified run directory.

        :param run_dir: The directory from which to remove the output files (default is 'test_run')
        :type run_dir: str
        """
        output_files = [os.path.join(run_dir, self.CMG_GEM_data_file.split('\\')[-1].replace('.dat', i)) for i in ['.out', '.rst', '.sr3']]
        for file in output_files:
            if os.path.isfile(file):
                os.remove(file)

    def _read_out_well(self, output_file, num_well=1,) :
        """
        Reads data from a specified output file and extracts well information.
        
        Args:
            output_file (str): The path to the output file to be read.
            num_well (int, optional): The number of wells to read. Defaults to 1.
        
        Returns:
            dict: A dictionary containing the extracted well information, including
            WBHP, WGBP, WWIR, WCWI, time, and date.
        """
        f = open(output_file, 'r')
        f_ = f.readlines()
        f.close()
                
        time = []
        date = []
        WBHP = []
        WGBP = [] 
        WWIR = [] 
        WCWI = []

        flag = 0        
        for i in range(len(f_)):   
            line = f_[i].split()        
            if flag == 0:
                if [''.join(line[3 : -2])][0] == 'GEMFIELDSUMMARY':
                    time.append(float(f_[i].split('TIME:')[-1].split('days')[0]))
                    date.append(datetime.strptime(f_[i].split('DATE:')[-1].split('\n')[0].replace(' ',''), '%Y:%b:%d'))
                    flag = 1                    
            if flag == 1:
                WBHP.append(np.array(f_[i+17].split('+')[1:5],dtype = float))   #WBHP, psia
                WGBP.append(np.array(f_[i+18].split('+')[1:5],dtype = float))   #WGBP, psia                
                WWIR.append(np.array(f_[i+29].split('+')[1:5],dtype = float))   #WWIR, STB/day
                WCWI.append(np.array(f_[i+31].split('+')[1:5],dtype = float))   #WCWI, MSTB       
                
                flag = 0;      
        WBHP, WGBP, WWIR, WCWI = np.array(WBHP), np.array(WGBP), np.array(WWIR), np.array(WCWI)
        time, date = np.array(time), np.array(date)
        well = {'WBHP':WBHP, 'WGBP':WGBP, 'WWIR':WWIR, 'WCWI':WCWI, 'time': time, 'date': date}
        return well


    def _read_out_grid(self, output_file):
        """
        Read the output grid from the specified file and parse the data into a dictionary 
        containing the grid properties for each property in the list_of_property. 
        It returns the result as a dictionary.
        """
        f = open(output_file, 'r')
        lines = f.readlines()
        f.close()
        nz, ny, nx = self.NZ,self.NY,self.NX
        result = {}
        for prop in self.list_of_property: 
            grid_property = - np.ones((len(self.Ptime)+1,nz, ny, nx))
            t_flag = 0

            for i in range(len(lines)):
                # Read Line:
                line = lines[i].split()
                ## Read grid property:
                if [''.join(line)][0] == prop:
                    start = i + 2
                    for j in range(start,len(lines)):
                        if '*' * 70 in lines[j].split(' '):
                            end = j
                            break
                        elif '*' * 131 + '\n' in lines[j].split(' '):
                            end = j
                            break
                    val_lines = lines[start:end]

                    if not any(['Plane K = ' in var for var in val_lines]):
                        value = [''.join(val_lines[np.where(['All values are' in val for val in val_lines])[0][0]].split(' '))][0]
                        value = value.split('Allvaluesare')[1].split('\n')[0]
                        grid_property[t_flag] = value 
                        t_flag += 1
                        
                    else:
                        value_line_lst = list(np.where(['Plane K =' in  val for val in val_lines])[0])
                        value_line_lst.append(len(val_lines))
                        for v in range(len(value_line_lst)-1):
                            s, e = value_line_lst[v], value_line_lst[v+1]
                            values =  val_lines[s:e]
                            k_flag = int(val_lines[s:e][0].split(' Plane K = ')[-1].split(' ')[0]) - 1 
                            if any(['All values are' in val for val in values]):
                                value = [''.join(values[np.where(['All values are' in val for val in values])[0][0]].split(' '))][0]
                                value = value.split('Allvaluesare')[1].split('\n')[0]
                                grid_property[t_flag, k_flag] = float(value)

                            else:
                                i_index = np.where(['    I = ' in val for val in values])[0]
                                for tick in i_index: 
                                    temp = [i for i in values[tick].split('I =')[-1].split('     ') if i != '']
                                    i_flags = np.array(temp, dtype = int) -1
                                    f = lambda x: x if not x in ['', '\n'] else None
                                    j_flags = [list(filter(f,val.split(' J= ')[-1].split(' ')))[0] for val in values[tick+1:tick+ny+1]]
                                    j_flags = np.array(j_flags, dtype = int) - 1
                                    value = [list(filter(f,val.split(' J= ')[-1].split(' ')))[1:] for val in values[tick+1:tick+ny+1]]
                                    value = [ [val[:-1] if val[-1] in ['p','i'] else val for val in vals] for vals in value]
                                    value = np.array(value, dtype = float)

                                    for j_, val_ in zip(j_flags,value):
                                        grid_property[t_flag, k_flag, j_, list(i_flags)[0]:list(i_flags)[-1]+1] = val_
                                
                        t_flag += 1

            result[prop] = grid_property  
        print(f'read output:{self.list_of_property}')
        return result

    def _make_new_schedule(self, 
                          run_dir = '', 
                          step = 1,   
                          ratio = 4*[0.25],):
        """
        Creates a new schedule file for the given run directory, step, and injection ratio.

        Parameters:
            run_dir (str): The directory where the schedule file will be created.
            step (int): The step number for the schedule file.
            ratio (list): A list of injection ratios for each well.

        Returns:
            None
        """
        inject_rates = np.array(ratio)*self.total_injection_ammount
        with open(f'{run_dir}/operation_{step}.inc', 'w') as file:
            for i, rate in enumerate(inject_rates, 1):
                file.write(f"INJECTOR  'Well {i}' \n")
                file.write(f"OPERATE  MAX  STG  {rate}  CONT REPEAT \n")
                if self.bhp_constraint is not None:
                    file.write(f"OPERATE  MAX  BHP  {self.bhp_constraint}  CONT REPEAT \n")

                
    def _search_and_replace(self, 
                            new_path,
                            search_word, 
                            replace_word):
        """
        Search and replace specified words in the file contents and save the modified contents to a new file.

        Args:
            new_path (str): The path of the new file to be created.
            search_word (list): The list of words to search for in the file contents.
            replace_word (list): The list of words to replace the corresponding search words with.

        Returns:
            None
        """
        with open(self.CMG_GEM_data_file, 'r') as file:
            file_contents = file.read()
            for search, replace in zip(search_word, replace_word):
                if type(replace) is not str:
                    replace = str(replace)
                file_contents = file_contents.replace(search, replace)
        with open(new_path, 'w') as file:
            file.write(file_contents)

    def _copy_and_paste_inc_files(self, run_dir = None):
        """
        Copies and pastes the CMG gem include files to a specified directory.

        :param run_dir: the directory where CMG will be run
        :type run_dir: str
        """
        assert run_dir is not None, 'run_dir is needed to be a specific path where to run CMG'
        targets = [os.path.join(run_dir, i.split('\\')[-1]) for i in CMG_GEM_inc_files]
        for source, target in zip(self.CMG_GEM_inc_files, targets):
           shutil.copyfile(source, target)

    def visual_3d(self, property, aspect_x_to_z = 1/3, show_edges=False, cmap = 'viridis'):
        """
        Visualize 3D data with optional parameters for aspect ratio, edge visibility, and color map.
        
        :param property: The 3D property to visualize
        :param aspect_x_to_z: The aspect ratio of x to z (default is 1/3)
        :param show_edges: Whether to show the edges (default is False)
        :param cmap: The color map to use (default is 'viridis')
        """                
        grid = pv.ImageData()
        grid.dimensions = np.array([self.NX, self.NY, self.NZ]) + 1
        grid.origin = (1, 1, 1)  # The bottom left corner of the data set
        grid.spacing = (1, 1, aspect_x_to_z)  # These are the cell sizes along each axis
        grid.cell_data["values"] =property[::-1].T.flatten(order="F")  # Flatten the array

        grid.plot(show_edges=False, cmap= cmap)
        
