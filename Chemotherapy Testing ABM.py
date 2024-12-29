# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import copy
from matplotlib.colors import ListedColormap
import matplotlib.colors as mcolors
from matplotlib.animation import FuncAnimation
import time
import pickle
class Cell():
    def __init__(self, row, column, size, is_cancer, args = None):
        """
        Initializes a cell with its attributes.

        Parameters
        ----------
        row : int
            The row position of the cell in the grid.
        column : int
            The column position of the cell in the grid.
        size : int
            The size of the grid.
        is_cancer : bool
            Whether the cell is a cancer cell or not.
        args : list of float, optional
            A list containing mutation, invasion, growth, resistance, recovery, and age parameters for cancer cells. The default is None.

        Returns
        -------
        None.
        """
        self.is_cancer = is_cancer
        self.is_dying = False  # Initially, the cell is not dying
        self.is_invading = False  # Initially, the cell is not invading
        self.y_position = column  # Set the column position
        self.x_position = row  # Set the row position
        self.cell_cycle_stage = "G0"  # Initial stage of the cell cycle is G0 (resting phase)
        self.stage_duration = None  # Duration of each stage (to be set later)
        self.stage_counter = 0  # Counter for tracking progression through stages
        self.oxygen_level = .3  # Initial oxygen level
        self.nutrient_level = .3  # Initial nutrient level
        self.chemo_level = 0.0  # Initial chemotherapy level
        self.neighbors = self.initialize_neighbors(row, column, size)  # Initialize neighboring cells
        self.num_mutation = 0  # Counter for mutations
        self.damage = 0.0  # Damage level, starts at 0
        
        # If the cell is a cancer cell, assign cancer-related attributes
        if is_cancer:
            number = int(np.random.normal(loc=5000, scale=200))  # Cancer cell lifespan
            self.lifespan = max(0, number)  # Ensure non-negative lifespan
            if args == None:
                # Default cancer parameters
                self.mutation_probability = 0.02
                self.invasion_prob = 0.2
                self.growth_acceleration_factor = 2
                self.resistance_apoptosis = 0.2
                self.recovery_rate = 0.3
                self.age = np.random.randint(0, int(self.lifespan * .5))
        else:
            # Normal cell lifespan
            number = int(np.random.normal(loc=1000, scale=20))
            self.lifespan = max(0, number)  # Ensure non-negative lifespan
            if args == None:
                # Default normal cell parameters
                self.mutation_probability = 0.000001
                self.invasion_prob = 0.0 
                self.growth_acceleration_factor = 1
                self.resistance_apoptosis = 0.0 
                self.recovery_rate = 0.4
                self.age = np.random.randint(0, self.lifespan)
                
        # If custom parameters (args) are provided, set them
        if args != None:
            self.num_mutation = args[0]
            self.mutation_probability = args[1]
            self.invasion_prob = args[2]
            self.growth_acceleration_factor = args[3]
            self.resistance_apoptosis = args[4]
            self.recovery_rate = args[5]
            self.age = 0

    def update_cell(self, environment_nutrient, environment_oxygen, environment_chemo, is_empty_neighbor):
        """
        Updates the cell's internal state based on the environmental conditions.

        Parameters
        ----------
        environment_nutrient : float
            The nutrient level in the environment.
        environment_oxygen : float
            The oxygen level in the environment.
        environment_chemo : float
            The chemotherapy level in the environment.
        is_empty_neighbor : bool
            Whether the cell has an empty neighboring space for division.

        Returns
        -------
        bool
            Whether the cell is ready to divide.
        """
        if self.is_dying == False:
            self.update_internal_concentrations(environment_nutrient, environment_oxygen, environment_chemo)
            
            # If oxygen or nutrient levels are too low, cell may die unless resistant to apoptosis
            if self.oxygen_level < 0.05 or self.nutrient_level < 0.05:
                if np.random.random() > self.resistance_apoptosis:
                    self.is_dying = True
                    self.stage_counter = 0
                    
            # Set damage multiplier based on cell cycle stage (higher damage during S and M phases)
            if self.cell_cycle_stage == "S" or self.cell_cycle_stage == "M":
                damage_multiplier = 1.5
            elif self.cell_cycle_stage == "G1" or self.cell_cycle_stage == "G2":
                damage_multiplier = 1 
            else: 
                damage_multiplier = .5
            self.damage += self.chemo_level * damage_multiplier - self.recovery_rate 
            
            # Cap damage between 0 and 1
            if self.damage < 0:
                self.damage = 0.0
            if self.damage > 1:
                self.damage = 1
            
            # If the cell reaches its maximum age, it dies
            if self.age == self.lifespan:
                self.is_dying = True
                self.stage_counter = 0
            
            # Random chance of mutation occurring
            if np.random.random() <= self.mutation_probability:
                self.mutate_cell()
            
            # Check if the damage exceeds the resistance to apoptosis
            if np.random.random() < self.damage - self.resistance_apoptosis:
                self.is_dying = True
                self.stage_counter = 0
                
            self.age += 1  # Increment the cell's age
            
            # If there's space for the cell to divide and it's in G0 phase, start division
            if is_empty_neighbor and self.cell_cycle_stage == "G0":
                self.stage_duration = np.zeros(4)
                self.stage_duration[0] = int(np.random.normal(loc= 11 / self.growth_acceleration_factor, scale=3))
                self.stage_duration[1] = int(np.random.normal(loc= 8 / self.growth_acceleration_factor, scale=3))
                self.stage_duration[2] = int(np.random.normal(loc= 4 / self.growth_acceleration_factor, scale=2))
                self.stage_duration[3] = 1 
                self.cell_cycle_stage = "G1"
                self.stage_counter = 0
                is_ready_to_split = False
            
            # Progress through the cell cycle stages
            elif self.cell_cycle_stage == "G1":
                if self.stage_counter == self.stage_duration[0]:
                    self.stage_counter = 0
                    self.cell_cycle_stage = "S"   
                else:
                    if self.oxygen_level > 0.1 and self.nutrient_level > 0.1:
                        self.stage_counter += 1
                is_ready_to_split = False
                
            elif self.cell_cycle_stage == "S":
                if self.stage_counter == self.stage_duration[1]:
                    self.stage_counter = 0
                    self.cell_cycle_stage = "G2"
                else:
                    self.stage_counter += 1 
                is_ready_to_split = False
                
            elif self.cell_cycle_stage == "G2":
                if self.stage_counter == self.stage_duration[2]:
                    self.stage_counter = 0
                    self.cell_cycle_stage = "M"
                else:
                    if self.oxygen_level > 0.1 and self.nutrient_level > 0.1:
                        self.stage_counter += 1
                is_ready_to_split = False
            elif self.cell_cycle_stage == "M":
                if self.stage_counter == self.stage_duration[3]:
                    self.stage_counter = 0
                    self.cell_cycle_stage = "G0"
                    if is_empty_neighbor or self.is_invading:
                        is_ready_to_split = True  # Cell is ready to divide
                        self.is_invading = False  # Reset invasion status
                    else: 
                        is_ready_to_split = False
                else:
                    self.stage_counter += 1
                    is_ready_to_split = False
                        
            # If cell is ready to invade, start the invasion process
            elif np.random.random() < self.invasion_prob:
                self.stage_duration = np.zeros(4)
                self.stage_duration[0] = int(np.random.normal(loc= 11 / self.growth_acceleration_factor, scale=3))
                self.stage_duration[1] = int(np.random.normal(loc= 8 / self.growth_acceleration_factor, scale=3))
                self.stage_duration[2] = int(np.random.normal(loc= 4 / self.growth_acceleration_factor, scale=2))
                self.stage_duration[3] = 1 
                self.cell_cycle_stage = "G1"
                self.is_invading = True  # Mark the cell as invading
                self.stage_counter = 0
                is_ready_to_split = False
            
            else:
                is_ready_to_split = False
            
            return is_ready_to_split  # Return whether the cell is ready to divide
        
        else:
            self.stage_counter += 1  # If the cell is dying, just increase the stage counter
            is_ready_to_split = False
            return is_ready_to_split
    
    def get_cell_attributes(self):
        """
        Returns the attributes of the cell.

        Returns
        -------
        tuple
            A tuple containing mutation, invasion, growth, resistance, recovery, and cancer status attributes.
        """
        return self.num_mutation, self.mutation_probability, self.invasion_prob, self.growth_acceleration_factor, self.resistance_apoptosis, self.recovery_rate, self.is_cancer
            
    def initialize_neighbors(self, row, column, size):
        """
        Initializes the neighbors of the cell.

        Parameters
        ----------
        row : int
            The row position of the cell in the grid.
        column : int
            The column position of the cell in the grid.
        size : int
            The size of the grid.

        Returns
        -------
        list of tuples
            A list of neighbor coordinates in the form (row, column).
        """
        neighbors = []
        # List of relative neighbor offsets
        neighbor_offsets = [
            (-1, -1), (-1, 0), (-1, 1),
            ( 0, -1),          ( 0, 1),
            ( 1, -1), ( 1, 0), ( 1, 1)]
        
        for dr, dc in neighbor_offsets:
            neighbor_row = row + dr
            neighbor_col = column + dc
            if 0 <= neighbor_row < size and 0 <= neighbor_col < size:
                neighbors.append((neighbor_row, neighbor_col))
        return neighbors
    
        
    def update_internal_concentrations(self, environment_nutrient, environment_oxygen, environment_chemo):
        """
        Updates the internal concentrations of nutrients, oxygen, and chemotherapy agents based on the external environment.
    
        Parameters
        ----------
        environment_nutrient : float
            The current concentration of nutrients in the environment.
        environment_oxygen : float
            The current concentration of oxygen in the environment.
        environment_chemo : float
            The current concentration of chemotherapy agents in the environment.
    
        Returns
        -------
        None
        """
        
        # Update oxygen level based on the difference between the environment and internal oxygen
        if environment_oxygen > self.oxygen_level:
            self.oxygen_level += (environment_oxygen - self.oxygen_level) * environment_oxygen
        else:
            self.oxygen_level += (environment_oxygen - self.oxygen_level) * self.oxygen_level
        
        # Cap the oxygen level at 1 (maximum saturation)
        if self.oxygen_level > 1:
            self.oxygen_level = 1
        
        # Update nutrient level based on the difference between the environment and internal nutrient
        if environment_nutrient > self.nutrient_level:
            self.nutrient_level += (environment_nutrient - self.nutrient_level) * environment_nutrient
        else:
            self.nutrient_level += (environment_nutrient - self.nutrient_level) * self.nutrient_level
        
        # Cap the nutrient level at 1 (maximum saturation)
        if self.nutrient_level > 1:
            self.nutrient_level = 1
        
        # Update chemotherapy level based on the difference between the environment and internal chemotherapy level
        if environment_chemo > self.chemo_level:
            self.chemo_level += (environment_chemo - self.chemo_level) * environment_chemo
        else:
            self.chemo_level += (environment_chemo - self.chemo_level) * self.chemo_level
        
        # Cap the chemotherapy level at 1 (maximum saturation) and set it to 0 if below threshold
        if self.chemo_level > 1:
            self.chemo_level = 1
        if self.chemo_level < 0.05:
            self.chemo_level = 0

    def mutate_cell(self):
        """
        Simulates a mutation process, altering the cell's attributes with certain probabilities.
    
        Parameters
        ----------
        None
    
        Returns
        -------
        None
        """
        
        # Generate a random number to decide which mutation occurs
        random_number = np.random.random()
        
        # Mutation to increase mutation probability
        if random_number < .2:
            if self.mutation_probability < .99:
                self.mutation_probability += .01
        
        # Mutation to increase invasion probability
        if random_number >= .2 and random_number < .4:
            if self.invasion_prob < .98:
                self.invasion_prob += .02
        
        # Mutation to increase growth acceleration factor
        if random_number >= .4 and random_number < .6:
            self.growth_acceleration_factor += .05
        
        # Mutation to increase resistance to apoptosis
        if random_number >= .6 and random_number < .8:
            if self.resistance_apoptosis < .98:
                self.resistance_apoptosis += .02
        
        # Mutation to increase recovery rate
        if random_number >= .8 and random_number < 1:
            if self.recovery_rate < .8:
                self.recovery_rate += .01
            
        # Increase the cell's damage from the mutation
        self.damage += .05
        
        # Cap damage at 1 (maximum damage)
        if self.damage > 1:
            self.damage = 1
        
        # Increment the mutation counter
        self.num_mutation += 1


class Extracellular_Matrix():
    def __init__(self, size):
        """
        Initializes the extracellular matrix (ECM) with nutrient, oxygen, and chemotherapy concentrations.

        Parameters
        ----------
        size : int
            The size of the matrix (both dimensions are equal to `size`).

        Returns
        -------
        None
        """
        
        # Size of the ECM and its matrix (3 channels for nutrients, oxygen, and chemo)
        self.size = size                              
        self.ecm = np.zeros((size, size, 3), dtype=float)  # Initialize the ECM with zeros
        self.vasculature = np.zeros((size, size))         # Initialize the vasculature matrix
        self.midline = size // 2                         # Midline index to divide the matrix
        
        # Set vasculature along the midline of the ECM
        self.vasculature[:, self.midline - 1] = 1
        
        # Set nutrient gradient across the midline (0 to 1 from the center)
        nutrient_step_size = 0.02
        self.ecm[:,:self.midline, 0] = np.arange(1 - (self.midline) * nutrient_step_size, 1, nutrient_step_size)
        self.ecm[:, self.midline + 1:, 0] = np.arange(1 - nutrient_step_size, 1 - (self.midline + 1) * nutrient_step_size, -nutrient_step_size)
        self.ecm[:, self.midline, 0] = 1  # Set nutrient level at the midline
        
        # Set oxygen gradient across the midline (0 to 1 from the center)
        oxygen_step_size = 0.01
        self.ecm[:,:self.midline, 1] = np.arange(1 - (self.midline) * oxygen_step_size, 1, oxygen_step_size)
        self.ecm[:, self.midline + 1:, 1] = np.arange(1 - oxygen_step_size, 1 - (self.midline + 1) * oxygen_step_size, -oxygen_step_size)
        self.ecm[:, self.midline, 1] = 1  # Set oxygen level at the midline
    
    def set_chemo_concentration(self, max_concentration):
        """
        Sets the chemotherapy agent concentration gradient in the extracellular matrix.

        Parameters
        ----------
        max_concentration : float
            The maximum concentration of chemotherapy at the center (midline).

        Returns
        -------
        None
        """
        
        # Set chemo concentration gradient across the midline (from center to edges)
        chemo_step_size = 0.04
        self.ecm[:, :self.midline, 2] = np.arange(max_concentration - (self.midline) * chemo_step_size, max_concentration, chemo_step_size)
        self.ecm[:, self.midline + 1:, 2] = np.arange(max_concentration - chemo_step_size, max_concentration - (self.midline + 1) * chemo_step_size, -chemo_step_size)
        self.ecm[:, self.midline, 2] = max_concentration  # Set chemo level at the midline
        
        # Apply threshold to set concentrations below 0.005 to zero
        self.ecm[:,:,2] = np.where(self.ecm[:,:,2] > 0.005, self.ecm[:,:,2], 0)
    
    def reduce_chemo_concentration(self, halflife):
        """
        Reduces the chemotherapy agent concentration in the extracellular matrix based on its half-life.

        Parameters
        ----------
        halflife : float
            The half-life of the chemotherapy agent, representing the time required for its concentration to decrease by half.

        Returns
        -------
        None
        """
        
        # Reduce the chemotherapy concentration by the half-life factor
        self.ecm[:,:,2] = self.ecm[:,:,2] * (.5 ** (1/halflife))
        
        # Apply threshold to set concentrations below 0.005 to zero
        self.ecm[:,:,2] = np.where(self.ecm[:,:,2] > 0.005, self.ecm[:,:,2], 0)

class Tissue():
    '''
    A class representing a tissue environment, containing a grid of cells, an extracellular matrix, and 
    methods to simulate cell behavior, track various metrics, and visualize the tissue structure.

    Parameters
    ----------
    size : int
        The size of the tissue grid (size x size).

    Returns
    -------
    None
    '''
    
    def __init__(self, size):
        '''
        Initializes the tissue grid with cells and an extracellular matrix.

        Parameters
        ----------
        size : int
            The size of the tissue grid (size x size).

        Returns
        -------
        None
        '''
        self.size = size
        
        # Initialize an empty grid to store cells. None represents an empty space.
        self.cell_matrix = np.full((size, size), None, dtype=object)
        
        # Populate the tissue with normal cells, and randomly place cancer cells (1% chance).
        for row in range(size):
            for column in range(size):
                rand_chance = np.random.uniform()  # Random chance to place a cancer cell.
                if rand_chance < .99:
                    # Place normal cells with 99% probability.
                    self.cell_matrix[row, column] = Cell(row, column, size, is_cancer=False)
                else:
                    # Leave the spot empty for cancer cells (None).
                    self.cell_matrix[row, column] = None
        
        # Initialize the extracellular matrix for environmental conditions.
        self.extracellular_matrix = Extracellular_Matrix(size)
        
        # Generate all possible indices (row, column) for the tissue grid.
        a = np.indices((size, size))
        self.all_indices = list(zip(a[0].ravel(), a[1].ravel()))

    def simulate_step(self, chemo_concentration=None):
        '''
        Simulates one step of cell behavior based on environmental factors such as nutrients, oxygen, 
        and chemotherapy concentration. Cells can divide if conditions are suitable.

        Parameters
        ----------
        chemo_concentration : float, optional
            The chemotherapy concentration in the extracellular matrix. If None, the concentration will 
            be reduced by a default amount.

        Returns
        -------
        None
        '''
        # Set or reduce chemotherapy concentration in the extracellular matrix.
        if chemo_concentration is not None:
            self.extracellular_matrix.set_chemo_concentration(chemo_concentration)
        else:
            self.extracellular_matrix.reduce_chemo_concentration(18)
        
        # Prepare an empty matrix for the next step of cell states.
        next_step_cells = np.full((self.size, self.size), None, dtype=object)
        
        # Shuffle the order of cell locations for simulating step-wise updates.
        sample_indices = np.random.choice(len(self.all_indices), len(self.all_indices), replace=False)
        
        # Process each cell in the randomized order.
        for cell_index in sample_indices:
            cell_location = self.all_indices[cell_index]
            cell_row, cell_column = cell_location
        
            # Create a copy of the current cell to avoid modifying the original matrix during the simulation.
            cell = copy.deepcopy(self.cell_matrix[cell_row, cell_column])
            
            if cell is None:
                continue  # Skip empty locations.
            
            # Check the neighboring cells of the current cell.
            cell_neighbors = cell.neighbors
            is_empty_neighbor = False
            neighbor_status = []
            
            # Check if there are empty neighboring spots for cell division.
            for neighbor in cell_neighbors:
                if self.cell_matrix[neighbor[0], neighbor[1]] == None:
                    is_empty_neighbor = True
                    neighbor_status.append(True)
                else:
                    neighbor_status.append(False)
                    
            # Get environmental factors like nutrient, oxygen, and chemotherapy concentration at the cell's location.
            environment_nutrient = self.extracellular_matrix.ecm[cell_row, cell_column, 0]
            environment_oxygen = self.extracellular_matrix.ecm[cell_row, cell_column, 1]
            environment_chemo = self.extracellular_matrix.ecm[cell_row, cell_column, 2]
            
            # Handle cell death (apoptosis) based on stage counter.
            if cell.is_dying and cell.stage_counter == 12:
                next_step_cells[cell_row, cell_column] = None
            else:
                # Update the cell's state (whether it divides or not).
                is_dividing = cell.update_cell(environment_nutrient, environment_oxygen, environment_chemo, is_empty_neighbor)
                
                # If the cell divides and there is an empty neighboring spot, replace it.
                if is_dividing and is_empty_neighbor:
                    replacement_index = np.random.choice(np.where(np.array(neighbor_status))[0])
                    replacement_location = cell_neighbors[replacement_index]
                    cell_atributes = cell.get_cell_attributes()
                    next_step_cells[replacement_location[0], replacement_location[1]] = Cell(replacement_location[0], replacement_location[1], self.size, cell_atributes[6], cell_atributes[:6])
                
                # If the cell divides but there are no empty neighbors, still replace with a random neighbor.
                elif is_dividing and is_empty_neighbor == False:
                    replacement_index = np.random.choice(range(len(cell_neighbors)))
                    replacement_location = cell_neighbors[replacement_index]
                    cell_atributes = cell.get_cell_attributes()
                    next_step_cells[replacement_location[0], replacement_location[1]] = Cell(replacement_location[0], replacement_location[1], self.size, cell_atributes[6], cell_atributes[:6])
                
                # Update the cell's state in the next step matrix.
                next_step_cells[cell_row, cell_column] = cell
        
        # After processing all cells, update the tissue's cell matrix with the next step.
        self.cell_matrix = next_step_cells
    
    def generate_plot_array(self):
        '''
        Generates a 2D array representing the tissue, where 0 indicates an empty cell, 1 indicates a normal 
        cell, and 2 indicates a cancer cell.

        Parameters
        ----------
        None

        Returns
        -------
        np.ndarray
            A 2D array representing the tissue structure for plotting purposes.
        '''
        plotting_array = np.full((self.size, self.size), 0, dtype=int)
        for row in range(self.size):
            for column in range(self.size):
                if self.cell_matrix[row,column] != None:
                    # If the cell is a cancer cell, mark it with 2.
                    if self.cell_matrix[row,column].is_cancer:
                        plotting_array[row, column] = 2
                    # Otherwise, mark it as a normal cell (1).
                    else:
                        plotting_array[row, column] = 1
        return plotting_array
    
    def place_cancer_cell(self, row, column):
        '''
        Places a cancer cell at a specific location in the tissue grid.

        Parameters
        ----------
        row : int
            The row index of the cell in the grid.
        column : int
            The column index of the cell in the grid.

        Returns
        -------
        None
        '''
        self.cell_matrix[row, column] = Cell(row, column, self.size, True)
    
    def get_metric(self, metric="Number of Cancer"):
        '''
        Retrieves specific metrics about the cells in the tissue to measure attributes over time.

        Parameters
        ----------
        metric : str, optional
            The metric to retrieve. Possible values include:
                - "Number of Cancer": Counts the number of cancer cells in the tissue.
                - "Number of Mutations": Counts the number of mutations in each cell.
                - "Apoptosis Resistance": Measures each cell's resistance to apoptosis.
                - "Recovery Rate": Measures each cell's recovery rate.
                - "Mutation Probability": Measures each cell's mutation probability.
                - "Invasion Probability": Measures each cell's invasion probability.
                - "Damage": Measures the damage in each cell.

        Returns
        -------
        np.ndarray
            A 2D array representing the selected metric for each cell.
        '''
        metric_array = np.full((self.size, self.size), 0, dtype=float)
        
        # Depending on the metric, update the metric_array for each cell.
        if metric == "Number of Cancer":
            for row in range(self.size):
                for column in range(self.size):
                    if self.cell_matrix[row,column] != None:
                        # Count cancer cells as 1.
                        if self.cell_matrix[row,column].is_cancer:
                            metric_array[row, column] = 1
                        else:
                            metric_array[row, column] = 0
                            
        elif metric == "Number of Mutations":
            for row in range(self.size):
                for column in range(self.size):
                    if self.cell_matrix[row,column] != None:
                        metric_array[row, column] = self.cell_matrix[row, column].num_mutation
        
        elif metric == "Apoptosis Resistance":
            for row in range(self.size):
                for column in range(self.size):
                    if self.cell_matrix[row,column] != None:
                        metric_array[row, column] = self.cell_matrix[row, column].resistance_apoptosis
        
        elif metric == "Recovery Rate":
            for row in range(self.size):
                for column in range(self.size):
                    if self.cell_matrix[row,column] != None:
                        metric_array[row, column] = self.cell_matrix[row, column].recovery_rate
        
        elif metric == "Mutation Probability":
            for row in range(self.size):
                for column in range(self.size):
                    if self.cell_matrix[row,column] != None:
                        metric_array[row, column] = self.cell_matrix[row, column].mutation_probability
        
        elif metric == "Invasion Probability":
            for row in range(self.size):
                for column in range(self.size):
                    if self.cell_matrix[row,column] != None:
                        metric_array[row, column] = self.cell_matrix[row, column].invasion_prob
        
        elif metric == "Damage":
            for row in range(self.size):
                for column in range(self.size):
                    if self.cell_matrix[row,column] != None:
                        metric_array[row, column] = self.cell_matrix[row, column].damage
        
        return metric_array


# Define parameters for the test simulation
file_name = "Test"  # Name of the file to save the results
chemo_admin_interval_days = 28  # Interval in days for chemotherapy administration
chemo_admin_interval = chemo_admin_interval_days * 24  # Convert to hours
chemo_max_concentration = 0.6  # Maximum chemotherapy concentration
treatment_delay_days = 60  # Delay before chemotherapy starts in days
treatment_delay = treatment_delay_days * 24  # Convert to hours

# Define the total duration of the simulation
test_duration_days = 250  # Total duration of the simulation in days
num_steps = test_duration_days * 24  # Number of simulation steps (in hours)
size = 41  # Size of the tissue grid (41x41)
tissue = Tissue(size)  # Initialize the tissue with the given size
tissue.place_cancer_cell(20, 30)  # Place a cancer cell at position (20, 30)

# Set up color map for visualizing the tissue grid (black for normal, salmon for cancer, saddle brown for mutations)
cmap = ListedColormap(["Black", mcolors.CSS4_COLORS['salmon'], mcolors.CSS4_COLORS['saddlebrown']])

# Initialize a list to store the plot frames and a dictionary for metrics
plotting_array = tissue.generate_plot_array()  # Generate the initial plot array of the tissue grid
frames = []  # List to store the frames for animation
metrics = {
    "Number of Cancer": [],
    "Number of Mutations": [],
    "Apoptosis Resistance": [],
    "Recovery Rate": [],
    "Mutation Probability": [],
    "Invasion Probability": [],
    "Damage": []
}  # Dictionary to store various metrics for each step

def get_all_metrics():
    '''
    Collects and stores all metrics for the current state of the tissue.
    
    Parameters
    ----------
    None.
    
    Returns
    -------
    None.
    '''
    for metric in metrics.keys():  # Iterate through each metric
        metrics[metric].append(tissue.get_metric(metric))  # Append the current metric value to the corresponding list

# Start the simulation
start_time = time.time()  # Record the start time for performance measurement

for step in range(num_steps):  # Loop through the simulation steps (each step represents one hour)
    
    # Administer chemotherapy at specified intervals if the step is beyond the treatment delay
    if step >= treatment_delay and (step - treatment_delay) % chemo_admin_interval == 0:
        print("Delivered")  # Print message indicating chemotherapy was delivered
        tissue.simulate_step(chemo_max_concentration)  # Simulate the step with maximum chemotherapy concentration
    else:
        tissue.simulate_step()  # Simulate the step without chemotherapy
    
    # Every 12 hours (0.5 days), generate a frame for the plot and record metrics
    if step % 12 == 0:
        print(step)  # Print the current step number (every 12 hours)
        plotting_array = tissue.generate_plot_array()  # Generate a new plot array for the current state
        frames.append(plotting_array)  # Append the generated frame to the list
        
        get_all_metrics()  # Collect and store all metrics for the current step
        
        # If there are no cancer cells left in the tissue, print a message and break the loop
        if np.sum(tissue.get_metric()) == 0:
            print("Cancer Eradicated")  # Print message indicating cancer has been eradicated
            break

# End the simulation and measure the time taken
end_time = time.time()
print(f"Time to simulate {test_duration_days} days: {end_time - start_time}")  # Print the total time taken for the simulation

# Save the collected metrics to a file using pickle
with open(f'{file_name}.pkl', 'wb') as f:
    pickle.dump(metrics, f)  # Serialize the metrics dictionary and save it to a .pkl file

# Create an animation of the tissue over time
fig, ax = plt.subplots()  # Create a new figure and axis for the plot
img_anim = ax.imshow(frames[0], cmap=cmap)  # Display the first frame in the animation
title = ax.set_title("Day 0.00")  # Set the initial title of the plot

def update(frame_idx):
    '''
    Update the frame in the animation.
    
    Parameters
    ----------
    frame_idx : int
        The index of the current frame to display.
    
    Returns
    -------
    list
        The list of artists to update in the animation (image and title).
    '''
    img_anim.set_data(frames[frame_idx])  # Update the image data with the current frame
    days = frame_idx * 12 / 24  # Calculate the corresponding day for the current frame
    title.set_text(f"Day {days:.2f}")  # Update the title with the current day
    return [img_anim, title]  # Return the artists to update

# Create the animation
ani = FuncAnimation(fig, update, frames=len(frames), interval=500)  # Create the animation with the specified frame update interval

# Save the animation as a video file
ani.save(f"{file_name}.mp4", fps=15, dpi=150, writer='ffmpeg')  # Save the animation as an MP4 video


