# Chemotherapy Simulation Agent-Based Model 

## About the Project 
This project is designed to create a flexible framework for simulating cancer growth and the effects of chemotherapy treatments. The simulation is built around three interconnected classes: Cell, Tissue, and Extracellular Matrix (ECM). The model is structured with tunable parameters that reflect real-world biological processes, allowing for detailed exploration of how cancer cells grow, mutate, and respond to treatments like chemotherapy. This model only simulates a 2d cross-section of tissue with a simple vasculature arrangement of a single vessel in the center of the tissue. 

### Overview of the Components
The simulation uses a multilayered network of the previously mentioned classes:

Cell: Represent the individual biological units within the tissue. Each cell is characterized by its position in the tissue matrix, whether it is normal or cancerous, and a set of biological properties such as mutation probability, resistance to death, recovery rate, and invasion probability. The cells interact with each other and their environment and can undergo processes like division or cell death.

Extracellular Matrix (ECM): The ECM is the environment surrounding the cells, and it plays a crucial role in maintaining tissue structure and facilitating cell communication. In the simulation, the ECM is responsible for providing cells with vital nutrients and oxygen, as well as regulating chemo-concentrations. The ECM has a grid-based structure that stores values such as nutrient levels, oxygen availability, and chemotherapy drug concentrations at each tissue location. These values can change over time, impacting cell behavior and survival.

Tissue: The tissue serves as the container for the cells and the ECM. It defines the overall structure of the simulated environment, with a matrix of cells arranged in a grid. The tissue class is responsible for managing the simulation of cell behavior across time steps, updating cell states, and applying chemotherapy treatments when appropriate. The tissue also generates visualizations and metrics that allow users to track cancer progression and treatments' effectiveness. Time steps are one hour for all simulations.

### Key Features and Processes
Cancer Cell Behavior: Cells can be cancerous or normal. Cancer cells can mutate, divide, or die based on several factors, such as their genetic properties, the nutrient and oxygen availability from the ECM, and the presence of chemotherapy drugs. Normal cells also divide and grow but have different characteristics than cancer cells, such as a lower mutation probability and higher susceptibility to cell death.

Chemotherapy Effects: The simulation models the impact of chemotherapy on cancer cells over time. Chemotherapy is introduced at specified intervals, and the concentration of the drug in the ECM affects cell behavior. Chemotherapy is more damaging to cells in the cell cycle, reflecting chemotherapy's actual mechanism. The chemotherapy concentration can be varied, providing flexibility for exploring different treatment schedules and dosages.

Mutation and Cell Death Resistance: Cells can accumulate mutations, which may affect their characteristics. For example, mutations can increase a cell’s ability to survive chemotherapy (apoptosis resistance) or make it more aggressive regarding its ability to invade neighboring tissues. Mutation probabilities and apoptosis resistance are tunable parameters that can be adjusted to simulate different cancer types or genetic mutations.

### Tissue Metrics
The simulation provides multiple metrics that track the state of the tissue over time. These metrics include:

Number of Cancer Cells: Tracks how many cancer cells are present in the tissue at any given time. 
Number of Mutations: Tracks the mutation rate of cells in the tissue. 
Apoptosis Resistance: Tracks how resistant cells are to apoptosis, which affects their ability to survive chemotherapy. 
Recovery Rate: Measures how fast cells recover from damage, impacting their survival and division rates. 
Invasion Probability: Represents the likelihood that cells will invade surrounding tissue, a key factor in cancer metastasis. 
Damage: Tracks how much damage each cell has sustained, reflecting their health and the impact of chemotherapy or other stress factors. 

Visualization: The simulation generates visualizations of the tissue, where each cell is represented by a color based on its state (normal, cancerous, or dead). These visualizations can be animated to show how the tissue changes over time. A series of frames is created throughout the simulation, which can be exported as a video showing the progression of cancer and the effects of chemotherapy treatments. 

### Tunable Parameters
The framework allows the following parameters to be adjusted, enabling users to explore different biological scenarios:

Cell properties: Mutation rate, apoptosis resistance, invasion probability, recovery rate, and the chemotherapy resistance of cancer cells.
Extracellular Matrix: Concentration of nutrients, oxygen, and chemotherapy drugs at different locations in the tissue.  
Tissue Size: The size of the grid representing the tissue affects the number of cells and their interactions.  
Treatment Parameters: The interval and concentration of chemotherapy treatments and the delay before treatment begins.

### Potential Applications
This simulation framework can be used for:

Exploring Cancer Dynamics: Understanding how cancer grows and spreads under different conditions, including the role of mutations and the influence of the ECM.      
Testing Chemotherapy Regimens: Investigating how different chemotherapy schedules, drug concentrations, and administration intervals affect cancer growth and cell survival. 
Investigating Cancer Evolution: Simulating how genetic mutations contribute to cancer progression, including the development of resistance to treatment and metastasis. 
Modeling Tumor Heterogeneity: Understanding how different genetic pathways or cell types coexist and compete within the tumor microenvironment. 

By adjusting the tunable parameters, users can model a variety of cancer types, treatment strategies, and genetic scenarios, making this framework a powerful tool for cancer research and drug development simulations.

Also included in this repository is a brief paper discussing the model in greater depth and an experiment performed using some of the tunable parameters. 
