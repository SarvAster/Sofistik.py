#+============================================================================
#| Company:   PINI Group
#|      pini_functions.py
#| Requires the file "sofistik_daten.py", obtainable on Sofistik website
#+============================================================================

import os 
import subprocess
from ctypes import * 
from sofistik_daten import * 
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from shapely.geometry import *

# Set DLL dir path
os.add_dll_directory(r"C:\Program Files\SOFiSTiK\2024\SOFiSTiK 2024")

class CDBData:
    
    def __init__(self, library_path="sof_cdb_w-2024.dll"):
        """
        Initializes the CDB manager with the path to the DLL library.
        
        :param library_path: Path to the DLL library for CDB management.
        """
        self.myDLL = cdll.LoadLibrary(library_path)
        self.cdbStat = None
        self.Index = None

    def open_cdb(self, file_name, cdb_index=99):
        """
        Opens the specified CDB file.
        
        :param file_name: Path to the CDB file.
        :param cdb_index: CDB index (default: 99).
        """
        self.Index = c_int()
        self.Index.value = self.myDLL.sof_cdb_init(file_name.encode('utf8'), cdb_index)
        self.cdbStat = c_int()
        self.cdbStat.value = self.myDLL.sof_cdb_status(self.Index.value)
        print("CDB opened successfully, CDB Status =", self.cdbStat.value)

    def close_cdb(self):
        """
        Closes the CDB.
        """
        self.myDLL.sof_cdb_close(0)
        self.cdbStat.value = self.myDLL.sof_cdb_status(self.Index.value)
        if self.cdbStat.value == 0:
            print("CDB closed successfully, CDB Status = 0")

    def extract_quad_forces(self, start_key=1, max_attempts=11):
        """
        Extracts the forces of quadrilateral elements from the CDB.
        
        :param start_key: Initial key of the group to be read (default: 1).
        :param max_attempts: Maximum number of attempts to increment the key (default: 11).

        :return: List of forces for each quadrilateral element.
        """
        forces_data = []  # List to store the forces of each element
        ie = c_int(0)
        RecLen = c_int(sizeof(CQUAD_FOR))  # Size of the forces structure
        forces = CQUAD_FOR()  # Instance of the forces structure
        current_key = start_key  # Current key initialized to start_key

        # Attempt to read the forces by incrementing the key up to max_attempts
        for attempt in range(1, max_attempts):
            while ie.value < 2:
                ie.value = self.myDLL.sof_cdb_get(self.Index, 210, current_key, byref(forces), byref(RecLen), 1)

                if ie.value == 0:  # No error, data has been read
                    forces_data.append({
                        "quad_number": forces.m_nr,  # Quadrilateral number
                        "mxx": forces.m_mxx,        # Bending moment mxx
                        "myy": forces.m_myy,        # Bending moment myy
                        "mxy": forces.m_mxy,        # Bending moment mxy
                        "vx": forces.m_vx,          # Shear force vx
                        "vy": forces.m_vy,          # Shear force vy
                        "nx": forces.m_nx,          # Membrane force nx
                        "ny": forces.m_ny,          # Membrane force ny
                        "nxy": forces.m_nxy         # Membrane force nxy
                    })

                # Reset RecLen before the next sof_cdb_get call
                RecLen = c_int(sizeof(CQUAD_FOR))

            # If data has been found, stop searching
            if forces_data:
                break
            else:
                # Reset ie for the next key and increment the key
                ie = c_int(0)
                current_key += 1

        # Check if data was found after all attempts
        if not forces_data:
            print(f"No data found after {max_attempts} attempts starting from key 210/{start_key}.")

        return forces_data  # Return the list with all force data

    def extract_quad_membran_forces(self,start_key=1,max_attempts=11):
        """
        Extracts the membrane forces of quadrilateral elements from the CDB.

        :param start_key: Initial key of the group to be read (default: 1).
        :param max_attempts: Maximum number of attempts to increment the key (default: 11).

        :return: DataFrame with membrane forces only.
        """
        df_forces=pd.DataFrame(self.extract_quad_forces(start_key,max_attempts))
        membran_forces = df_forces[['quad_number', 'nx', 'ny', 'nxy']] #Select the membran forces
        membran_forces = membran_forces.iloc[2:].reset_index(drop=True)  # Delete first 2 rows that give useless data
        return membran_forces
    
    def extract_node_coords(self):
        """
        Extracts node coordinates from the CDB.

        :return: Dictionary with node coordinates.
        """
        ie = c_int(0)
        RecLen = c_int(sizeof(CNODE))  # Size of the CNODE structure
        node = CNODE()  # Instance of the CNODE structure
        node_coords = {}  # Dictionary to store node coordinates

        while ie.value < 2:
            ie.value = self.myDLL.sof_cdb_get(self.Index, 20, 0, byref(node), byref(RecLen), 1)  # Read nodes from CDB
            if ie.value == 0:  # No error, node read correctly
                node_coords[node.m_nr] = (node.m_xyz[0], node.m_xyz[1], node.m_xyz[2])  # Store X, Y, Z coordinates of the node
            RecLen = c_int(sizeof(CNODE))  # Reset the record size
        return node_coords  # Return node coordinates

    def extract_quad_elements(self):
        """
        Extracts quadrilateral elements (with associated nodes) from the CDB.

        :return: Tuple with list of nodes and list of element numbers.
        """
        ie = c_int(0)
        RecLen = c_int(sizeof(CQUAD))  # Size of the CQUAD structure
        quad_element = CQUAD()  # Instance of the CQUAD structure
        elements = []  # List to store quadrilateral elements
        element_numbers = []  # List to store element numbers

        while ie.value < 2:
            ie.value = self.myDLL.sof_cdb_get(self.Index, 200, 0, byref(quad_element), byref(RecLen), 1)  # Read quadrilateral elements
            if ie.value == 0:  # No error, element read correctly
                elements.append(quad_element.m_node[:])  # Add the 4 nodes of the element to the list
                element_numbers.append(quad_element.m_nr)
            RecLen = c_int(sizeof(CQUAD))  # Reset the record size

        return elements, element_numbers  # Return the nodes and the number of each element
    
    def extract_structure_bounds(self):
        """
        Extracts the bounds of the structure from the CDB.

        :return: List of bounds of this form: [xmin,ymin,zmin,xmax,ymax,zmax].
        """
        syst = CSYST()
        RecLen = c_int(sizeof(CSYST))
        ie = c_int(0)
        x = []
        while ie.value < 2:
            # Read SYST data (10/00)
            ie.value = self.myDLL.sof_cdb_get(self.Index, 10, 00, byref(syst), byref(RecLen), 1)
            if ie.value == 0:
                # Extract values from m_box (2x3)
                for i in range(3):
                    for j in range(2):
                        x.append(syst.m_box[i][j])
                return x
            
    def extract_fctm_conc(self):
        """
        Extracts average maximum traction force needed for the concrete to break.

        :return: float number.
        """
        ie = c_int(0)
        mat_conc = CMAT_CONC()
        RecLen = c_int(sizeof(CMAT_CONC))

        while ie.value < 2:
            ie.value = self.myDLL.sof_cdb_get(self.Index, 1, 1, byref(mat_conc), byref(RecLen), 2)
            if ie.value == 0:  # No error, node read correctly
                if mat_conc.m_id == 1:
                    fctm = mat_conc.m_ftm
            # Always read the length of the record before sof_cdb_get is called
            RecLen = c_int(sizeof(mat_conc))

        return fctm / 1000  # Conversion to MPa
    
    def height(self, bounds, direction='y'):
        """
        Returns the height of the plane according to the specified direction ('y' by default).
        
        :param bounds: List of bounds of the structure [xmin, ymin, zmin, xmax, ymax, zmax]
        :param direction: 'y' or 'z' to indicate the direction of height calculation
        :return: Height (ymax - ymin if direction is 'y', zmax - zmin otherwise)
        """
        _, ymin, zmin = bounds[:3]  # First row (min values)
        _, ymax, zmax = bounds[3:]  # Second row (max values)
        
        # Calculate height based on the direction
        if direction == 'y':
            height = ymax - ymin
        else:
            height = zmax - zmin

        return height
    
    def count_meshes_height(self, node_coords,bounds,direction='y'):
        """
        Counts the number of meshes along the height according to the specified direction (y or z).

        :param node_coords: List of node coordinates
        :param bounds: List of bounds of the structure [xmin, ymin, zmin, xmax, ymax, zmax]
        :param direction: Plane direction ('y' for XY, 'z' for XZ)
        :return: Number of meshes along the height 
        """

        if direction == 'y':
            # XY plane: find all nodes at `xmin`
            side_nodes = [node for node in node_coords if np.isclose(node_coords[node][0], bounds[0])]  # bounds[0] = x_min
            # Extract Y coordinates of nodes located at `xmin`
            side_coords = [node_coords[node][1] for node in side_nodes]  # Y-coordinates

        elif direction == 'z':
            # XZ plane: find all nodes at `xmin`
            side_nodes = [node for node in node_coords if np.isclose(node_coords[node][0], bounds[0])]  # bounds[0] = x_min
            # Extract Z coordinates of nodes located at `xmin`
            side_coords = [node_coords[node][2] for node in side_nodes]  # Z-coordinates
        else:
            raise ValueError("The direction must be 'y' for the XY plane or 'z' for the XZ plane.")

        # Sort the coordinates and remove duplicates
        side_coords_sorted = sorted(set(side_coords))

        # The number of meshes is the number of distinct intervals
        number_of_meshes = len(side_coords_sorted) - 1  # Number of meshes = intervals between distinct coordinates

        return number_of_meshes

    def average_mesh_height(self,direction='y'):
            return self.height(direction) / self.count_meshes_height(direction)
    
    def tensile_threshold(self, node_coords,bounds,fctm, direction='y', e=1):
        """
        Calculates the tensile force threshold. 

        :param node_coords: List of node coordinates
        :param bounds: List of bounds of the structure [xmin, ymin, zmin, xmax, ymax, zmax]
        :param fctm: Average maximum traction force needed for rupture.
        :param direction: Plane direction ('y' for XY, 'z' for XZ)
        :param e: Thickness of a mesh element (1m in 2D)
        :return: Force threshold F in kN
        """

        height = self.height(bounds,direction)
        number_of_meshes=self.count_meshes_height(node_coords,bounds,direction)
        h_mesh = height/number_of_meshes
        F = fctm * h_mesh * e

        return F
    

class Ties:
    
    def __init__(self):
        """
        Initializes the Ties class with empty lists for ties and forces.
        """
        self.ties = []  # List to store tie coordinates
        self.forces = []  # List to store forces acting on the ties
        self.direction = 'y'  # Default direction of the ties
        self.spacing = (5, 5)  # Default spacing in X and Y/Z directions

    def add_structure_ties(self, bounds, dir='y', n=5, m=5):
        """
        Generates a list of ties based on the structure bounds.
        
        WARNING: Only applies to rectangular structures.
        
        :param bounds: List of bounds of the structure [xmin, ymin, zmin, xmax, ymax, zmax]
        :param dir: Direction of the ties ('y' for horizontal/vertical in XY plane, 'z' for XZ plane)
        :param n: Number of ties to add in X direction
        :param m: Number of ties to add in Y or Z direction
        :return: List of added ties
        """
        # Get the vertical and horizontal bounds
        x_min, y_min, z_min = bounds[:3]
        x_max, y_max, z_max = bounds[3:]

        # Add margin
        margin = 0.1
        x_min, y_min, z_min = x_min + margin, y_min + margin, z_min + margin
        x_max, y_max, z_max = x_max - margin, y_max - margin, z_max - margin

        # Calculate tie spacing
        tie_spacing_x = (x_max - x_min) / n
        tie_spacing_y = (y_max - y_min) / m
        tie_spacing_z = (z_max - z_min) / m

        self.direction = dir
        self.spacing = (n, m)

        if dir == 'y':
            # Horizontal ties    
            for i in range(m+1):
                y = y_min + i * tie_spacing_y
                self.ties.append(((x_min, y), (x_max, y)))
            # Vertical ties
            for i in range(n+1):
                x = x_min + i * tie_spacing_x
                self.ties.append(((x, y_min), (x, y_max)))        
        else:
            # Horizontal ties
            for i in range(m+1):
                z = z_min + i * tie_spacing_z
                self.ties.append(((x_min, z), (x_max, z)))
            # Vertical ties
            for i in range(n+1):
                x = x_min + i * tie_spacing_x
                self.ties.append(((x, z_min), (x, z_max)))
        
        return self.ties
    
    def elements_crossed_by_ties(self, node_coords, elements):
        """
        Finds the elements crossed by each tie.
        
        :param node_coords: List of mesh node coordinates
        :param elements: List of mesh elements (each element is a list of nodes)
        :param dir: Direction for the intersection calculation
        :return: Dictionary with ties as keys and list of crossed elements as values
        """
        crossed_elements = {}  # Store the sub-ties and their crossed elements
        
        for tie in self.ties:

            tie_line = LineString([tie[0], tie[1]])
            min_x, min_y, max_x, max_y = tie_line.bounds 

            crossed_elements_tie = []

            for element in elements:
                if self.direction == 'y':
                    quad_coords = np.array([node_coords[node][:2] for node in element])
                else:
                    quad_coords = np.array([(node_coords[node][0], node_coords[node][2]) for node in element])

                quad_polygon = Polygon(quad_coords)
                if tie_line.intersects(quad_polygon):
                    crossed_elements_tie.append(element)  # Add the crossed element

            crossed_elements[tie] = crossed_elements_tie  # Associate the crossed elements with this tie

        return crossed_elements
    
    def membran_forces_in_ties(self, crossed_elements, elements, element_numbers, membran_forces):
        """
        Calculates the average force for each sub-tie considering the orientation of the tie.

        :param crossed_elements: Elements crossed by each tie
        :param elements: List of mesh elements (each element is a list of nodes)
        :param element_numbers: List of mesh elements'numbers
        :param membran_forces: Dataframe with membran_forces for each element
        :return: Forces for each tie
        """
        for tie in self.ties:
            if crossed_elements[tie]:

                total_force=0
                count=0
                
                # Determine the orientation of the tie: horizontal (nx) or vertical (ny)
                (x1, y1), (x2, y2) = tie
                if abs(y1 - y2) < 1e-5:  # Horizontal tie (in X)
                    direction = 'nx'
                elif abs(x1 - x2) < 1e-5:  # Vertical tie (in Y or Z)
                    direction = 'ny'
                else:
                    # For inclined ties, choose a default direction or add extra logic
                    direction = 'nx'  # Example: choose 'nx' by default, or adjust as needed

                # Sum the forces for all elements crossed by this sub-tie
                for element in crossed_elements[tie]:

                    # Get the element's number for membran_forces
                    index=elements.index(element)
                    element_number=element_numbers[index] 

                    # Filter the DataFrame to retrieve the row corresponding to 'quad_number'
                    element_forces = membran_forces[membran_forces['quad_number'] == element_number]
                    # Retrieve the force in the specified direction
                    if direction == 'nx':
                        total_force+=element_forces['nx'].values[0]
                        count+=1
                    else:
                        total_force+=element_forces['ny'].values[0]
                        count+=1
                        
                # Add the force corresponding to the tie
                self.forces.append(total_force/count)
            else:
                self.forces.append(0)
                
        return self.forces
    
    def remove_ties_above_threshold(self, threshold):
        """
        Removes ties which average force is above the threshold.
        
        :param threshold: Threshold
        :return: List of retained ties
        """
        new_ties=[]
        new_forces=[]
        for i in range(len(self.ties)):
            if self.forces[i] >= threshold:  # Keep only ties with forces below the threshold
                new_ties.append(self.ties[i])
                new_forces.append(self.forces[i])

        self.ties = new_ties
        self.forces = new_forces

    def generate_reinforcement_ties_from_concrete(self, node_coords, elements, element_numbers, bounds, membran_forces, threshold, dir='y', n=5, m=5):
        """
        Generate a ties' structure based on the structure of the concrete.
        
        :param node_coords: List of node coordinates
        :param elements : List of elements
        :param element_numbers: List of mesh elements'numbers
        :param bounds: List of bounds of the structure [xmin, ymin, zmin, xmax, ymax, zmax]
        :param threshold: Threshold
        :param membran_forces: Dataframe with membran_forces for each element
        :param dir: Indicates the orientation of the reinforcement ('y' or 'z')
        :param n: Number of ties to add in the X direction
        :param m: Number of ties to add in the Y or Z direction
        :return: List of added ties
        """

        # Initially add a full frame
        self.add_structure_ties(bounds,dir,n,m) 

        # Find the crossed elements by every tie
        crossed_elements=self.elements_crossed_by_ties(node_coords,elements)

        # Update the forces for each tie
        self.membran_forces_in_ties(crossed_elements,elements,element_numbers,membran_forces)
    
        # Remove the ties below the threshold
        self.remove_ties_above_threshold(threshold)
    
        return self.ties
    
class Plot:
    def __init__(self):
        """
        Initializes the Plot class as an empty Plotly figure.
        """
        self.fig = go.Figure()  # Create an empty Plotly figure

    def add_fig_ties(self, ties):
        """
        Adds ties (green lines) to the figure based on the extremities extracted from the CDB.

        :param ties: Ties to add to the figure
        :param dir: Ties direction ('y' for the XY plane, 'z' for the XZ plane)
        :param n: Number of ties in the X direction
        :param m: Number of ties in the Y or Z direction
        :return: Updated figure and list of ties with coordinates of the extremities
        """

        # Add ties to the figure
        for tie in ties:
            # Extract coordinates of tie extremities
            (x1, y1), (x2, y2) = tie
            # Add the green line to the figure
            self.fig.add_trace(go.Scatter(
                x=[x1, x2], y=[y1, y2], mode='lines',
                line=dict(color='green', width=2, dash='dash'),
                showlegend=False
            ))

        # Return the updated figure 
        return self.fig
    
    def plot_structure_with_ties(self, node_coords,elements,ties, view_axis='xy'):
        """
        Plots the structure with green lines (ties) and returns the updated figure with ties.

        :param node_coords: List of node coordinates
        :param elements : List of elements
        :param ties : List of ties
        :param view_axis: Selected view ('xy', 'xz', 'yz')
        :param n: Number of ties to add in the X direction
        :param m: Number of ties to add in the Y or Z direction
        :return: Updated Plotly figure with ties and list of ties
        """
        
        # Define axes based on the chosen view
        axis_map = {
            'xy': (0, 1),  # Use X and Y axes
            'xz': (0, 2),  # Use X and Z axes
            'yz': (1, 2)   # Use Y and Z axes
        }
        ax1, ax2 = axis_map.get(view_axis, (0, 1))  # Default to XY view

        # Plot the mesh without values
        for element in elements:
            quad_coords = np.array([node_coords[node] for node in element])
            quad_coords = np.vstack([quad_coords, quad_coords[0]])  # Close the quadrilateral
            x_coords = quad_coords[:, ax1]
            y_coords = quad_coords[:, ax2]

            # Add the quadrilateral to the plot
            self.fig.add_trace(go.Scatter(
                x=x_coords, y=y_coords, mode='lines',
                line=dict(color='black', width=1),
                showlegend=False
            ))

        # Use the add_fig_ties function to add the ties to the figure
        self.fig = self.add_fig_ties(ties)

        # Configure layout for interactivity
        axis_labels = {
            'xy': ('X', 'Y'),
            'xz': ('X', 'Z'),
            'yz': ('Y', 'Z')
        }
        x_label, y_label = axis_labels.get(view_axis, ('X', 'Y'))

        self.fig.update_layout(
            title=f"Structure with Ties (View: {view_axis})",
            xaxis_title=x_label,
            yaxis_title=y_label,
            xaxis=dict(scaleanchor="y", scaleratio=1),  # Force equal scale
            yaxis=dict(scaleanchor="x", scaleratio=1),
            showlegend=False,
            autosize=True,
        )

        # Return the updated figure and list of ties
        return self.fig
    
        # Function to plot the structure, hatch crossed elements, and display ties
    def plot_structure_with_hachured_elements(self, node_coords, elements, crossed_elements, ties, view_axis='xy'):
        """
        Plots the structure, hatches the crossed elements in blue, and displays the ties.
        
        :param node_coords: List of mesh node coordinates
        :param elements: List of elements 
        :param crossed_elements: Dictionary with ties and their respective crossed elements
        :param ties: List of ties (green lines)
        :param view_axis: Visualization plane ('xy', 'xz', 'yz')
        """


        # Define axes based on the chosen view
        axis_map = {'xy': (0, 1), 'xz': (0, 2), 'yz': (1, 2)}
        ax1, ax2 = axis_map.get(view_axis, (0, 1))

        # Extract all crossed elements from the dictionary into a set for faster lookup
        all_crossed_elements = set()
        for elements_list in crossed_elements.values():
            # Convert each element list to a tuple to make it hashable
            all_crossed_elements.update(tuple(sorted(element)) for element in elements_list)
        # Plot the mesh with black edges and hatch crossed elements
        for element in elements:
            quad_coords = np.array([node_coords[node] for node in element])
            quad_coords = np.vstack([quad_coords, quad_coords[0]])  # Close the quadrilateral
            x_coords = quad_coords[:, ax1]
            y_coords = quad_coords[:, ax2]

            # Draw edges in black
            self.fig.add_trace(go.Scatter(
                x=x_coords, y=y_coords, mode='lines',
                line=dict(color='black', width=2),
                showlegend=False
            ))

            # If the element is crossed, hatch it in blue
            if tuple(sorted(element)) in all_crossed_elements:
                # Fill the element with a translucent blue color
                self.fig.add_trace(go.Scatter(
                    x=x_coords, y=y_coords, mode='lines', fill='toself',
                    fillcolor='rgba(0, 0, 255, 0.2)',  # Translucent blue
                    line=dict(color='black', width=2),
                    showlegend=False
                ))

        # Add green ties to the figure
        for tie in ties:
            tie_coords = np.array(tie)
            x_tie = tie_coords[:, 0]
            y_tie = tie_coords[:, 1]

            self.fig.add_trace(go.Scatter(
                x=x_tie, y=y_tie, mode='lines',
                line=dict(color='green', width=2, dash='dash'),
                showlegend=False
            ))

        # Configure layout for interactivity
        self.fig.update_layout(
            title="Structure with Hachured Elements and Ties",
            xaxis_title='X',
            yaxis_title='Y',
            xaxis=dict(scaleanchor="y", scaleratio=1),
            yaxis=dict(scaleanchor="x", scaleratio=1),
            showlegend=False,
            autosize=True,
        )

        return self.fig

class SofiFileHandler:
    def __init__(self):
        """
        Initializes the SofiFileHandler class with the path to the .dat file.

        :param dat_file_path: Path to the .dat file to be modified
        """
        self.dat_file_path = None
        self.cdb_file_path = None

    def add_dat(self, dat_file_path):
        """
        Sets the path to the .dat file.

        :param dat_file_path: Path to the .dat file to be modified
        """
        self.dat_file_path = dat_file_path
        print(f".dat file path set to: {self.dat_file_path}")

    def add_cdb(self, cdb_file_path):
        """
        Sets the path to the .cdb file.

        :param cdb_file_path: Path to the .cdb file
        """
        self.cdb_file_path = cdb_file_path

    def calculate_with_sps(self):
        """
        Executes the calculation of the current .dat file using SOFiSTiK in batch mode via sps.exe.
        """
        if not self.dat_file_path:
            print("Error: .dat file path is not set. Use add_dat() to set the file path.")
            return

        try:
            # Command to run sps.exe with the specified .dat file
            sps_command = r"C:\Program Files\SOFiSTiK\2024\SOFiSTiK 2024\sps.exe" + f' "{self.dat_file_path}"'

            # Launch sps.exe with the .dat file and wait for it to complete
            process = subprocess.Popen(sps_command)

            # Wait for the process to complete
            process.wait()

            # Check if the process finished successfully
            if process.returncode == 0:
                print("Calculation successfully completed in SOFiSTiK.")
            else:
                print(f"Calculation failed with exit code {process.returncode}.")

        except FileNotFoundError:
            print("Error: sps.exe not found. Ensure that SOFiSTiK is installed and sps.exe is accessible in your system's PATH.")

        except Exception as e:
            print(f"Error during SOFiSTiK execution: {e}")

    def add_ties(self,ties,dir='y'):
        """
        Modifies the .dat file by adding SLN and SLNB lines for each tie,
        with a unique marker to facilitate future removal.

        :param ties: List of ties with their coordinates (each tie is a tuple of two points)
        """
        # Read the existing .dat file
        with open(self.dat_file_path, 'r') as file:
            lines = file.readlines()

        last_sln_number = 0
        end_sofimshc_index = None

        # Iterate through the lines to identify the PROG SOFIMSHC section and find the last SLN
        in_sofimshc_section = False
        for i, line in enumerate(lines):
            if "PROG SOFIMSHC" in line:
                in_sofimshc_section = True
            elif "END" in line and in_sofimshc_section:
                end_sofimshc_index = i
                break  # End of the SOFIMSHC section found
            elif line.startswith("SLN") and in_sofimshc_section:
                # Extract the number after SLN
                parts = line.split()
                if len(parts) > 1 and parts[1].isdigit():
                    last_sln_number = int(parts[1])

        # Check that the END of SOFIMSHC was found
        if end_sofimshc_index is None:
            print("Error: End of SOFIMSHC section (END) not found.")
            return

        # Insert SLN and SLNB lines for each sub-tirant
        new_lines = []
        for tie in ties:
            last_sln_number += 1
            point1, point2 = tie  # The two points of the tirant: (x1, z1), (x2, z2)

            if dir=='y':
                # Create the SLN line with a marker !ADDED
                new_sln_line = f"SLN       {last_sln_number}  GRP 1 STYP 'CE' SNO 1 DRX 0 -1 0 TITL \"Line\" !ADDED\n"
                new_lines.append(new_sln_line)

                # Create the SLNB line with the coordinates of the two endpoints and a marker !ADDED
                x1, y1 = point1
                x2, y2 = point2
                new_slnb_line = f"SLNB X1   {x1:.6f} {y1:.6f} 0.000000 X2   {x2:.6f} {y2:.6f} 0.000000  !ADDED\n"
                new_lines.append(new_slnb_line)

            elif dir=='z':
                # Create the SLN line with a marker !ADDED
                new_sln_line = f"SLN       {last_sln_number}  GRP 1 STYP 'CE' SNO 1 DRX 0 0 -1 TITL \"Line\" !ADDED\n"
                new_lines.append(new_sln_line)

                # Create the SLNB line with the coordinates of the two endpoints and a marker !ADDED
                x1, z1 = point1
                x2, z2 = point2
                new_slnb_line = f"SLNB X1   {x1:.6f} 0.000000 {z1:.6f} X2   {x2:.6f} 0.000000 {z2:.6f} !ADDED\n"
                new_lines.append(new_slnb_line)

        # Insert the new lines just before the END of SOFIMSHC
        lines = lines[:end_sofimshc_index] + new_lines + lines[end_sofimshc_index:]

        # Write the modifications to the .dat file
        with open(self.dat_file_path, 'w') as file:
            file.writelines(lines)

        print(f"The ties have been added to the file {self.dat_file_path}.")

    def remove_added_ties(self, ties,dir='y'):
        """
        Removes tirants marked '!ADDED' from the .dat file that are not in the provided list of ties.

        :param ties: List of ties to keep (each tie is a tuple of two points)
        """
        with open(self.dat_file_path, 'r') as file:
            lines = file.readlines()

        # Convert the tirants from the list into a comparable form for SLNB lines
        # Use of a set to optimize efficiency, round(_,6) is used to get 6 decimals.
        ties_set = {((round(x1, 6), round(y1, 6)), 
                            (round(x2, 6), round(y2, 6))) 
                            for x1, y1, x2, y2 in 
                            [(tie[0][0], tie[0][1], tie[1][0], tie[1][1]) 
                            for tie in ties]}

        # Iterate through the lines to detect tirants marked for removal
        lines_to_keep = []
        skip_next = False

        for i, line in enumerate(lines):
            if skip_next:  # Do not add the next line if SLN should be removed
                skip_next = False
                continue

            if line.startswith("SLN") and "!ADDED" in line:
                # The next line is probably the SLNB line, so we check it
                slnb_line = lines[i + 1] if i + 1 < len(lines) else ""
                if slnb_line.startswith("SLNB") and "!ADDED" in slnb_line:
                    # Extract the coordinates from the SLNB
                    parts = slnb_line.split()
                    if len(parts) >= 9:  # Ensure we have enough data

                        if dir=='y':
                            x1, y1 = float(parts[2]), float(parts[3])
                            x2, y2 = float(parts[6]), float(parts[7])

                            # Check if this tie is in the list of ties to keep
                            if ((round(x1, 6), round(y1, 6)), 
                                (round(x2, 6), round(y2, 6))) not in ties_set:
                                # Do not add this SLN and SLNB
                                skip_next = True
                                continue    
                        elif dir=='z':
                            x1, z1 = float(parts[2]), float(parts[4])
                            x2, z2 = float(parts[6]), float(parts[8])

                            # Check if this tie is in the list of ties to keep
                            if ((round(x1, 6), round(z1, 6)), 
                                (round(x2, 6), round(z2, 6))) not in ties_set:
                                # Do not add this SLN and SLNB
                                skip_next = True
                                continue

            # If we reach here, keep the line
            lines_to_keep.append(line)

        # Write the updated lines to the .dat file
        with open(self.dat_file_path, 'w') as file:
            file.writelines

    def add_linear_calculation_block(self, n):
            """
            Adds a linear calculation block after the last 'END' associated with the last 'SOFILOAD' block in the .dat file.

            :param n: Integer used to define the line 'LC n'.
            """

            # Read the content of the .dat file
            with open(self.dat_file_path, 'r') as file:
                lines = file.readlines()

            # Find the position of the last "SOFILOAD" block and its associated "END"
            last_sofiload_index = None
            end_after_sofiload = None
            for i, line in enumerate(lines):
                if "PROG SOFILOAD" in line:
                    last_sofiload_index = i
                    end_after_sofiload = None  # Reset each time a new "SOFILOAD" is found
                elif "END" in line and last_sofiload_index is not None and end_after_sofiload is None:
                    # Record the position of the first "END" after a "SOFILOAD"
                    end_after_sofiload = i

            # Check if a "SOFILOAD" block and its "END" were found
            if last_sofiload_index is None or end_after_sofiload is None:
                print("Error: No complete 'SOFILOAD' block found with an associated 'END'.")
                return

            # Create the linear calculation block to add
            linear_block = [
                "$Linear analysis\n",
                "+PROG ASE urs:9\n",
                "HEAD Calculation of forces and moments\n",
                "CTRL OPT WARP VAL 0\n",
                f"LC {n}\n",
                "END\n"
            ]

            # Insert the block after the "END" associated with the last "SOFILOAD"
            lines = lines[:end_after_sofiload + 1] + linear_block + lines[end_after_sofiload + 1:]

            # Write the modifications to the .dat file
            with open(self.dat_file_path, 'w') as file:
                file.writelines(lines)

            print(f"Linear calculation block successfully added after the last 'END' associated with 'SOFILOAD' in {self.dat_file_path}.")

    def remove_ase_block(self):
        """
        Removes the linear calculation block PROG ASE from the .dat file.
        """

        # Read the content of the .dat file
        with open(self.dat_file_path, 'r') as file:
            lines = file.readlines()

        # Variables to store the start and end indices of the block to be removed
        start_index = None
        end_index = None

        # Search for the start of the "PROG ASE" block and the associated END
        for i, line in enumerate(lines):
            if "PROG ASE" in line and "urs" in line:  # Identify the start of the block
                start_index = i
            elif start_index is not None and "END" in line:  # Find the associated END
                end_index = i
                break  # Exit as soon as the first END after the ASE block is found

        # Check if the block was found and remove the corresponding lines
        if start_index is not None and end_index is not None:
            del lines[start_index:end_index + 1]

            # Write the modifications to the .dat file
            with open(self.dat_file_path, 'w') as file:
                file.writelines(lines)

            print(f"PROG ASE block successfully removed from {self.dat_file_path}.")
            
        else:
            print("Error: PROG ASE block not found or associated END missing.")

    def add_nonlinear_calculation_block(self, n):
        """
        Adds a nonlinear calculation block after the last 'END' associated with the last 'SOFILOAD' block in the .dat file.

        :param n: Integer used to define the line 'LC n'.
        """
        if self.dat_file_path is None:
            print("Error: .dat file path is not set. Use add_dat() to set the file path.")
            return

        # Read the content of the .dat file
        with open(self.dat_file_path, 'r') as file:
            lines = file.readlines()

        # Find the position of the last "SOFILOAD" block and its associated "END"
        last_sofiload_index = None
        end_after_sofiload = None
        for i, line in enumerate(lines):
            if "PROG SOFILOAD" in line:
                last_sofiload_index = i
                end_after_sofiload = None  # Reset each time a new "SOFILOAD" is found
            elif "END" in line and last_sofiload_index is not None and end_after_sofiload is None:
                # Record the position of the first "END" after a "SOFILOAD"
                end_after_sofiload = i

        # Check if a "SOFILOAD" block and its "END" were found
        if last_sofiload_index is None or end_after_sofiload is None:
            print("Error: No complete 'SOFILOAD' block found with an associated 'END'.")
            return

        # Create the nonlinear calculation block to add
        nonlinear_block = [
            "$Non-linear analysis\n",
            "+PROG ASE urs:17.1\n",
            "ECHO FULL FULL\n",
            "CTRL CONC VAL 0.2\n",
            "SYST PROB NONL  ITER 500 NMAT YES\n",
            "NSTR KMOD k1  KSV ULD\n",
            "GRP -\n",
            f"LC {n}\n",
            "END\n"
        ]

        # Insert the block after the "END" associated with the last "SOFILOAD"
        lines = lines[:end_after_sofiload + 1] + nonlinear_block + lines[end_after_sofiload + 1:]

        # Write the modifications to the .dat file
        with open(self.dat_file_path, 'w') as file:
            file.writelines(lines)

        print(f"Nonlinear calculation block successfully added after the last 'END' associated with 'SOFILOAD' in {self.dat_file_path}.")