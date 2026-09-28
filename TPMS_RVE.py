# coding=utf-8
# This program generates 11 types of TPMS cells
from abaqus import *
from abaqusConstants import *
import regionToolset
import mesh
import numpy as np
import math
import os
from math import pi, sin, cos


def schwarz_p(x, y, z):
    return cos(x) + cos(y) + cos(z)


def gyroid(x, y, z):
    return sin(x) * cos(y) + sin(y) * cos(z) + sin(z) * cos(x)


def diamond(x, y, z):
    return sin(x) * sin(y) * sin(z) + sin(x) * cos(y) * cos(z) + cos(x) * sin(y) * cos(z) + cos(x) * cos(y) * sin(z)


def iwp(x, y, z):
    return 2 * (cos(x) * cos(y) + cos(y) * cos(z) + cos(z) * cos(x)) - (cos(2 * x) + cos(2 * y) + cos(2 * z))


def neovius(x, y, z):
    return 3 * (cos(x) + cos(y) + cos(z)) + 4 * cos(x) * cos(y) * cos(z)


def fischer_koch_s(x, y, z):
    return cos(2 * x) * sin(y) * cos(z) + cos(2 * y) * sin(z) * cos(x) + cos(2 * z) * sin(x) * cos(y)


def lidinoid(x, y, z):
    return 0.5 * (sin(2 * x) * cos(y) * sin(z) + sin(2 * y) * cos(z) * sin(x) + sin(2 * z) * cos(x) * sin(y)) - 0.5 * (
                cos(2 * x) * cos(2 * y) + cos(2 * y) * cos(2 * z) + cos(2 * z) * cos(2 * x)) + 0.15


def schoen_frd(x, y, z):
    return 4 * cos(x) * cos(y) * cos(z) - (cos(2 * x) * cos(2 * y) + cos(2 * y) * cos(2 * z) + cos(2 * z) * cos(2 * x))


def split_p(x, y, z):
    return cos(x) * cos(y) + cos(y) * cos(z) + cos(z) * cos(x) - cos(x) * cos(y) * cos(z)


def double_gyroid(x, y, z):
    return 2.0 * (cos(x) * sin(y) + cos(y) * sin(z) + cos(z) * sin(x)) - (cos(2 * x) + cos(2 * y) + cos(2 * z))


def double_diamond(x, y, z):
    return sin(x) * sin(y) * sin(z) + sin(x) * cos(y) * cos(z) + cos(x) * sin(y) * cos(z) + cos(x) * cos(y) * sin(
        z) - 0.4


# Dictionary mapping TPMS names to their functions
TPMS_FUNCTIONS = {
    'Schwarz P': schwarz_p,
    'Gyroid': gyroid,
    'Diamond': diamond,
    'Schoen IWP': iwp,
    'Neovius': neovius,
    'Fischer-Koch S': fischer_koch_s,
    'Lidinoid': lidinoid,
    'Schoen FRD': schoen_frd,
    'Split P': split_p,
    'Double Gyroid': double_gyroid,
    'Double Diamond': double_diamond
}

# Gradient material configuration
POISSON = 0.3


def get_user_input():
    """Get user input through dialog"""
    # Create radio button options for TPMS types
    tpms_types = TPMS_FUNCTIONS.keys()
    tpms_type_options = '\n'.join(tpms_types)

    # Show dialog with fields and parameter ranges
    fields = (
        ('TPMS Type (Choose one from below):\n' + tpms_type_options, 'Gyroid'),
        ('Cell Size (mm) [Range: 0.1-50.0, recommended: 1.0-10.0]:', '1.0'),
        ('Wall Thickness [Range: 0.01-0.8, recommended: 0.1-0.3]:', '0.2'),
        ('Mesh Size [Range: 0.005-0.2, recommended: 0.02-0.1]:', '0.05'),
        #('Enable Gradient Material (YES/NO):', 'YES'),
        #('Gradient Direction (X/Y/Z):', 'Z'),
        #('Number of Gradient Layers [Range: 2-20, recommended: 5-15]:', '10'),
    )

    title = 'TPMS Generator Settings'

    # Show the dialog
    values = getInputs(fields=fields, label='Enter parameters (hover over fields to see ranges):', dialogTitle=title)

    if values is None:
        print('User cancelled the operation')
        return None

    # Parse and validate inputs
    tpms_type = values[0].strip()
    cell_size = float(values[1])
    thickness = float(values[2])
    mesh_size = float(values[3])
    enable_gradient = values[4].strip().upper() == 'YES'
    gradient_direction = values[5].strip().upper()
    num_layers = int(values[6])

    # Validate parameter ranges
    if cell_size < 0.1 or cell_size > 50.0:
        print('Warning: Cell size should be between 0.1 and 50.0 mm (recommended: 1.0-10.0 mm)')
        return None

    if thickness < 0.01 or thickness > 0.8:
        print('Warning: Wall thickness should be between 0.01 and 0.8 (recommended: 0.1-0.3)')
        return None

    if mesh_size < 0.005 or mesh_size > 0.2:
        print('Warning: Mesh size should be between 0.005 and 0.2 (recommended: 0.02-0.1)')
        return None

    if gradient_direction not in ['X', 'Y', 'Z']:
        print('Warning: Gradient direction must be X, Y, or Z')
        return None

    if num_layers < 2 or num_layers > 20:
        print('Warning: Number of layers should be between 2 and 20')
        return None

    # Add warning for potentially problematic parameter combinations
    if mesh_size > cell_size * 0.5:
        print('Warning: Mesh size might be too large relative to cell size, which could result in poor mesh quality')
        retry = getWarningReply(
            message='Mesh size is large relative to cell size. Continue anyway?',
            buttons=(YES, NO)
        )
        if retry == NO:
            return None

    if thickness > cell_size * 0.5:
        print('Warning: Wall thickness might be too large relative to cell size')
        retry = getWarningReply(
            message='Wall thickness is large relative to cell size. Continue anyway?',
            buttons=(YES, NO)
        )
        if retry == NO:
            return None

    if tpms_type not in TPMS_FUNCTIONS:
        print('Invalid TPMS type selected')
        return None

    return {
        'tpms_func': TPMS_FUNCTIONS[tpms_type],
        'cell_size': cell_size,
        'thickness': thickness,
        'mesh_size': mesh_size,
        'tpms_type': tpms_type,
        'enable_gradient': enable_gradient,
        'gradient_direction': gradient_direction,
        'num_layers': num_layers
    }


def generate_voxel_model(tpms_func, cell_size, thickness, mesh_size):
    """Generate TPMS structure using voxel method"""
    # 修改这里：计算体素数量（而不是网格点数量）
    n_voxels = int(cell_size / mesh_size)  # 体素数量
    n_points = n_voxels + 1  # 节点数量 = 体素数量 + 1

    # Generate coordinate grid (scaled to 2π)
    x = np.linspace(0, 2 * pi, n_points)
    y = np.linspace(0, 2 * pi, n_points)
    z = np.linspace(0, 2 * pi, n_points)
    X, Y, Z = np.meshgrid(x, y, z)

    # Calculate TPMS values at each grid point
    tpms_values = np.zeros((n_points, n_points, n_points))
    for i in range(n_points):
        for j in range(n_points):
            for k in range(n_points):
                tpms_values[i, j, k] = tpms_func(X[i, j, k], Y[i, j, k], Z[i, j, k])

    # 修改这里：现在为每个体素（不是每个网格点）分配材料
    # 体素索引从 0 到 n_voxels-1
    tpms_core = np.zeros((n_voxels, n_voxels, n_voxels))
    interface = np.zeros((n_voxels, n_voxels, n_voxels))
    fill_core = np.zeros((n_voxels, n_voxels, n_voxels))

    # 对于每个体素，检查其8个角点的TPMS值来决定材料类型
    interface_thickness = thickness * 0.2  # 20% of wall thickness for interface

    for i in range(n_voxels):
        for j in range(n_voxels):
            for k in range(n_voxels):
                # 获取体素8个角点的TPMS值
                values = [
                    abs(tpms_values[i, j, k]),
                    abs(tpms_values[i + 1, j, k]),
                    abs(tpms_values[i + 1, j + 1, k]),
                    abs(tpms_values[i, j + 1, k]),
                    abs(tpms_values[i, j, k + 1]),
                    abs(tpms_values[i + 1, j, k + 1]),
                    abs(tpms_values[i + 1, j + 1, k + 1]),
                    abs(tpms_values[i, j + 1, k + 1])
                ]

                # 使用平均值或最大/最小值来判断
                avg_value = np.mean(values)

                if avg_value <= (thickness - interface_thickness):
                    tpms_core[i, j, k] = 1
                elif avg_value <= thickness:
                    interface[i, j, k] = 1
                else:
                    fill_core[i, j, k] = 1

    return tpms_core, interface, fill_core


def create_combined_part_from_voxels(model, tpms_core, interface, fill_core, mesh_size):
    """Create a single Abaqus part from all three voxel regions with shared nodes"""
    # Create new part
    part = model.Part(name='TPMS_Connected_Structure',
                      dimensionality=THREE_D,
                      type=DEFORMABLE_BODY)

    # Get dimensions of the voxel grid
    n_voxels_x, n_voxels_y, n_voxels_z = tpms_core.shape
    n_nodes_x = n_voxels_x + 1
    n_nodes_y = n_voxels_y + 1
    n_nodes_z = n_voxels_z + 1

    # 修正这里：使用传统字符串格式化
    print("Voxel grid: {} x {} x {}".format(n_voxels_x, n_voxels_y, n_voxels_z))
    print("Node grid: {} x {} x {}".format(n_nodes_x, n_nodes_y, n_nodes_z))

    # ====== STEP 1: Create all nodes first ======
    # Create a 3D array to store node labels for easy lookup
    node_labels = np.zeros((n_nodes_x, n_nodes_y, n_nodes_z), dtype=int)
    node_counter = 0

    print("Creating nodes...")
    for i in range(n_nodes_x):
        for j in range(n_nodes_y):
            for k in range(n_nodes_z):
                # Calculate node coordinates
                x = i * mesh_size
                y = j * mesh_size
                z = k * mesh_size

                # Create node
                node_counter += 1
                node = part.Node(coordinates=(x, y, z), label=node_counter)
                node_labels[i, j, k] = node_counter

    print("Created {} nodes total".format(node_counter))

    # ====== STEP 2: Create elements with shared nodes ======
    tpms_elements = []
    interface_elements = []
    fill_elements = []
    element_counter = 0

    print("Creating elements with shared nodes...")
    for i in range(n_voxels_x):
        for j in range(n_voxels_y):
            for k in range(n_voxels_z):
                # Check which region this voxel belongs to
                is_tpms = tpms_core[i, j, k] > 0
                is_interface = interface[i, j, k] > 0
                is_fill = fill_core[i, j, k] > 0

                if is_tpms or is_interface or is_fill:
                    # Get the 8 nodes for this voxel
                    # Note: HEX8 element node ordering in Abaqus:
                    # Bottom face: 1-2-3-4 (counter-clockwise)
                    # Top face: 5-6-7-8 (counter-clockwise)
                    node_indices = [
                        (i, j, k),  # Node 1
                        (i + 1, j, k),  # Node 2
                        (i + 1, j + 1, k),  # Node 3
                        (i, j + 1, k),  # Node 4
                        (i, j, k + 1),  # Node 5
                        (i + 1, j, k + 1),  # Node 6
                        (i + 1, j + 1, k + 1),  # Node 7
                        (i, j + 1, k + 1)  # Node 8
                    ]

                    # Get actual node objects using stored labels
                    node_objects = []
                    for idx in node_indices:
                        node_label = node_labels[idx]
                        node_objects.append(part.nodes.getFromLabel(node_label))

                    # Create element
                    element_counter += 1
                    element = part.Element(nodes=node_objects,
                                           elemShape=HEX8,
                                           label=element_counter)

                    # Assign to appropriate region
                    if is_tpms:
                        tpms_elements.append(element)
                    elif is_interface:
                        interface_elements.append(element)
                    elif is_fill:
                        fill_elements.append(element)

    print("Created {} elements total".format(element_counter))
    print("- TPMS elements: {}".format(len(tpms_elements)))
    print("- Interface elements: {}".format(len(interface_elements)))
    print("- Fill elements: {}".format(len(fill_elements)))

    # ====== STEP 3: Create element sets ======
    print("Creating element sets...")

    if tpms_elements:
        tpms_element_labels = [elem.label for elem in tpms_elements]
        tpms_elem_sequence = part.elements.sequenceFromLabels(tpms_element_labels)
        part.Set(elements=tpms_elem_sequence, name='TPMS_Set')
        print("Created TPMS_Set with {} elements".format(len(tpms_elements)))

    if interface_elements:
        interface_element_labels = [elem.label for elem in interface_elements]
        interface_elem_sequence = part.elements.sequenceFromLabels(interface_element_labels)
        part.Set(elements=interface_elem_sequence, name='Interface_Set')
        print("Created Interface_Set with {} elements".format(len(interface_elements)))

    if fill_elements:
        fill_element_labels = [elem.label for elem in fill_elements]
        fill_elem_sequence = part.elements.sequenceFromLabels(fill_element_labels)
        part.Set(elements=fill_elem_sequence, name='Fill_Set')
        print("Created Fill_Set with {} elements".format(len(fill_elements)))

    # ====== STEP 4: Verify connectivity ======
    print("\nVerifying connectivity...")

    # Check for shared nodes between elements
    element_connections = {}
    for elem in part.elements:
        for node in elem.getNodes():
            node_label = node.label
            if node_label not in element_connections:
                element_connections[node_label] = []
            element_connections[node_label].append(elem.label)

    # Count how many nodes are shared
    shared_nodes = 0
    max_shared = 0
    for node_label, elem_list in element_connections.items():
        if len(elem_list) > 1:
            shared_nodes += 1
            if len(elem_list) > max_shared:
                max_shared = len(elem_list)

    print("Shared nodes: {}/{}".format(shared_nodes, len(element_connections)))
    print("Maximum elements sharing a single node: {}".format(max_shared))

    if shared_nodes == 0:
        print("ERROR: No nodes are shared between elements! Model will have unconnected regions.")
    else:
        print("SUCCESS: Nodes are shared between elements. Model should be connected.")

    return part


def generate_modulus_values(num_layers):
    """Generate default modulus values for gradient layers"""
    # Default modulus range: 100000 to 200000 MPa
    min_modulus = 100000.0
    max_modulus = 200000.0
    # Linear gradient from min to max
    return [min_modulus + (max_modulus - min_modulus) * i / (num_layers - 1) for i in range(num_layers)]


def assign_gradient_materials_to_combined_part(model, part, gradient_direction, num_layers):
    """Assign gradient materials to TPMS elements in the combined part"""
    print('\nAssigning gradient materials to TPMS elements in combined part...')

    # Check if TPMS_Set exists
    if 'TPMS_Set' not in part.sets.keys():
        print('Warning: No TPMS_Set found in part. Gradient materials not applied.')
        return

    # Get TPMS elements
    tpms_elements = part.sets['TPMS_Set'].elements
    if len(tpms_elements) == 0:
        print('Warning: TPMS_Set is empty. Gradient materials not applied.')
        return

    # Generate default modulus values
    modulus_values = generate_modulus_values(num_layers)

    print('Default modulus values: {}'.format([format(v, '.1f') for v in modulus_values]))
    print('Note: You can modify these values manually in Abaqus after generation')

    # Get coordinate range for gradient direction from TPMS elements only
    if gradient_direction == 'X':
        coord_index = 0
    elif gradient_direction == 'Y':
        coord_index = 1
    else:  # Z direction
        coord_index = 2

    # Get coordinates from TPMS elements
    coords = []
    for elem in tpms_elements:
        nodes = elem.getNodes()
        elem_coord = sum([n.coordinates[coord_index] for n in nodes]) / len(nodes)
        coords.append(elem_coord)

    coord_min = min(coords)
    coord_max = max(coords)
    height = coord_max - coord_min
    layer_height = height / float(num_layers)

    print('{} range for TPMS elements: {:.3f} to {:.3f}'.format(gradient_direction, coord_min, coord_max))
    print('Height: {:.3f}'.format(height))
    print('Layer height: {:.3f}'.format(layer_height))

    # Process elements layer by layer
    print('\nProcessing gradient layers...')
    for i in range(num_layers):
        # Calculate layer bounds
        layer_min = coord_min + i * layer_height
        layer_max = layer_min + layer_height

        # Get modulus value for this layer
        E_current = modulus_values[i]

        # Create material
        mat_name = 'TPMS_Gradient_Material_{}'.format(i)
        material = model.Material(name=mat_name)
        material.Elastic(table=((E_current, POISSON),))

        # Create section
        sect_name = 'TPMS_Gradient_Section_{}'.format(i)
        model.HomogeneousSolidSection(
            name=sect_name,
            material=mat_name,
            thickness=None
        )

        # Select elements in this layer using element labels
        layer_element_labels = []
        for elem in tpms_elements:
            # Get element center coordinate in gradient direction
            nodes = elem.getNodes()
            elem_coord = sum([n.coordinates[coord_index] for n in nodes]) / len(nodes)
            # For the last layer, include elements at the maximum boundary
            if i == num_layers - 1:  # Last layer
                if layer_min <= elem_coord <= layer_max:
                    layer_element_labels.append(elem.label)
            else:  # Other layers
                if layer_min <= elem_coord < layer_max:
                    layer_element_labels.append(elem.label)

        if layer_element_labels:
            # Use element labels to create region
            try:
                # Create element sequence from labels
                elem_sequence = part.elements.sequenceFromLabels(layer_element_labels)
                region = regionToolset.Region(elements=elem_sequence)
                part.SectionAssignment(
                    region=region,
                    sectionName=sect_name,
                    offset=0.0,
                    offsetType=MIDDLE_SURFACE,
                    offsetField=''
                )
                print(
                    'Successfully assigned section for layer {} with {} elements'.format(i, len(layer_element_labels)))
            except Exception as e:
                print('Warning: Could not assign section for layer {}. Error: {}'.format(i, str(e)))
                # Try alternative method
                try:
                    # Create set first, then assign
                    set_name = 'Temp_Layer_{}'.format(i)
                    part.Set(elements=part.elements.sequenceFromLabels(layer_element_labels), name=set_name)
                    region = regionToolset.Region(elements=part.sets[set_name].elements)
                    part.SectionAssignment(
                        region=region,
                        sectionName=sect_name,
                        offset=0.0,
                        offsetType=MIDDLE_SURFACE,
                        offsetField=''
                    )
                    print('Successfully assigned section for layer {} using alternative method'.format(i))
                except Exception as e2:
                    print('Alternative method also failed for layer {}. Error: {}'.format(i, str(e2)))

            print('Layer {}: {} elements, E = {:.1f} MPa'.format(
                i + 1, len(layer_element_labels), E_current))

    # Final validation: check that all TPMS elements have been assigned to gradient layers
    print('\n=== Gradient Assignment Validation ===')
    total_tpms_elements = len(tpms_elements)
    total_assigned_to_gradient = 0

    for i in range(num_layers):
        layer_min = coord_min + i * layer_height
        layer_max = layer_min + layer_height
        layer_count = 0

        for elem in tpms_elements:
            nodes = elem.getNodes()
            elem_coord = sum([n.coordinates[coord_index] for n in nodes]) / len(nodes)
            if i == num_layers - 1:  # Last layer
                if layer_min <= elem_coord <= layer_max:
                    layer_count += 1
            else:  # Other layers
                if layer_min <= elem_coord < layer_max:
                    layer_count += 1

        total_assigned_to_gradient += layer_count
        print('Layer {}: {} elements assigned'.format(i + 1, layer_count))

    print('Total TPMS elements: {}'.format(total_tpms_elements))
    print('Total assigned to gradient: {}'.format(total_assigned_to_gradient))

    if total_assigned_to_gradient != total_tpms_elements:
        print('ERROR: {} TPMS elements not assigned to any gradient layer!'.format(
            total_tpms_elements - total_assigned_to_gradient))
        print('This will cause missing property definitions in Abaqus!')
    else:
        print('✓ All TPMS elements successfully assigned to gradient layers')


def create_model_and_assembly(tpms_core, interface, fill_core, mesh_size, enable_gradient, gradient_direction,
                              num_layers):
    """Create the Abaqus model and assembly using a single combined part"""
    # Create a new model
    model = mdb.Model(name='TPMS_Model')

    # Create a single combined part from all three regions
    print('\nCreating combined part with all regions...')
    combined_part = create_combined_part_from_voxels(model, tpms_core, interface, fill_core, mesh_size)

    # Create materials
    print('\nCreating materials...')

    # Create TPMS material
    if enable_gradient:
        print('Applying gradient materials to TPMS region...')
        assign_gradient_materials_to_combined_part(model, combined_part, gradient_direction, num_layers)
    else:
        # Create homogeneous TPMS material
        tpms_material = model.Material(name='TPMS_Material')
        tpms_material.Elastic(table=((210000.0, 0.3),))
        print('Created TPMS material: E = 210000 MPa')

    # Create interface material
    interface_material = model.Material(name='Interface_Material')
    interface_material.Elastic(table=((105000.0, 0.3),))
    print('Created Interface material: E = 105000 MPa')

    # Create fill material
    fill_material = model.Material(name='Fill_Material')
    fill_material.Elastic(table=((21000.0, 0.3),))
    print('Created Fill material: E = 21000 MPa')

    # Create sections
    print('\nCreating sections...')

    if not enable_gradient:
        model.HomogeneousSolidSection(name='TPMS_Section', material='TPMS_Material')
    model.HomogeneousSolidSection(name='Interface_Section', material='Interface_Material')
    model.HomogeneousSolidSection(name='Fill_Section', material='Fill_Material')

    # Assign sections to different element sets
    print('\nAssigning sections to element sets...')

    try:
        if 'TPMS_Set' in combined_part.sets.keys() and not enable_gradient:
            tpms_set_elements = combined_part.sets['TPMS_Set'].elements
            if len(tpms_set_elements) > 0:
                region = regionToolset.Region(elements=tpms_set_elements)
                combined_part.SectionAssignment(region=region, sectionName='TPMS_Section')
                print('Assigned TPMS_Section to {} TPMS elements'.format(len(tpms_set_elements)))
            else:
                print('Warning: TPMS_Set exists but contains no elements')
    except Exception as e:
        print('ERROR: Could not assign TPMS section: {}'.format(str(e)))
        print('This will cause missing property definitions in Abaqus!')

    try:
        if 'Interface_Set' in combined_part.sets.keys():
            interface_set_elements = combined_part.sets['Interface_Set'].elements
            if len(interface_set_elements) > 0:
                region = regionToolset.Region(elements=interface_set_elements)
                combined_part.SectionAssignment(region=region, sectionName='Interface_Section')
                print('Assigned Interface_Section to {} Interface elements'.format(len(interface_set_elements)))
            else:
                print('Warning: Interface_Set exists but contains no elements')
    except Exception as e:
        print('ERROR: Could not assign Interface section: {}'.format(str(e)))
        print('This will cause missing property definitions in Abaqus!')

    try:
        if 'Fill_Set' in combined_part.sets.keys():
            fill_set_elements = combined_part.sets['Fill_Set'].elements
            if len(fill_set_elements) > 0:
                region = regionToolset.Region(elements=fill_set_elements)
                combined_part.SectionAssignment(region=region, sectionName='Fill_Section')
                print('Assigned Fill_Section to {} Fill elements'.format(len(fill_set_elements)))
            else:
                print('Warning: Fill_Set exists but contains no elements')
    except Exception as e:
        print('ERROR: Could not assign Fill section: {}'.format(str(e)))
        print('This will cause missing property definitions in Abaqus!')

    # Create assembly with single instance
    print('\nCreating assembly...')
    assembly = model.rootAssembly
    combined_instance = assembly.Instance(name='TPMS_Combined_Instance', part=combined_part, dependent=ON)

    print('Successfully created combined TPMS model with single part and instance')

    # Final validation: Check that all elements have section assignments
    print('\n=== FINAL COMPREHENSIVE VALIDATION ===')
    total_part_elements = len(combined_part.elements)

    # Count elements in each set
    tpms_set_count = 0
    interface_set_count = 0
    fill_set_count = 0

    if 'TPMS_Set' in combined_part.sets.keys():
        tpms_set_count = len(combined_part.sets['TPMS_Set'].elements)
        print('TPMS_Set contains: {} elements'.format(tpms_set_count))
    else:
        print('WARNING: TPMS_Set not found!')

    if 'Interface_Set' in combined_part.sets.keys():
        interface_set_count = len(combined_part.sets['Interface_Set'].elements)
        print('Interface_Set contains: {} elements'.format(interface_set_count))
    else:
        print('WARNING: Interface_Set not found!')

    if 'Fill_Set' in combined_part.sets.keys():
        fill_set_count = len(combined_part.sets['Fill_Set'].elements)
        print('Fill_Set contains: {} elements'.format(fill_set_count))
    else:
        print('WARNING: Fill_Set not found!')

    assigned_elements = tpms_set_count + interface_set_count + fill_set_count

    print('\n--- Element Assignment Summary ---')
    print('Total elements in part: {}'.format(total_part_elements))
    print('Elements with set assignments: {}'.format(assigned_elements))
    print('- TPMS elements: {}'.format(tpms_set_count))
    print('- Interface elements: {}'.format(interface_set_count))
    print('- Fill elements: {}'.format(fill_set_count))

    if assigned_elements != total_part_elements:
        missing_elements = total_part_elements - assigned_elements
        print('\n🚨 CRITICAL ERROR: {} elements are missing set assignments!'.format(missing_elements))
        print('This WILL cause the "missing property definitions" error in Abaqus!')
        print('Recommendation: Check voxel generation logic and element assignment conditions.')
    else:
        print('\n✓ All elements have proper set assignments')

    # Additional check: verify section assignments exist
    print('\n--- Section Assignment Verification ---')
    sections_assigned = len(combined_part.sectionAssignments)
    expected_sections = 0
    if enable_gradient and tpms_set_count > 0:
        expected_sections = num_layers + (1 if interface_set_count > 0 else 0) + (1 if fill_set_count > 0 else 0)
    else:
        expected_sections = (1 if tpms_set_count > 0 else 0) + (1 if interface_set_count > 0 else 0) + (
            1 if fill_set_count > 0 else 0)

    print('Section assignments created: {}'.format(sections_assigned))
    print('Expected section assignments: {}'.format(expected_sections))

    if sections_assigned < expected_sections:
        print('WARNING: Fewer section assignments than expected!')
    else:
        print('Section assignments appear complete')

    # Detailed element-by-element verification if there are missing assignments
    if assigned_elements != total_part_elements:
        print('\n=== DETAILED DIAGNOSTIC ===')
        print('Performing element-by-element check...')

        # Get all element labels that should be assigned
        all_assigned_labels = set()

        if 'TPMS_Set' in combined_part.sets.keys():
            tpms_labels = [elem.label for elem in combined_part.sets['TPMS_Set'].elements]
            all_assigned_labels.update(tpms_labels)
            print('TPMS_Set element labels: {}'.format(len(tpms_labels)))

        if 'Interface_Set' in combined_part.sets.keys():
            interface_labels = [elem.label for elem in combined_part.sets['Interface_Set'].elements]
            all_assigned_labels.update(interface_labels)
            print('Interface_Set element labels: {}'.format(len(interface_labels)))

        if 'Fill_Set' in combined_part.sets.keys():
            fill_labels = [elem.label for elem in combined_part.sets['Fill_Set'].elements]
            all_assigned_labels.update(fill_labels)
            print('Fill_Set element labels: {}'.format(len(fill_labels)))

        # Get all part element labels
        all_part_labels = set([elem.label for elem in combined_part.elements])

        # Find missing labels
        missing_labels = all_part_labels - all_assigned_labels

        print('Total part element labels: {}'.format(len(all_part_labels)))
        print('Total assigned element labels: {}'.format(len(all_assigned_labels)))
        print('Missing element labels: {}'.format(len(missing_labels)))

        if missing_labels:
            print('Missing element label list: {}'.format(sorted(list(missing_labels))[:10]))  # Show first 10
            if len(missing_labels) > 10:
                print('... and {} more missing elements'.format(len(missing_labels) - 10))


def verify_model_connectivity(part):
    """Verify that all elements are connected in the model"""
    print("\n=== Connectivity Verification ===")

    # Get all elements
    all_elements = list(part.elements)
    n_elements = len(all_elements)

    if n_elements == 0:
        print("No elements in the model")
        return False

    # Build adjacency graph
    element_to_nodes = {}
    node_to_elements = {}

    for elem in all_elements:
        nodes = [n.label for n in elem.getNodes()]
        element_to_nodes[elem.label] = nodes
        for node_label in nodes:
            if node_label not in node_to_elements:
                node_to_elements[node_label] = []
            node_to_elements[node_label].append(elem.label)

    # Find connected components
    visited = set()
    components = []

    for elem_label in element_to_nodes.keys():
        if elem_label not in visited:
            # Start a new component
            component = set()
            stack = [elem_label]

            while stack:
                current = stack.pop()
                if current not in visited:
                    visited.add(current)
                    component.add(current)

                    # Get neighboring elements through shared nodes
                    for node_label in element_to_nodes[current]:
                        for neighbor_elem in node_to_elements.get(node_label, []):
                            if neighbor_elem not in visited:
                                stack.append(neighbor_elem)

            components.append(component)

    print("Found {} connected component(s)".format(len(components)))

    if len(components) > 1:
        print("WARNING: Model is not fully connected! Found {} separate regions:".format(len(components)))
        for i, comp in enumerate(components):
            print("  Component {}: {} elements".format(i + 1, len(comp)))

        # Check if this might be expected (e.g., isolated fill regions)
        # You might want to check if each component belongs to a different material
    else:
        print("SUCCESS: All elements are connected in a single component")

    return len(components) == 1


def main():
    """Main function to create TPMS model"""
    while True:
        # Get user input
        user_input_data = get_user_input()
        if user_input_data is None:
            retry = getWarningReply(
                message='Invalid input parameters. Would you like to try again?',
                buttons=(YES, NO)
            )
            if retry == YES:
                continue
            else:
                return
        break

    tpms_func = user_input_data['tpms_func']
    cell_size = user_input_data['cell_size']
    thickness = user_input_data['thickness']
    mesh_size = user_input_data['mesh_size']
    tpms_type = user_input_data['tpms_type']
    enable_gradient = user_input_data['enable_gradient']
    gradient_direction = user_input_data['gradient_direction']
    num_layers = user_input_data['num_layers']

    # Generate voxel model for different regions
    print('\nGenerating TPMS structure...')
    tpms_core, interface, fill_core = generate_voxel_model(tpms_func, cell_size, thickness, mesh_size)

    # Validate the generated model
    print('\nValidating generated model...')
    tpms_elements = np.sum(tpms_core)
    interface_elements = np.sum(interface)
    fill_elements = np.sum(fill_core)
    total_elements = tpms_elements + interface_elements + fill_elements

    print('Model statistics:')
    print('- TPMS elements: {}'.format(int(tpms_elements)))
    print('- Interface elements: {}'.format(int(interface_elements)))
    print('- Fill elements: {}'.format(int(fill_elements)))
    print('- Total elements: {}'.format(int(total_elements)))

    # Check for potential issues
    if total_elements == 0:
        print('ERROR: No elements generated. Check your parameters.')
        return

    if tpms_elements == 0:
        print('WARNING: No TPMS elements generated. The thickness might be too small.')

    if fill_elements == 0:
        print('WARNING: No fill elements generated. This might cause instability.')
        retry = getWarningReply(
            message='No fill material detected. This may cause numerical issues. Continue anyway?',
            buttons=(YES, NO)
        )
        if retry == NO:
            return

    # Check element size quality
    # Check element size quality
    n_voxels = int(cell_size / mesh_size)  # 体素数量
    n_points = n_voxels + 1  # 节点数量
    if n_voxels < 5:
        print('WARNING: Very coarse mesh detected ({} voxels per cell). Consider reducing mesh size for better quality.'.format(n_voxels))

    # Create model and assembly
    print('\nCreating Abaqus model...')
    create_model_and_assembly(tpms_core, interface, fill_core, mesh_size, enable_gradient, gradient_direction,
                              num_layers)

    # ====== 添加连接性验证 ======
    print('\n=== Running connectivity check ===')
    try:
        # 获取模型中的部件
        model = mdb.models['TPMS_Model']
        part = model.parts['TPMS_Connected_Structure']  # 注意使用新的部件名称

        # 验证连接性
        is_connected = verify_model_connectivity(part)

        if not is_connected:
            print('\nWARNING: Model has unconnected regions. This may cause convergence issues.')
            print('Consider:')
            print('1. Reducing wall thickness')
            print('2. Increasing mesh density')
            print('3. Using a simpler TPMS type')
        else:
            print('\nSUCCESS: Model is fully connected.')
    except KeyError as e:
        print('\nERROR: Could not find model or part. {}'.format(e))
        print('Model name might have changed. Available models:', list(mdb.models.keys()))
        if 'TPMS_Model' in mdb.models:
            model = mdb.models['TPMS_Model']
            print('Available parts in this model:', list(model.parts.keys()))
    except Exception as e:
        print('\nERROR during connectivity check: {}'.format(e))
    # ====== 连接性验证结束 ======

    print('\nTPMS model generation completed:')
    print('- Type: {}'.format(tpms_type))
    print('- Cell size: {} mm'.format(cell_size))
    print('- Wall thickness: {}'.format(thickness))
    print('- Mesh size: {}'.format(mesh_size))
    print('- Elements per cell: {}'.format(n_points))
    if enable_gradient:
        print('- Gradient material: ENABLED')
        print('- Gradient direction: {}'.format(gradient_direction))
        print('- Number of gradient layers: {}'.format(num_layers))
        print('- Default modulus range: 100000 to 200000 MPa (can be modified in Abaqus)')
    else:
        print('- Gradient material: DISABLED')

    print('\n=== IMPORTANT NOTES FOR ANALYSIS ===')
    print('1. The model now uses a SINGLE COMBINED PART with different element sets')
    print('2. This eliminates connectivity issues between separate parts')
    print('3. Material properties are assigned to different element sets:')
    print('   - TPMS_Set: {} MPa'.format(210000.0 if not enable_gradient else 'Gradient 100k-200k'))
    print('   - Interface_Set: 105000 MPa')
    print('   - Fill_Set: 21000 MPa')
    print('4. Assembly contains a single instance: TPMS_Combined_Instance')
    print('5. If you still encounter numerical issues:')
    print('   - Check your boundary conditions (ensure proper constraints)')
    print('   - Use smaller load increments in your analysis step')
    print('   - Consider using implicit analysis with automatic stabilization')
    print('   - Apply loads and constraints to the combined instance')
    print('6. For better convergence:')
    print('   - Use displacement-controlled loading instead of force-controlled')
    print('   - Enable geometric nonlinearity if large deformations are expected')
    print('   - The unified mesh should eliminate most numerical singularities')


# Run main function
main()