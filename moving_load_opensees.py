"""
===========================================================================
  Title: Simple 2D Beam Bridge Model in OpenSees with Moving Loads
  Author:  Mohammad Talebi-Kalaleh
           Research Assistant at University of Alberta, Structural Engineering
  Email:   talebika@ualberta.ca

  Description:
    This is a comprehensive example of how to model a 2D beam bridge
    subjected to moving loads (e.g. truck axles) using OpenSees in Python.
    Here, we use NumPy's built-in 1D interpolation (np.interp) to 
    interpolate the trailing axle position over time rather than writing 
    a custom linear interpolation routine.

  Usage:
    - The code is organized into functions for clarity:
        1) build_bridge_model(...)
        2) interpolate_position_at_time(...)
        3) get_nodal_forces_for_all_times(...)
        4) run_dynamic_analysis_with_path_ts(...)
        5) post_process_and_plot(...)
        6) main(...)
    - Provide your own input parameters for bridge geometry, material
      properties, damping ratios, moving loads, etc. in main() or from
      external sources.

  Citation:
    If you find this code helpful for your research or professional work,
    please cite:
      "M. Talebi Kalaleh, Research Assistant, University of Alberta"
    and/or reference this source code. For inquiries, feel free to contact:
      talebika@ualberta.ca

  Prerequisites:
    - numpy
    - matplotlib
    - OpenSeesPy (Python version of OpenSees)
      (See https://openseespydoc.readthedocs.io for installation details)

===========================================================================
"""

import openseespy.opensees as ops
import openseespy.postprocessing.Get_Rendering as opsv
import numpy as np
import matplotlib.pyplot as plt


def build_bridge_model(bridge_length,
                       support_positions,
                       E,
                       rho,
                       A,
                       I,
                       damping_ratio,
                       mesh_size=0.5,
                       g=9.81):
    """
    Build a 2D beam (3-dof at each node) bridge model in OpenSees using
    elasticBeamColumn elements.

    Parameters
    ----------
    bridge_length : float
        Total length of the bridge (meters).
    support_positions : list of float
        Positions of the supports along the bridge measured from the left end.
    E : float
        Elastic modulus (Pa or N/m^2).
    rho : float
        Mass density (kg/m^3) -- or linear mass density if desired.
    A : list of float
        Cross-sectional area(s). For prismatic spans, supply one A per span.
    I : list of float
        Second moment(s) of area. Same indexing as A.
    damping_ratio : float
        Desired damping ratio for the first two modes.
    mesh_size : float, optional
        Finite element mesh size (m) used to divide each span. Default 0.5.
    g : float, optional
        Acceleration due to gravity (m/s^2). Default is 9.81.

    Returns
    -------
    bridge_dict : dict
        Dictionary storing metadata about the constructed model, including:
        - node tags, element tags, node coordinates, mode shapes, etc.
    """

    # Wipe any existing OpenSees model
    ops.wipe()

    # Define a 2D model, 3 dofs/node
    ops.model('BasicBuilder', '-ndm', 2, '-ndf', 3)

    # Ensure 0.0 and bridge_length are included in supports
    all_support_positions = sorted(set([0.0] + support_positions + [bridge_length]))

    current_node_tag = 0
    span_node_tags = []
    all_node_positions = []

    # Create nodes across each span
    for i in range(len(all_support_positions) - 1):
        span_start = all_support_positions[i]
        span_end   = all_support_positions[i + 1]
        span_length = span_end - span_start

        n_seg = int(np.ceil(span_length / mesh_size))
        dx = span_length / n_seg

        this_span_nodes = []
        for j in range(n_seg + 1):
            xcoord = span_start + j * dx
            if i == 0:
                current_node_tag += 1
                ops.node(current_node_tag, xcoord, 0.0)
                all_node_positions.append(xcoord)
            else:
                if j != 0:
                    current_node_tag += 1
                    ops.node(current_node_tag, xcoord, 0.0)
                    all_node_positions.append(xcoord)
            this_span_nodes.append(current_node_tag)

        span_node_tags.append(this_span_nodes)

    # Flatten node tags, ensure uniqueness
    all_node_tags = [tag for span_nodes in span_node_tags for tag in span_nodes]
    all_node_tags = list(set(all_node_tags))
    number_of_nodes = len(all_node_tags)

    # Assign boundary conditions (hinge for first, roller for others)
    support_node_tags = []
    for sPosIdx, sPos in enumerate(all_support_positions):
        idx = np.argmin(np.abs(np.array(all_node_positions) - sPos))
        support_node_tag = all_node_tags[idx]
        support_node_tags.append(support_node_tag)
        if sPosIdx == 0:
            # Hinge
            ops.fix(support_node_tag, 1, 1, 0)
        else:
            # Roller
            ops.fix(support_node_tag, 0, 1, 0)

    # Geometric transformation
    transf_tag = 1
    ops.geomTransf('Linear', transf_tag)

    # If A or I are single values, replicate them across all spans
    if len(A) == 1:
        A = A * (len(all_support_positions) - 1)
    if len(I) == 1:
        I = I * (len(all_support_positions) - 1)

    current_element_tag = 1
    span_element_tags = []

    for i, span_nodes in enumerate(span_node_tags):
        Ei = E
        Ai = A[i]
        Ii = I[i]

        this_span_elems = []
        for j in range(len(span_nodes) - 1):
            nd_i = span_nodes[j]
            nd_j = span_nodes[j + 1]
            ops.element('elasticBeamColumn', current_element_tag,
                        nd_i, nd_j,
                        Ai, Ei, Ii, transf_tag)
            this_span_elems.append(current_element_tag)
            current_element_tag += 1

        span_element_tags.append(this_span_elems)

    # Assign mass to nodes
    for i, span_nodes in enumerate(span_node_tags):
        linear_mass_i = rho * A[i]  # linear mass
        for j in range(len(span_nodes) - 1):
            nd_i = span_nodes[j]
            nd_j = span_nodes[j + 1]
            xi = ops.nodeCoord(nd_i)[0]
            xj = ops.nodeCoord(nd_j)[0]
            le = xj - xi
            ops.mass(nd_i, 0.0,
                     ops.nodeMass(nd_i)[1] + 0.5 * linear_mass_i * le,
                     0.0)
            ops.mass(nd_j, 0.0,
                     ops.nodeMass(nd_j)[1] + 0.5 * linear_mass_i * le,
                     0.0)

    # Eigen analysis (first 2 modes)
    num_eigen = 2
    eigenvals = ops.eigen(num_eigen)
    omega1 = np.sqrt(eigenvals[0])
    omega2 = np.sqrt(eigenvals[1])

    # Plot mode shapes (optional)
    opsv.plot_mode_shape(1)
    plt.title(f'Mode 1, f_1={round(omega1/(2*np.pi), 2)} Hz')
    plt.show()

    opsv.plot_mode_shape(2)
    plt.title(f'Mode 2, f_2={round(omega2/(2*np.pi), 2)} Hz')
    plt.show()

    # Return a dictionary with metadata
    bridge_dict = {
        'all_node_tags': all_node_tags,
        'all_node_positions': all_node_positions,
        'span_node_tags': span_node_tags,
        'span_element_tags': span_element_tags,
        'all_support_positions': all_support_positions,
        'support_node_tags': support_node_tags,
        'number_of_nodes': number_of_nodes,
        'num_eigen': num_eigen,
        'omega': [omega1, omega2]
    }
    return bridge_dict


def interpolate_position_at_time(t, time_position_table):
    """
    Interpolate the trailing axle position at time 't' using np.interp.
    If t < t_min, returns x_min; if t > t_max, returns x_max.
    """
    # Extract separate arrays for times and positions
    times = np.array([row[0] for row in time_position_table])
    positions = np.array([row[1] for row in time_position_table])

    # Use np.interp for linear interpolation
    x_t = np.interp(t, times, positions)
    return x_t


def get_nodal_forces_for_all_times(time_steps,
                                   trailing_axle_path,
                                   vehicle_config,
                                   node_tags,
                                   node_positions):
    """
    Pre-compute the full time history of vertical load at each node due to
    a vehicle whose trailing axle position is given by 'trailing_axle_path',
    and whose axle loads/distances are given in 'vehicle_config'.

    vehicle_config = {
      'loads': [axle1_load, axle2_load, ...],
      'distances': [0.0, dist2, dist3, ...]
    }

    Parameters
    ----------
    time_steps : np.array
        Array of times (s).
    trailing_axle_path : list of [time, position]
        The (time, x-position) history of the trailing axle.
    vehicle_config : dict
        'loads' -> list of axle load magnitudes (N),
        'distances' -> list of distances from trailing axle (m).
    node_tags : list of int
        Node tags in order along the bridge.
    node_positions : list of float
        The x-coordinates of the node_tags (same order).

    Returns
    -------
    node_load_history : dict
        key=node_tag, value = array of shape (len(time_steps),)
        containing the downward load for that node.
    """

    axle_loads = vehicle_config['loads']
    axle_distances = vehicle_config['distances']
    n_axles = len(axle_loads)  # number of axles

    node_load_history = {nd: np.zeros(len(time_steps)) for nd in node_tags}

    for it, t in enumerate(time_steps):
        # Trailing axle position by interpolation
        trailing_pos = interpolate_position_at_time(t, trailing_axle_path)

        # For each axle, compute x-position
        for ax_idx in range(n_axles):
            x_axle = trailing_pos + axle_distances[ax_idx]
            load_val = axle_loads[ax_idx]

            # Skip if off the bridge
            if x_axle < min(node_positions) or x_axle > max(node_positions):
                continue

            # Find segment to do linear interpolation
            for iSeg in range(len(node_positions) - 1):
                x_i = node_positions[iSeg]
                x_j = node_positions[iSeg + 1]
                nd_i = node_tags[iSeg]
                nd_j = node_tags[iSeg + 1]

                if x_i <= x_axle <= x_j:
                    seg_length = x_j - x_i
                    if seg_length < 1e-12:
                        node_load_history[nd_i][it] += load_val
                    else:
                        r = (x_axle - x_i)/seg_length
                        node_load_history[nd_i][it] += load_val * (1.0 - r)
                        node_load_history[nd_j][it] += load_val * r
                    break

    return node_load_history


def run_dynamic_analysis_with_path_ts(bridge_dict,
                                      trailing_axle_path,
                                      vehicle_config,
                                      dt):
    """
    Constructs Path TimeSeries for each node, capturing the node's 
    time-varying load. Then performs a transient analysis step-by-step
    to collect nodal/element responses at each time step.

    Parameters
    ----------
    bridge_dict : dict
        Output of build_bridge_model(...).
    trailing_axle_path : list of [time, position]
        Time-position history for the trailing axle.
    vehicle_config : dict
        'loads': [...], 'distances': [...]
    dt : float
        Time step (s).

    Returns
    -------
    results : dict
        Contains time array, midspan disp/vel/acc, shear/moment over time,
        and envelope values.
    """

    node_tags          = bridge_dict['all_node_tags']
    node_positions     = bridge_dict['all_node_positions']
    span_element_tags  = bridge_dict['span_element_tags']
    all_support_pos    = bridge_dict['all_support_positions']
    omega1, omega2     = bridge_dict['omega']

    # Determine start/end times
    sorted_path = sorted(trailing_axle_path, key=lambda row: row[0])
    t_start = sorted_path[0][0]
    t_end   = sorted_path[-1][0]
    total_time = t_end - t_start if t_end > t_start else 1e-3

    nSteps = int(np.ceil(total_time / dt)) + 1
    time_steps = np.linspace(t_start, t_end, nSteps)

    # Build nodal load history
    node_load_history = get_nodal_forces_for_all_times(
        time_steps, trailing_axle_path, vehicle_config,
        node_tags, node_positions
    )

    # Create Path TimeSeries & Patterns
    for nd in node_tags:
        ts_tag = 1000 + nd
        ops.timeSeries('Path', ts_tag,
                       '-time',  *time_steps,
                       '-values', *node_load_history[nd],
                       '-factor', 1.0,
                       '-useLast',
                       '-prependZero')
        pat_tag = 2000 + nd
        ops.pattern('Plain', pat_tag, ts_tag)
        # Vertical downward => (0, -1, 0)
        ops.load(nd, 0.0, -1.0, 0.0)

    # Flatten elements
    all_element_tags = [elem for sp_elems in span_element_tags for elem in sp_elems]

    # Analysis configuration
    ops.constraints('Transformation')
    ops.numberer('Plain')
    ops.system('BandGeneral')
    ops.algorithm('Linear')
    ops.integrator('Newmark', 0.5, 0.25)
    ops.analysis('Transient')
    ops.setTime(t_start)

    # Identify midspan node (or other node) for response tracking
    x_left  = min(node_positions)
    x_right = max(node_positions)
    x_mid   = 0.5*(x_left + x_right)
    idx_mid = np.argmin(np.abs(np.array(node_positions) - x_mid))
    mid_node_tag = node_tags[idx_mid]

    # Lists to store time, disp, vel, acc
    time_list = []
    disp_list = []
    vel_list  = []
    acc_list  = []

    # Track shear & moment time histories
    shear_records  = []
    moment_records = []

    # Step-by-step analysis
    for iStep in range(nSteps):
        current_time = t_start + iStep*dt
        ops.setTime(current_time)

        ok = ops.analyze(1, dt)
        if ok != 0:
            print(f"Analysis failed at step={iStep}, time={ops.getTime():.3f}.")
            break

        # Record time & midspan
        tNow = ops.getTime()
        time_list.append(tNow)
        disp_list.append(ops.nodeDisp(mid_node_tag, 2))
        vel_list.append(ops.nodeVel(mid_node_tag, 2))
        acc_list.append(ops.nodeAccel(mid_node_tag, 2))

        # Shear & moment in each element (i-end)
        # eleForce => [Fx_i, Fz_i, My_i, Fx_j, Fz_j, My_j]
        shear_step  = []
        moment_step = []
        for eleTag in all_element_tags:
            forces = ops.eleForce(eleTag)
            shear_i  = forces[1]
            moment_i = forces[2]
            shear_step.append(shear_i)
            moment_step.append(moment_i)
        shear_records.append(shear_step)
        moment_records.append(moment_step)

    # Convert to arrays
    time_array = np.array(time_list)
    disp_array = np.array(disp_list)
    vel_array  = np.array(vel_list)
    acc_array  = np.array(acc_list)

    shear_time_history  = np.array(shear_records)
    moment_time_history = np.array(moment_records)

    # Envelopes
    shear_envelope_max  = np.max(shear_time_history, axis=0)
    shear_envelope_min  = np.min(shear_time_history, axis=0)
    moment_envelope_max = np.max(moment_time_history, axis=0)
    moment_envelope_min = np.min(moment_time_history, axis=0)

    results = {
        'time': time_array,
        'midspan_disp': disp_array,
        'midspan_vel': vel_array,
        'midspan_acc': acc_array,
        'shear_time_history': shear_time_history,
        'moment_time_history': moment_time_history,
        'shear_envelope_max': shear_envelope_max,
        'shear_envelope_min': shear_envelope_min,
        'moment_envelope_max': moment_envelope_max,
        'moment_envelope_min': moment_envelope_min,
        'element_tags': all_element_tags
    }
    return results


def post_process_and_plot(results_dict):
    """
    Simple post-processing:
      - Plots midspan displacement, velocity, acceleration vs. time.
      - Plots shear & moment envelope diagrams.
    """
    t    = results_dict['time']
    disp = results_dict['midspan_disp']
    vel  = results_dict['midspan_vel']
    acc  = results_dict['midspan_acc']

    fig, axs = plt.subplots(3, 1, figsize=(9, 10), sharex=True)

    axs[0].plot(t, 1000*disp, 'b-', label='Midspan Disp (mm)')
    axs[0].legend()
    axs[0].grid(True)

    axs[1].plot(t, 1000*vel, 'r-', label='Midspan Vel (mm/s)')
    axs[1].legend()
    axs[1].grid(True)

    axs[2].plot(t, acc / 9.81, 'g-', label='Midspan Acc (g)')
    axs[2].legend()
    axs[2].grid(True)
    axs[2].set_xlabel('Time (s)')
    plt.suptitle('Midspan Response Time Histories')
    plt.tight_layout()
    plt.show()

    # Envelope
    shear_max  = results_dict['shear_envelope_max']
    shear_min  = results_dict['shear_envelope_min']
    moment_max = results_dict['moment_envelope_max']
    moment_min = results_dict['moment_envelope_min']
    elem_tags  = results_dict['element_tags']

    x_elem = []
    for et in elem_tags:
        nd_i, nd_j = ops.eleNodes(et)
        xi = ops.nodeCoord(nd_i)[0]
        x_elem.append(xi)

    x_elem = np.array(x_elem)
    sort_idx = np.argsort(x_elem)
    x_sorted = x_elem[sort_idx]

    sh_max_s = shear_max[sort_idx] / 1000.0
    sh_min_s = shear_min[sort_idx] / 1000.0
    m_max_s  = moment_max[sort_idx] / 1000.0
    m_min_s  = moment_min[sort_idx] / 1000.0

    fig2, ax2 = plt.subplots(2, 1, figsize=(9, 8), sharex=True)

    ax2[0].plot(x_sorted, sh_max_s, 'r-', label='V max')
    ax2[0].plot(x_sorted, sh_min_s, 'b-', label='V min')
    ax2[0].set_ylabel('Shear (kN)')
    ax2[0].legend()
    ax2[0].grid(True)

    ax2[1].plot(x_sorted, m_max_s, 'r-', label='M max')
    ax2[1].plot(x_sorted, m_min_s, 'b-', label='M min')
    ax2[1].set_xlabel('Position (m)')
    ax2[1].set_ylabel('Moment (kN.m)')
    ax2[1].legend()
    ax2[1].grid(True)

    plt.suptitle('Shear & Moment Envelope Diagrams')
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    """
    Example usage:
      1) Build a 90 m bridge with interior supports at 30, 60 => three 30 m spans.
      2) Vehicle with 2 axles (20 kN each), 3 m apart.
      3) Trailing axle from x=0 at t=0 to x=bridge_length at t=(bridge_length/speed).
      4) Use np.interp for trailing axle position interpolation.
      5) Step-by-step analysis, plot results.
    """

    # 1) Bridge geometry & properties
    bridge_length     = 90.0
    support_positions = [30, 60]
    E   = 2.0e10
    rho = 2500.0
    A   = [4.5]  # cross-sectional area (m^2)
    I   = [1.5**4 / 12.0]  # second moment of area (m^4)
    damping_ratio = 0.02
    mesh_size = 0.5

    # Build the bridge
    bridge_dict = build_bridge_model(
        bridge_length,
        support_positions,
        E,
        rho,
        A,
        I,
        damping_ratio,
        mesh_size
    )

    # 2) Vehicle config
    vehicle_config = {
        'loads': [20000.0, 20000.0],  # N
        'distances': [0.0, 3.0]       # m
    }

    # 3) Trailing axle path: from (t=0, x=0) to (t=end, x=bridge_length)
    speed = 10.0  # m/s
    trailing_axle_path = [
        [0.0, 0.0],
        [bridge_length/speed, bridge_length]
    ]

    # 4) Run analysis
    dt = 0.01
    results_dict = run_dynamic_analysis_with_path_ts(
        bridge_dict,
        trailing_axle_path,
        vehicle_config,
        dt
    )

    # 5) Post-process & plot
    post_process_and_plot(results_dict)

    # 6) Clean up
    ops.wipe()
