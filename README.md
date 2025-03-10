# 2D Beam Bridge Model in OpenSees with Moving Loads

This repository contains a **simple yet comprehensive** Python script demonstrating how to build a **2D beam bridge model** in [OpenSeesPy](https://openseespydoc.readthedocs.io/en/latest/) and apply **moving loads** (e.g., truck axles) using **Path TimeSeries**. The process covers:

1. **Model Construction**  
   - Nodes and elements (`elasticBeamColumn`) for a beam bridge with multiple supports.  
   - Assigned mass, boundary conditions (hinge, roller), and Rayleigh damping for the first two modes.

2. **Moving Load Computation**  
   - A time-position table for the trailing axle of a vehicle.  
   - Distribution of axle loads to nodes using linear interpolation along the beam.

3. **Path TimeSeries Approach**  
   - Each node receives a time series of vertical loads.  
   - A single transient analysis (`ops.analyze(1, dt)` in a loop) captures displacement, velocity, acceleration, and internal forces (shear, moment) at each time step.

4. **Post-Processing**  
   - Extraction of midspan responses (displacement, velocity, acceleration).  
   - Shear and moment envelopes plotted along the bridge length.

---

## Features

- **Nodal Time Histories**: Pre-computes how each moving axle distributes load to the beam nodes across time.  
- **NumPy Interpolation**: Uses `np.interp` for linear interpolation of trailing axle positions.  
- **Step-by-Step Transient Analysis**: Captures detailed element and node responses without external recorder files.  
- **Straightforward Post-Processing**: Provides envelope diagrams for shear and moment, plus midspan response graphs.  

---

## File Overview

- **`build_bridge_model(...)`**  
  Creates a 2D beam model with user-specified geometry, mesh size, boundary conditions, and Rayleigh damping.

- **`interpolate_position_at_time(t, time_position_table)`**  
  Uses `np.interp` to return the vehicle axle position at any requested time.

- **`get_nodal_forces_for_all_times(time_steps, trailing_axle_path, vehicle_config, node_tags, node_positions)`**  
  Pre-computes vertical loads for each node at all analysis times.

- **`run_dynamic_analysis_with_path_ts(bridge_dict, trailing_axle_path, vehicle_config, dt)`**  
  Sets up a Path TimeSeries for each node, performs a step-by-step transient analysis, and extracts nodal and element responses.

- **`post_process_and_plot(results_dict)`**  
  Plots midspan displacement, velocity, and acceleration time histories. Also plots shear and moment envelopes.

- **`main()`**  
  Demonstrates building a multi-span bridge, defining a vehicle, generating a trailing axle time-position table, running the analysis, and plotting results.

---

## Requirements

- **Python 3.x**  
- **[OpenSeesPy](https://openseespydoc.readthedocs.io/en/latest/)**  
- **NumPy**  
- **Matplotlib**

Install dependencies via pip:
```bash
pip install numpy matplotlib openseespy
