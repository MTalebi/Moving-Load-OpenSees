# 2D Beam Bridge Model in OpenSees with Moving Loads

This repository provides a **thorough** example of modeling a 2D beam bridge in [OpenSeesPy](https://openseespydoc.readthedocs.io) and analyzing its **dynamic response** to moving vehicular loads. The approach combines **finite element discretization**, **Rayleigh damping**, and **Path TimeSeries** to apply time-varying forces from multiple axles.

---

## Overview

1. **Model Construction**  
   - Nodes and `elasticBeamColumn` elements for a 2D Euler–Bernoulli beam.
   - Boundary conditions (hinge, roller), mass assignment, and Rayleigh damping based on the first two modes.

2. **Moving Load Definition**  
   - A time-position table for the trailing axle of the vehicle.
   - Axle loads distributed to beam nodes by linear interpolation at each time step.

3. **Path TimeSeries Approach**  
   - Each node gets a Path TimeSeries with its unique vertical load history.
   - A single transient analysis (`ops.analyze(1, dt)` step-by-step) provides nodal displacements, velocities, accelerations, and internal forces (shear, moment).

4. **Post-Processing**  
   - Extraction of midspan response (displacement, velocity, acceleration).
   - Computation of shear and moment envelopes across the entire bridge length.

---

## File & Function Descriptions

- **`build_bridge_model(...)`**  
  Creates a 2D beam model in OpenSees:  
  + Meshes each span.  
  + Defines elements, supports, nodal mass, and damping (via first two modes).

- **`interpolate_position_at_time(t, time_position_table)`**  
  Leverages `np.interp` for time \(\to\) position interpolation.

- **`get_nodal_forces_for_all_times(...)`**  
  Builds the time history of node loads by distributing each axle load among the nearest nodes at each time step.

- **`run_dynamic_analysis_with_path_ts(...)`**  
  Assigns a Path TimeSeries to each node, then executes a transient analysis step-by-step. Records nodal responses, shear, and moment.

- **`post_process_and_plot(...)`**  
  Plots the midspan time histories (disp, vel, acc) and shear/moment envelopes.

- **`main()`**  
  Demonstrates building a multi-span bridge, defining axle loads, providing the trailing axle time-position table, running the analysis, and visualizing results.

---

## Usage

1. **Install Dependencies**  
   ```bash
   pip install openseespy numpy matplotlib
   ```

2. **Adjust Parameters** in the `main()` function:
   + Bridge length, supports, cross-sectional properties, damping ratio.
   + Axle loads, spacing, time-position data.
   + Time step `dt`.

3. **Run** the script (e.g. `python code.py`).

4. **Plot & Analyze** the resulting midspan response and envelope diagrams.

---

## Citation

If you find this code beneficial to your research or practical projects, please use the following reference:

> **Talebi-Kalaleh, M. (2025).** *Modeling of 2D Moving Loads using OpenSeesPy [Source code].* GitHub.  
> <https://github.com/MTalebi/Moving_Load_OpenSees> *(Published March 10, 2025).*

---

**Thank you** for using and sharing this 2D moving-load OpenSees example!
