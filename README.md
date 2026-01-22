This is data and code for the paper "Biglobal input–output analysis of separated flow over a periodic hill".

Dedalus-VPM is a structured, 3D CFD solver that solves the incompressible navier stokes equations, using spectral discretization in X and Z, while using Chebyshev discretization in Y direction. The volume penalty method has been implemented with the Dedalus DNS solver to compute flow over a periodic hill. Dedalus-VPM code for periodic hill: main_channelflow_volumepenalty_symmetry_nek5000_1b_bulk_velocity2D.py
To install Dedalus, please go to the following link and follow the steps: https://dedalus-project.readthedocs.io/en/latest/pages/installation.html

Dedalus-VPM post processing code: dedalus_flow_output_data_channeflow_velocityNe5000_2D_jfm.m

Mean velocity of Periodic hill functio: mean_velocity_dedalus_postprocess.m

Nek Output Mean flow computation code and relative error with Dedalus : hillp_outlet_flow2_jfmlegend_Re100.m

Input-Output code for flow over periodic hill: linear_stability_Hillp_2Dmeanflow_HPC_quick2.m. fourdiff.m, chebdiff.m and clencurt.m are needed to run the input-output code

Signed distance function code for input-output analysis over periodic hill: main_signeddistfunctn_inputoutput_hillp, This code needs the matlab function normalized_mask_optimized_peclet.m , to compute the numerical interface thickness for VPM in input-output analysis. 
data_x and data_y arrays or .csv files needed to read the dedalus grid. Other inputs include Re, Ny and Nx.

Input-Output post processing file for periodic hill:  input_output_post_process_2D_hillp_3Dec2026.m

Input-output post processing codes to compute grid independence: input_output_post_process_2D_hillp_mesh_independence.m , read_mesh_indepence_files_2.m
