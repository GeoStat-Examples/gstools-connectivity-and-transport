# -*- coding: utf-8 -*-
"""
Transport example using GSTools.

This uses a heterogeneous transmissivity field
with Zinn&Harvey transformation which connects
the highest/lowest values.

Authors: Alraune Zech and Sebastian Müller
"""
import os
import numpy as np
import gstools as gs
from gstools import transform as tf
from ogs5py import OGS, by_id

# set the connectivity for Zinn&Harvey
connectivity = ["mean", "low", "high"]
# rescaling correlation length for zinn&harvey transform (Gong et al. 2013)
length_scales = [1, 1.67, 1.67]
# time settings
hour = 3600
day = 24 * hour
no_days = 30
end_time = no_days * day
step = 2 * hour
# ogs configuration
RES = os.path.join("..", "results")
task_root = os.path.abspath(os.path.join(RES, "tracer_test_conn"))
# initiate new models
model = OGS(task_root=task_root, task_id="model")
# generate a rectangular 2D mesh (x-z cross section):
# ranging from x in [-1,9] with dx=0.025m & z in [-5,-1] with dz=0.1m
model.msh.generate(
    "rectangular",
    dim=2,
    mesh_origin=[-1, -5],
    element_no=[400, 40],
    element_size=[0.025, 0.1],
)
# specify points and lines for boundary conditions and output
model.gli.generate("rectangular", dim=2, ori=(-1, -5), size=(10, 4))
model.gli.add_polyline(points=[0, 1], name="bottom_BC")
model.gli.add_polyline(points=[1, 2], name="right_BC")
model.gli.add_polyline(points=[2, 3], name="top_BC")
model.gli.add_polyline(points=[3, 0], name="left_BC")
model.gli.add_polyline(points=[[0, -2.0, 0], [0, -4.0, 0]], name="source")
# add polylines for observation locations
model.gli.add_polyline(points=[[2.5, -1, 0], [2.5, -5, 0]], name="ow1")
model.gli.add_polyline(points=[[5.5, -1, 0], [5.5, -5, 0]], name="ow2")
model.gli.add_polyline(points=[[8.5, -1, 0], [8.5, -5, 0]], name="ow3")
# create a x-z mesh by swaping y and z axis
model.msh.swap_axis("y", "z")
model.gli.swap_axis("y", "z")
model.pcs.add_block(  # set the process type for flow
    PCS_TYPE="GROUNDWATER_FLOW", NUM_TYPE="NEW", TIM_TYPE="STEADY"
)
model.pcs.add_block(  # set the process type for transport
    PCS_TYPE="MASS_TRANSPORT"
)
for bound in ["bottom_BC", "top_BC"]:
    model.st.add_block(  # neumann boundary condition flow
        PCS_TYPE="GROUNDWATER_FLOW",
        PRIMARY_VARIABLE="HEAD",
        GEO_TYPE=["POLYLINE", bound],
        DIS_TYPE=["CONSTANT_NEUMANN", 0.0],
    )
model.bc.add_block(  # boundary condition flow
    PCS_TYPE="GROUNDWATER_FLOW",
    PRIMARY_VARIABLE="HEAD",
    GEO_TYPE=["POLYLINE", "left_BC"],
    DIS_TYPE=["CONSTANT", 0.01],
)
model.bc.add_block(  # boundary condition flow
    PCS_TYPE="GROUNDWATER_FLOW",
    PRIMARY_VARIABLE="HEAD",
    GEO_TYPE=["POLYLINE", "right_BC"],
    DIS_TYPE=["CONSTANT", 0.0],
)
for bound in ["left_BC", "bottom_BC", "top_BC"]:
    model.bc.add_block(  # boundary condition transport
        PCS_TYPE="MASS_TRANSPORT",
        PRIMARY_VARIABLE="CONCENTRATION1",
        GEO_TYPE=["POLYLINE", bound],
        DIS_TYPE=["CONSTANT", 0.0],
    )
model.bc.add_block(  # source term transport
    PCS_TYPE="MASS_TRANSPORT",
    PRIMARY_VARIABLE="CONCENTRATION1",
    GEO_TYPE=["POLYLINE", "source"],
    DIS_TYPE=["CONSTANT", 1.0],
    TIM_TYPE="CURVE 1",  # from the rfd file
)
model.rfd.add_block(  # Pulse injection
    PROJECT="tracer_test",
    CURVES=zip([0, step, 2 * step, end_time], [0.0, 1.0, 0.0, 0.0]),
)
model.ic.add_block(  # initial condition flow
    PCS_TYPE="GROUNDWATER_FLOW",
    PRIMARY_VARIABLE="HEAD",
    GEO_TYPE="DOMAIN",
    DIS_TYPE=["CONSTANT", 0.0],
)
model.ic.add_block(  # initial condition transport
    PCS_TYPE="MASS_TRANSPORT",
    PRIMARY_VARIABLE="CONCENTRATION1",
    GEO_TYPE="DOMAIN",
    DIS_TYPE=["CONSTANT", 0.0],
)
model.mpd.add(name="transmissivity")
model.mpd.add_block(  # edit recent mpd file
    MSH_TYPE="GROUNDWATER_FLOW",
    MMP_TYPE="PERMEABILITY",
    DIS_TYPE="ELEMENT",
)
model.mmp.add_block(  # medium properties
    GEOMETRY_DIMENSION=2,
    POROSITY=[1, 0.3],
    TORTUOSITY=[1, 1.0],
    STORAGE=[1, 1.0e-04],
    PERMEABILITY_DISTRIBUTION=model.mpd.file_name,
    MASS_DISPERSION=[1, 0.001, 0.001],
)
model.mcp.add_block(  # Component Properties
    NAME="CONCENTRATION1", MOBILE=1, DIFFUSION=[1, 1.0e-08]
)
model.msp.add_block(DENSITY=[1, 2000])  # medium properties
model.tim.add_block(  # set the timesteps for flow
    PCS_TYPE="GROUNDWATER_FLOW",
    TIME_START=0,
    TIME_END=end_time,
    TIME_STEPS=[1, end_time],
)
model.tim.add_block(  # set the timesteps for concentration (2 hourly)
    PCS_TYPE="MASS_TRANSPORT",
    TIME_START=0,
    TIME_END=end_time,
    TIME_STEPS=[end_time // step, step],
)
model.out.add_block(  # domain observation in vtk (15 daily)
    NOD_VALUES=[["HEAD"], ["CONCENTRATION1"]],
    ELE_VALUES=[["VELOCITY1_X"], ["VELOCITY1_Y"], ["VELOCITY1_Z"]],
    GEO_TYPE="DOMAIN",
    DAT_TYPE="VTK",
    TIM_TYPE=["STEPS", end_time // step // 2],
)
model.out.add_block(  # observations at OW1
    PCS_TYPE="MASS_TRANSPORT",
    NOD_VALUES="CONCENTRATION1",
    GEO_TYPE=["POLYLINE", "ow1"],
    DAT_TYPE="TECPLOT",
    TIM_TYPE=["STEPS", 1],
)
model.out.add_block(  # observations at OW2
    PCS_TYPE="MASS_TRANSPORT",
    NOD_VALUES="CONCENTRATION1",
    GEO_TYPE=["POLYLINE", "ow2"],
    DAT_TYPE="TECPLOT",
    TIM_TYPE=["STEPS", 1],
)
model.out.add_block(  # observations at OW3
    PCS_TYPE="MASS_TRANSPORT",
    NOD_VALUES="CONCENTRATION1",
    GEO_TYPE=["POLYLINE", "ow3"],
    DAT_TYPE="TECPLOT",
    TIM_TYPE=["STEPS", 1],
)
model.num.add_block(  # numerical solver flow
    PCS_TYPE="GROUNDWATER_FLOW",
    LINEAR_SOLVER=[2, 5, 1.0e-14, 1000, 1.0, 100, 4],
)
model.num.add_block(  # numerical solver transport
    PCS_TYPE="MASS_TRANSPORT",
    LINEAR_SOLVER=[2, 5, 1.0e-10, 1000, 1.0, 100, 4],
    ELE_GAUSS_POINTS=3,
)
model.write_input()
seed = gs.random.MasterRNG(0)
# iterate over connectivity types
for conn, len_scale in zip(connectivity, length_scales):
    # create the transmissivity field
    cov_model = gs.Gaussian(dim=2, var=2, len_scale=len_scale, anis=0.5)
    srf = gs.SRF(model=cov_model, mean=np.log(1.0e-3), seed=seed())
    # 2d spatial random field in x-z direction on the mesh
    srf.mesh(model.msh, direction="xz")
    # apply Zinn&Harvey transformation
    if conn != "mean":
        tf.zinnharvey(srf, conn=conn)
    # transform to log-normal
    tf.normal_to_lognormal(srf)
    # update the transmissivity the field
    model.mpd.update_block(DATA=by_id(srf.field))
    # write the new mpd file
    model.mpd.write_file()
    model.output_dir = "out_" + conn
    model.msh.export_mesh(  # export field
        os.path.join(model.task_root, "trans_field_" + conn + ".vtk"),
        file_format="vtk",
        cell_data_by_id={"transmissivity": srf.field},
    )
    print("run model for connectivity:", conn)
    success = model.run_model()
    print("success:", success)
