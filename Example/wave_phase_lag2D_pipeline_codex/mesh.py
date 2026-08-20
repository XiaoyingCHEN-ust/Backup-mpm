import argparse
import json
from pathlib import Path

PARSER = argparse.ArgumentParser(
    description="Generate the shallow-buried rigid-pipeline MPM geometry")
PARSER.add_argument("--mesh-size", type=float, default=0.02)
PARSER.add_argument("--particle-spacing", type=float, default=0.01)
PARSER.add_argument("--cover-over-diameter", type=float, default=0.25)
PARSER.add_argument("--output-dir", type=Path, default=Path("."))
PARSER.add_argument("--show", action="store_true")
ARGS = PARSER.parse_args()

# Keep server-side generation independent of Matplotlib. Pass --show on a
# workstation only when a geometry preview is wanted.
SHOW_PLOT = ARGS.show
import numpy

if SHOW_PLOT:
    try:
        import matplotlib
        import matplotlib.pyplot as plt
    except ModuleNotFoundError as exception:
        if exception.name != "matplotlib":
            raise
        raise RuntimeError(
            "--show requires Matplotlib; install it or run without --show"
        ) from exception


SCRIPT_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = ARGS.output_dir
if not OUTPUT_DIR.is_absolute():
    OUTPUT_DIR = SCRIPT_DIR / OUTPUT_DIR
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

def mpm_mesh2d(xmin=0.0, xmax=1.0, ymin=0.0, ymax=1.0, sizex=0.02, sizey=0.02):
    #  4-node Quadrilateral Element
    # ! 3 0----------0 2
    # !   |          |
    # !   |          |
    # !   |          |
    # !   |          |
    # ! 0 0----------0 1

    # Generate suitable ranges for parametrization
    x_range = numpy.arange(xmin, xmax+1.e-15, sizex)
    y_range = numpy.arange(ymin, ymax+1.e-15, sizey)

    nx=x_range.shape[0]
    ny=y_range.shape[0]

    # Create the vertices.
    x, y = numpy.meshgrid(x_range, y_range, indexing="ij")
    nodes = numpy.array([x, y]).T.reshape(-1, 2)

    v0 = numpy.add.outer(numpy.array(range(nx - 1)), nx * numpy.array(range(ny - 1))).T.flatten()

    elem = numpy.array([v0,v0,v0,v0],dtype='int').transpose()

    elem[:, 0] += 0
    elem[:, 1] += 1
    elem[:, 2] += nx+1
    elem[:, 3] += nx

    return nodes, elem



def gimp_mesh2d(xmin=0.0, xmax=1.0, ymin=0.0, ymax=1.0, sizex=0.1, sizey=0.1):
    # #16 Nodes GIMP
    # //!   13----------12----------11----------10
    # //!   |           |           |           |
    # //!   |           |           |           |
    # //!   |           |           |           |
    # //!   |        (-1, 1)      (1,1)         |
    # //!   14----------3-----------2-----------9
    # //!   |           |           |           |
    # //!   |           | particle  |           |
    # //!   |           | location  |           |
    # //!   |           |           |           |
    # //!   15----------0-----------1-----------8
    # //!   |        (-1,-1)      (1,-1)        |
    # //!   |           |           |           |
    # //!   |           |           |           |
    # //!   |           |           |           |
    # //!   4-----------5-----------6-----------7

    # Creat2 additional node for gimp, before start and after end
    x_range = numpy.arange(xmin-sizex, xmax+sizex+1.e-15, sizex)
    y_range = numpy.arange(ymin-sizex, ymax+sizey+1.e-15, sizey)

    nx=x_range.shape[0]
    ny=y_range.shape[0]

    # Create the vertices.
    x, y = numpy.meshgrid(x_range, y_range , indexing="ij")
    nodes = numpy.array([x, y]).T.reshape(-1, 2)

    # Get the index of the start point of each group
    v0 = numpy.add.outer(numpy.array(range(1, nx - 2)), nx * numpy.array(range(1, ny - 2))).T.flatten()
    v1 = numpy.add.outer(numpy.array(range(0, nx - 3)), nx * numpy.array(range(0, ny - 3))).T.flatten()

    # Initialize each goup with start point id
    group0 = numpy.repeat([v0],4,axis=0)
    group1 = numpy.repeat([v1],12,axis=0)

    # Concatenate them
    elem=numpy.concatenate((group0,group1), axis=0).transpose()

    #First group
    elem[:,0] += 0
    elem[:,1] += 1
    elem[:,2] += nx+1
    elem[:,3] += nx

    #Second group
    elem[:,4] += 0
    elem[:,5] += 1
    elem[:,6] += 2
    elem[:,7] += 3
    elem[:,8] += nx+3
    elem[:,9] += 2*nx +3
    elem[:,10] += 3*nx +3
    elem[:,11] += 3*nx +2
    elem[:,12] += 3*nx +1
    elem[:,13] += 3*nx
    elem[:,14] += 2*nx
    elem[:,15] += 1*nx

    return nodes, elem


def mpm_mesh3d(xmin=0.0, xmax=1.0, ymin=0.0, ymax=1.0, zmin=0.0, zmax=1.0, sizex=1, sizey=1, sizez=1):
    # Hexahedron:
    #        v
    # 3----------2
    # |\     ^   |\
    # | \    |   | \
    # |  \   |   |  \
    # |   7------+---6
    # |   |  +-- |-- | -> u
    # 0---+---\--1   |
    #  \  |    \  \  |
    #   \ |     \  \ |
    #    \|      w  \|
    #     4----------5

    # Generate suitable ranges for parametrization
    x_range = numpy.arange(xmin, xmax+1.e-15, sizex)
    y_range = numpy.arange(ymin, ymax+1.e-15, sizey)
    z_range = numpy.arange(zmin, zmax+1.e-15, sizez)

    nx=x_range.shape[0]
    ny=y_range.shape[0]
    nz=z_range.shape[0]

    # Create the vertices.
    x, y, z = numpy.meshgrid(x_range, y_range, z_range, indexing="ij")
    nodes = numpy.array([x, y, z]).T.reshape(-1, 3)

    # Compute vertex id
    temp_v0 = numpy.add.outer(numpy.array(range(nx - 1)), nx * numpy.array(range(ny - 1)))
    v0 = numpy.add.outer(temp_v0, nx * ny * numpy.array(range(nz - 1))).T.flatten()

    temp_v3 = numpy.add.outer(numpy.array(range(nx-1)), nx * numpy.array(range(1, ny)))
    v3 = numpy.add.outer(temp_v3, nx * ny * numpy.array(range(nz-1))).T.flatten()

    temp_v4 = numpy.add.outer(numpy.array(range(nx-1)), nx * numpy.array(range(ny-1)))
    v4 = numpy.add.outer(temp_v4, nx * ny * numpy.array(range(1,nz))).T.flatten()

    temp_v7 = numpy.add.outer(numpy.array(range(nx-1)), nx * numpy.array(range(1,ny)))
    v7 = numpy.add.outer(temp_v7, nx * ny * numpy.array(range(1,nz))).T.flatten()

    elem = numpy.array([v0,v0,v3,v3,v4,v4,v7,v7],dtype='int').transpose()
    elem[:,1::4] +=1
    elem[:,2::4] +=1

    return nodes, elem


def gimp_mesh3d(xmin=0.0, xmax=1.0, ymin=0.0, ymax=1.0, zmin=0.0, zmax=1.0, sizex=1, sizey=1, sizez=1):

    # Creat additional node for gimp, before start and after end
    x_range = numpy.arange(xmin-sizex, xmax+sizex+1.e-15, sizex)
    y_range = numpy.arange(ymin-sizex, ymax+sizey+1.e-15, sizey)
    z_range = numpy.arange(zmin-sizez, zmax+sizez+1.e-15, sizez)

    nx=x_range.shape[0]
    ny=y_range.shape[0]
    nz=z_range.shape[0]

    # Create the vertices.
    x, y, z = numpy.meshgrid(x_range, y_range, z_range, indexing="ij")
    nodes = numpy.array([x, y, z]).T.reshape(-1, 3)

    # Get the index of the start point of each group
    temp_v0 = numpy.add.outer(numpy.array(range(1,nx-2)), nx * numpy.array(range(2,ny-1)))
    v0 = numpy.add.outer(temp_v0, nx * ny * numpy.array(range(1,nz-2))).T.flatten()

    temp_v1 = numpy.add.outer(numpy.array(range(0,nx-3)), nx * numpy.array(range(3,ny)))
    v1 = numpy.add.outer(temp_v1, nx * ny * numpy.array(range(0,nz-3))).T.flatten()

    temp_v2 = numpy.add.outer(numpy.array(range(0,nx-3)), nx * numpy.array(range(2,ny-1)))
    v2 = numpy.add.outer(temp_v2, nx * ny * numpy.array(range(0,nz-3))).T.flatten()

    temp_v3 = numpy.add.outer(numpy.array(range(0,nx-3)), nx * numpy.array(range(1,ny-2)))
    v3 = numpy.add.outer(temp_v3, nx * ny * numpy.array(range(0,nz-3))).T.flatten()

    temp_v4 = numpy.add.outer(numpy.array(range(0,nx-3)), nx * numpy.array(range(0,ny-3)))
    v4 = numpy.add.outer(temp_v4, nx * ny * numpy.array(range(0,nz-3))).T.flatten()


    # Initialize each goup with start point id
    group0 = numpy.repeat([v0],8,axis=0)
    group1 = numpy.repeat([v1],16,axis=0)
    group2 = numpy.repeat([v2],12,axis=0)
    group3 = numpy.repeat([v3],12,axis=0)
    group4 = numpy.repeat([v4],16,axis=0)

    # Concatenate them
    elem=numpy.concatenate((group0,group1,group2,group3,group4), axis=0).transpose()

    #First group
    elem[:,1] +=1
    elem[:,2] +=nx*ny+1
    elem[:,3] +=nx*ny
    elem[:,4] += -nx
    elem[:,5] += -nx+1
    elem[:,6] += -nx+1+nx*ny
    elem[:,7] += -nx+nx*ny

    #Second group
    elem[:,8] +=0
    elem[:,9] +=1
    elem[:,10] +=2
    elem[:,11] +=3
    elem[:,12] += nx*ny
    elem[:,13] += nx*ny+1
    elem[:,14] += nx*ny+2
    elem[:,15] += nx*ny+3
    elem[:,16] += 2*nx*ny
    elem[:,17] += 2*nx*ny+1
    elem[:,18] += 2*nx*ny+2
    elem[:,19] += 2*nx*ny+3
    elem[:,20] += 3*nx*ny
    elem[:,21] += 3*nx*ny+1
    elem[:,22] += 3*nx*ny+2
    elem[:,23] += 3*nx*ny+3

    #Third group
    elem[:,24] +=0
    elem[:,25] +=1
    elem[:,26] +=2
    elem[:,27] +=3
    elem[:,28] += nx*ny
    elem[:,29] += nx*ny+3
    elem[:,30] += 2*nx*ny
    elem[:,31] += 2*nx*ny+3
    elem[:,32] += 3*nx*ny
    elem[:,33] += 3*nx*ny+1
    elem[:,34] += 3*nx*ny+2
    elem[:,35] += 3*nx*ny+3

    #Fourth group
    elem[:,36] +=0
    elem[:,37] +=1
    elem[:,38] +=2
    elem[:,39] +=3
    elem[:,40] += nx*ny
    elem[:,41] += nx*ny+3
    elem[:,42] += 2*nx*ny
    elem[:,43] += 2*nx*ny+3
    elem[:,44] += 3*nx*ny
    elem[:,45] += 3*nx*ny+1
    elem[:,46] += 3*nx*ny+2
    elem[:,47] += 3*nx*ny+3

    #Fifth group
    elem[:,48] +=0
    elem[:,49] +=1
    elem[:,50] +=2
    elem[:,51] +=3
    elem[:,52] += nx*ny
    elem[:,53] += nx*ny+1
    elem[:,54] += nx*ny+2
    elem[:,55] += nx*ny+3
    elem[:,56] += 2*nx*ny
    elem[:,57] += 2*nx*ny+1
    elem[:,58] += 2*nx*ny+2
    elem[:,59] += 2*nx*ny+3
    elem[:,60] += 3*nx*ny
    elem[:,61] += 3*nx*ny+1
    elem[:,62] += 3*nx*ny+2
    elem[:,63] += 3*nx*ny+3

    return nodes, elem

#Save nodes coordinate up to 5 dicimal digit
#Save elem node

def save_mesh(nodes, elem, name):
    shape={4:'quadrilateral',16:'quadrilateral',8:'hexahedron',64:'hexahedron'}
    mesh_file = OUTPUT_DIR / (name+".txt")
    f_m = open(mesh_file, "w")

    head = '#! elementShape {}\n#! elementNumPoints {}\n{} {}\n'.format(shape[elem.shape[1]],
                                 elem.shape[1], nodes.shape[0], elem.shape[0])
    f_m.write(head)
    numpy.savetxt(f_m, nodes, fmt="%.5f", delimiter="\t")
    numpy.savetxt(f_m, elem, fmt="%i", delimiter="\t")
    f_m.close()


def cube_particle2d(xmin=0.0, xmax=1.0, ymin=0.0, ymax=1.0, size=1.):

    # Generate suitable ranges for parametrization
    x_range = numpy.arange(xmin+0.5*size, xmax, size)
    y_range = numpy.arange(ymin+0.5*size, ymax, size)

    # Create the vertices.
    x, y = numpy.meshgrid(x_range, y_range, indexing="ij")
    particles = numpy.array([x, y]).T.reshape(-1, 2)

    return particles


def cube_particle3d(xmin=0.0, xmax=1.0, ymin=0.0, ymax=1.0, zmin=0.0, zmax=1.0, size=1.):

    # Generate suitable ranges for parametrization
    x_range = numpy.arange(xmin+0.5*size, xmax, size)
    y_range = numpy.arange(ymin+0.5*size, ymax, size)
    z_range = numpy.arange(zmin+0.5*size, zmax, size)

    # Create the vertices.
    x, y, z = numpy.meshgrid(x_range, y_range, z_range, indexing="ij")
    particles = numpy.array([x, y, z]).T.reshape(-1, 3)

    return particles


def cylinder_particle3d(center=[0.,0.,0.],r=1,h=1, axis=2, size=1.):
    xmin = center[0]-r
    xmax = center[0]+r
    ymin = center[1]-r
    ymax = center[1]+r
    zmin = center[2]
    zmax = center[2]+h

    # Generate suitable ranges for parametrization
    x_range = numpy.arange(xmin+0.5*size, xmax, size)
    y_range = numpy.arange(ymin+0.5*size, ymax, size)
    z_range = numpy.arange(zmin+0.5*size, zmax, size)

    # Create the vertices.
    x, y, z = numpy.meshgrid(x_range, y_range, z_range, indexing="ij")
    particles_temp = numpy.array([x, y, z]).T.reshape(-1, 3)

    if axis == 0:
        index = numpy.where(numpy.linalg.norm(particles_temp[:,[1,2]],axis=1)<(numpy.sqrt(r*r)+1.e-15))
    elif axis == 1:
        index = numpy.where(numpy.linalg.norm(particles_temp[:,[0,2]],axis=1)<(numpy.sqrt(r*r)+1.e-15))
    elif axis == 2:
        index = numpy.where(numpy.linalg.norm(particles_temp[:,[0,1]],axis=1)<(numpy.sqrt(r*r)+1.e-15))
    else:
        print ("axis should be 0,1 or 2.")
    return particles_temp[index]

def cylinder_quarter_particle3d(center=[0.,0.,0.],r=1,h=1, axis=2, size=1.):
    particles_temp = cylinder_particle(center=center,r=r,h=h, axis=axis, size=size)

    if axis == 0:
        index = numpy.where((particles_temp[:,1]>-1.e-15) & (particles_temp[:,2]>-1.e-15) )
    elif axis == 1:
        index = numpy.where((particles_temp[:,0]>-1.e-15) & (particles_temp[:,2]>-1.e-15) )
    elif axis == 2:
        index = numpy.where((particles_temp[:,0]>-1.e-15) & (particles_temp[:,1]>-1.e-15) )
    else:
        print ("axis should be 0,1 or 2.")

    return particles_temp[index]



#Save particle coordinate up to 5 dicimal digit
def save_particle(particles, name, fmt="%.6f"):
    file_name = OUTPUT_DIR / (name+".txt")
    f = open(file_name, "w")
    head = '{}\n'.format(particles.shape[0])
    f.write(head)
    numpy.savetxt(f, particles, fmt=fmt, delimiter="\t")
    f.close()

def save_scalars(values, name):
    """Write the one-value-per-particle format expected by IOMeshAscii."""
    values = numpy.asarray(values, dtype=float)
    if values.ndim != 1 or not numpy.all(numpy.isfinite(values)):
        raise ValueError("{} must be a finite one-dimensional array".format(name))
    file_name = OUTPUT_DIR / (name + ".txt")
    with open(file_name, "w") as f:
        f.write('{}\n'.format(values.shape[0]))
        numpy.savetxt(f, values, fmt="%.10f")

def save_stresses(stresses, name):
    file_name = OUTPUT_DIR / (name + ".txt")
    f = open(file_name, "w")
    head = '{}\n'.format(stresses.shape[0])
    f.write(head)
    numpy.savetxt(f, stresses, fmt="%.10f", delimiter="\t")
    f.close()


def save_entity_sets(node_sets, particle_sets):
    """Write generated ids with each ``set`` kept on one horizontal line."""
    groups = (
        ("node_sets", node_sets),
        ("particle_sets", particle_sets),
    )
    with open(OUTPUT_DIR / "entity_sets.json", "w") as f:
        f.write("{\n")
        for group_index, (group_name, entries) in enumerate(groups):
            f.write('    "{}": [\n'.format(group_name))
            for entry_index, (set_id, ids) in enumerate(entries):
                values = numpy.asarray(ids, dtype=int).reshape(-1).tolist()
                f.write("        {\n")
                f.write('            "id": {},\n'.format(int(set_id)))
                f.write('            "set": {}\n'.format(json.dumps(values)))
                suffix = "," if entry_index + 1 < len(entries) else ""
                f.write("        }}{}\n".format(suffix))
            suffix = "," if group_index + 1 < len(groups) else ""
            f.write("    ]{}\n".format(suffix))
        f.write("}\n")


# Function to rotate particles by theta degrees
def rotate_particles(particles, theta):
    # Convert theta from degrees to radians
    theta_rad = numpy.deg2rad(theta)
    # Rotation matrix
    rotation_matrix = numpy.array([
        [numpy.cos(theta_rad), -numpy.sin(theta_rad)],
        [numpy.sin(theta_rad), numpy.cos(theta_rad)]
    ])
    # Apply rotation to each particle
    rotated_particles = numpy.dot(particles, rotation_matrix.T)
    return rotated_particles


# Model discretisation. Keep two particles per background cell in each
# direction when changing the resolution.
mesh_size = ARGS.mesh_size
particle_spacing = ARGS.particle_spacing
if mesh_size <= 0.0 or particle_spacing <= 0.0:
    raise ValueError("mesh and particle spacing must be positive")
particles_per_direction = mesh_size / particle_spacing
if not numpy.isclose(particles_per_direction, round(particles_per_direction)):
    raise ValueError("mesh size must be an integer multiple of particle spacing")

# MPM background mesh
nodes,elem = mpm_mesh2d(xmin=0, xmax=1.5, ymin=0, ymax=0.6,
                        sizex=mesh_size, sizey=mesh_size)
save_mesh(nodes, elem, "mesh")

# Pipeline outer geometry.  The primary study uses cover/D = 0.25 so the pipe
# invert intersects the approximately 0.14--0.15 m phase-lag-sensitive layer
# identified in the preceding manuscript section.  The circular rigid body is
# created independently by the solver JSON; these files document/check the
# identical geometry and provide the initial cavity/contact particle sets.
tunnel_radius = 0.06
seabed_elevation_geometry = 0.5
pipeline_diameter = 2.0 * tunnel_radius
pipeline_cover = ARGS.cover_over_diameter * pipeline_diameter
if pipeline_cover < 0.0:
    raise ValueError("cover/D must be non-negative")
tunnel_center = numpy.array([
    0.7, seabed_elevation_geometry - tunnel_radius - pipeline_cover])
pipeline_cover_over_diameter = pipeline_cover / pipeline_diameter
pipeline_invert_depth = seabed_elevation_geometry - (
    tunnel_center[1] - tunnel_radius)
pipeline_surface_nsegments = 64
pipeline_wall_thickness = 0.12 / 17.0
pipeline_density = 959.0
pipeline_friction = 0.15
pipeline_particle_radius = 0.5 * particle_spacing

pipeline_angles = numpy.linspace(
    0.0, 2.0 * numpy.pi, pipeline_surface_nsegments, endpoint=False)
pipeline_surface_nodes = tunnel_center + tunnel_radius * numpy.column_stack(
    (numpy.cos(pipeline_angles), numpy.sin(pipeline_angles)))

with open(OUTPUT_DIR / "pipeline_geometry.json", "w") as f:
    json.dump({
        "description": "Coupled circular rigid PE100 SDR17 pipeline",
        "center": tunnel_center.tolist(),
        "outer_radius": tunnel_radius,
        "diameter": pipeline_diameter,
        "cover": pipeline_cover,
        "cover_over_diameter": pipeline_cover_over_diameter,
        "invert_depth": pipeline_invert_depth,
        "wall_thickness": pipeline_wall_thickness,
        "density": pipeline_density,
        "contents": "empty",
        "friction_coefficient": pipeline_friction,
        "surface_segments": pipeline_surface_nsegments,
        "background_cell_size": mesh_size,
        "soil_particle_spacing": particle_spacing,
        "connected_to_mpm_solver": True
    }, f, indent=4)
    f.write("\n")

# # gimp 2d mesh
# nodes,elem = gimp_mesh2d(xmin=0.0, xmax=3.0, ymin=0.0, ymax=0.6, sizex=0.02, sizey=0.02)
# save_mesh(nodes, elem, "gimp_mesh2d")

# Find boundary node sets
bottom_node_ids = numpy.where(numpy.isclose(nodes[:, 1], 0.0))[0]

top_node_ids = numpy.where(numpy.isclose(nodes[:, 1], 0.6))[0]

left_node_ids = numpy.where(numpy.isclose(nodes[:, 0], 0.0))[0]

right_node_ids = numpy.where(numpy.isclose(nodes[:, 0], 1.5))[0]

# Create soil particles and remove the empty pipe cavity. The half-particle
# clearance makes each particle domain initially tangent to the pipe surface.
all_particles = cube_particle2d(
    xmin=0, xmax=1.5, ymin=0, ymax=0.5, size=particle_spacing)
all_pipeline_distances = numpy.linalg.norm(
    all_particles - tunnel_center, axis=1)
initial_contact_distance = tunnel_radius + pipeline_particle_radius
particles = all_particles[
    all_pipeline_distances >= initial_contact_distance - 1.e-15]
pipeline_distances = numpy.linalg.norm(particles - tunnel_center, axis=1)
tunnel_boundary_particle_ids = numpy.where(
    pipeline_distances <= initial_contact_distance + particle_spacing + 1.e-15
)[0]
save_particle(particles, "particles")

# Find particle sets. The top set contains two particle rows, matching the
# relative support used by the original 0.01/0.005 discretisation.
top_particle_ids = numpy.where(
    particles[:, 1] >= 0.5 - 2.0 * particle_spacing - 1.e-15)[0]

bottom_particle_ids_1 = numpy.where(
    numpy.isclose(particles[:, 1], 0.5 * particle_spacing))[0]

bottom_particle_ids_2 = numpy.where(
    numpy.isclose(particles[:, 1], 1.5 * particle_spacing))[0]

# The moving circular no-flux condition is applied by the rigid-pipeline solver
# directly to adjacent phase particles. No fixed stair-step background-node
# wall is generated here.
fixed_node_ids = bottom_node_ids
save_entity_sets(
    node_sets=[
        (0, fixed_node_ids),
        (1, left_node_ids),
        (2, top_node_ids),
        (3, right_node_ids),
    ],
    particle_sets=[
        (0, top_particle_ids),
        (1, bottom_particle_ids_1),
        (2, bottom_particle_ids_2),
        (3, tunnel_boundary_particle_ids),
    ])

temperatures = numpy.zeros(particles.shape[0])
save_scalars(temperatures, "initial_temperature")

# Hydrostatic pressure is stored as liquid pressure. ThreePhaseParticleLag
# combines it with the saturation-dependent suction after its liquid material
# has been initialised. Keeping every particle on this profile avoids the
# artificial bottom-only pressure jump that destabilises nearly saturated soil.
soil_density = 2650.0
water_density = 1000.0
porosity = 0.485
gravity = 9.81
sea_level = 1.0
depth_left = 0.5
depth_right = 0.5
liquid_pressure = water_density * gravity * numpy.maximum(sea_level - particles[:, 1], 0.0)
save_scalars(liquid_pressure, "initial_liquid_pressures")

# Initial skeleton stresses for the PIC=1, LinearElastic2D stabilisation stage.
# The solver stores compression as negative and subtracts pore pressure when it
# assembles total mixture stress, so these files must contain EFFECTIVE stress,
# not total stress.
water_depth = depth_left + (depth_right - depth_left) * (particles[:, 0] / 1.5)
seabed_elevation = sea_level - water_depth
burial_depth = numpy.maximum(seabed_elevation - particles[:, 1], 0.0)
gas_molar_mass = 0.029
gas_constant = 8.314
reference_pressure = 100000.0
reference_temperature = 273.15
gas_density = (
    gas_molar_mass * reference_pressure / gas_constant / reference_temperature
)
liquid_saturation = 0.94
K0_effective = 0.72
gas_saturation = 1.0 - liquid_saturation
mixture_density = (
    (1.0 - porosity) * soil_density
    + porosity * (liquid_saturation * water_density + gas_saturation * gas_density)
)
effective_unit_weight = (mixture_density - water_density) * gravity
if effective_unit_weight <= 0.0:
    raise RuntimeError("Initial effective unit weight must be positive")
sigma_yy = -effective_unit_weight * burial_depth
stresses = numpy.zeros((particles.shape[0], 6))
stresses[:, 0] = K0_effective * sigma_yy
stresses[:, 1] = sigma_yy
save_stresses(stresses, "initial_effective_stresses_HS")

# Basic consistency checks catch stale ids before a long server run.
for name, ids, upper_bound in (
        ("fixed nodes", fixed_node_ids, nodes.shape[0]),
        ("left nodes", left_node_ids, nodes.shape[0]),
        ("top nodes", top_node_ids, nodes.shape[0]),
        ("right nodes", right_node_ids, nodes.shape[0]),
        ("top particles", top_particle_ids, particles.shape[0]),
        ("bottom particles 1", bottom_particle_ids_1, particles.shape[0]),
        ("bottom particles 2", bottom_particle_ids_2, particles.shape[0]),
        ("pipeline boundary particles", tunnel_boundary_particle_ids,
         particles.shape[0])):
    ids = numpy.asarray(ids, dtype=int).reshape(-1)
    if ids.size == 0 or ids.min() < 0 or ids.max() >= upper_bound:
        raise RuntimeError("Invalid or empty {} set".format(name))

pipeline_radii = numpy.linalg.norm(
    pipeline_surface_nodes - tunnel_center, axis=1)
maximum_radius_error = float(numpy.max(numpy.abs(
    pipeline_radii - tunnel_radius)))
if maximum_radius_error > 1.e-12:
    raise RuntimeError("Smooth pipeline surface radius check failed")

summary = {
    "mesh_size": mesh_size,
    "particle_spacing": particle_spacing,
    "output_directory": str(OUTPUT_DIR),
    "nodes": int(nodes.shape[0]),
    "cells": int(elem.shape[0]),
    "particles": int(particles.shape[0]),
    "particles_per_cell": int(round((mesh_size / particle_spacing) ** 2)),
    "soil_porosity": porosity,
    "empty_pipeline_cavity_particles_removed": int(
        all_particles.shape[0] - particles.shape[0]),
    "pipeline_contact_particles": int(tunnel_boundary_particle_ids.shape[0]),
    "pipeline_wall_flow_condition": (
        "moving circular particle-level zero normal liquid/gas velocity"
    ),
    "pipeline_wall_thickness": pipeline_wall_thickness,
    "pipeline_density": pipeline_density,
    "pipeline_diameter": pipeline_diameter,
    "pipeline_cover": pipeline_cover,
    "pipeline_cover_over_diameter": pipeline_cover_over_diameter,
    "pipeline_invert_depth": pipeline_invert_depth,
    "pipeline_friction_coefficient": pipeline_friction,
    "smooth_pipeline_surface_segments": pipeline_surface_nsegments,
    "smooth_pipeline_max_radius_error": maximum_radius_error,
    "initial_liquid_pressure_file": "initial_liquid_pressures.txt",
    "initial_effective_stress_file": "initial_effective_stresses_HS.txt",
    "initial_effective_unit_weight_n_per_m3": effective_unit_weight,
    "initial_effective_stress_K0": K0_effective,
}
with open(OUTPUT_DIR / "generation_summary.json", "w") as f:
    json.dump(summary, f, indent=4)
    f.write("\n")

# Save a preview only when explicitly requested with `python mesh.py --show`.
if SHOW_PLOT:
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    tunnel_boundary_particles = particles[tunnel_boundary_particle_ids]
    for ax in axes:
        ax.scatter(particles[:, 0], particles[:, 1], s=1, color="0.75",
                   label="MPM particles")
        ax.scatter(tunnel_boundary_particles[:, 0],
                   tunnel_boundary_particles[:, 1], s=8,
                   color="tab:orange", label="contact particles")
        closed_surface = numpy.vstack(
            (pipeline_surface_nodes, pipeline_surface_nodes[0]))
        ax.plot(closed_surface[:, 0], closed_surface[:, 1], color="tab:blue",
                linewidth=2, label="smooth pipeline surface")
        ax.set_aspect("equal", adjustable="box")
        ax.grid(True, linewidth=0.3)

    axes[0].set_xlim(0.0, 1.5)
    axes[0].set_ylim(0.0, 0.6)
    axes[0].set_title("Coarse MPM domain")
    axes[1].set_xlim(tunnel_center[0] - 0.1, tunnel_center[0] + 0.1)
    axes[1].set_ylim(tunnel_center[1] - 0.1, tunnel_center[1] + 0.1)
    axes[1].set_title("Pipeline geometry close-up")
    axes[1].legend(loc="upper right", fontsize=8)
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "mesh_preview.png", dpi=200)
    plt.show()
    plt.close(fig)

print(json.dumps(summary, indent=2))
