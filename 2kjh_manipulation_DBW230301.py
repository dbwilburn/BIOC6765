
## BIOC6765
## Damien Beau Wilburn 
## https://github.com/dbwilburn/BIOC6765
## Molecular visualization
## Example script using Python API for PyMOL


# Library import
import pymol
from pymol import cmd

def save_png( file_name, width=2400, height=2400, dpi=600, ray=False, ):
    '''
    PyMOL is configured to run commands as fast as possible, such that
    when running slow commands (e.g. writing image files that have to
    access the hard drive), PyMOL can skip ahead and run subsequent
    commands before an image has had time to be properly rendered + written.
    As such, this is a simple wrapper function that will force it to slow
    down and complete a job before going to the next one
    '''
    cmd.refresh()
    cmd.png( file_name, width, height, dpi, ray, quiet=0, )
    cmd.refresh()

# Initiate a PyMOL GUI session
pymol.finish_launching()

# The majority of basic PyMOL commands are implemented within the cmd module
# List of all commands can be found in online documentation:
# https://pymolwiki.org/index.php/Category:Commands
# https://pymol.org/dokuwiki/doku.php?id=api 

# Fetch the model from the PDB
cmd.fetch( '2kjh' )
save_png( '0_fetch-data.png',  )

# Create selections for Ub and UbcH8
cmd.select( 'ubch8', 'chain A', )
cmd.select( 'ub', 'chain B', )

# Select UbcH8 cysteines
cmd.select( 'ubch8_allcys', 'ubch8 and resn cys', )

# Color each chain a different color
cmd.color( 'firebrick', 'ub', )
cmd.color( 'density', 'ubch8')
save_png( '1a_colored.png' )
cmd.bg_color( 'white' )
save_png( '1b_colored_bgwhite.png' )


# Show sticks for Ub
cmd.show( 'sticks', 'ub', )
save_png( '2a_ub-with-sticks.png' )
cmd.select( 'ub_c', 'ub and resi 76' )
cmd.hide( 'sticks', 'ub and not ub_c', )
save_png( '2b_ub-with-C76-sticks.png' )

# Zoom in on Ub C-terminal Cys
cmd.zoom( 'ub and resi 76', 8, )
save_png( '3a_ub-c76.png' )
cmd.select( 'ubch8_c', 'ubch8 and resn cys and (ub_c around 6)' )
cmd.select( 'het_disulf', 'ub_c + ubch8_c' )
cmd.show( 'spheres', 'het_disulf', )
cmd.color( 'yellow', 'het_disulf', )
save_png( '3b_hetero-disulfide.png', )

# De-select all atoms
cmd.deselect()

# Rotate the image to see more clearly
cmd.rotate( 'y', '30', )
save_png( '4a_rotate-around-hetdisulf.png', )
cmd.hide( 'sticks' )
cmd.hide( 'spheres' )
cmd.show( 'sticks', 'het_disulf and not name n+c+o and not elem h' )
save_png( '4b_clean-disulf.png', )
cmd.set( 'all_states', 1, )
save_png( '4b_clean-disulf-allstates.png', )

# Add surface visualization
cmd.set( 'all_states', 0, )
cmd.show( 'surface' )
save_png( '5a_surface-naive.png', )
cmd.set( 'surface_quality', 1, )
save_png( '5b_surface-hq.png', )
cmd.set( 'transparency', 0.6, )
save_png( '5c_surface-60p-transparent.png', )

# Finishing touches
save_png( '6a_ray-default.png', ray=True, )
cmd.set( 'ray_shadows', 0, )
cmd.set( 'spec_reflect', 0, )
save_png( '6b_ray-noreflect.png', ray=True, )
cmd.set( 'ray_trace_mode', 1 )
cmd.set( 'ray_trace_color', 'black', )
cmd.set( 'ray_trace_gain', 2.0, )
save_png( '6c_ray-outline.png', ray=True, )



