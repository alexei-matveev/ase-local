"""
Module for povray file format support.

See http://www.povray.org/ for details on the format.
"""
import os

import numpy as npy

from ase.io.eps import EPS
from ase.data import chemical_symbols

def pa(array):
    """Povray array syntax"""
    return '<% 6.2f, % 6.2f, % 6.2f>' % tuple(array)

def pc(array):
    """Povray color syntax"""
    if type(array) == str:
        return 'color ' + array
    if type(array) == float:
        return 'rgb <%.2f>*3' % array
    if len(array) == 3:
        return 'rgb <%.2f, %.2f, %.2f>' % tuple(array)
    if len(array) == 4: # filter
        return 'rgbf <%.2f, %.2f, %.2f, %.2f>' % tuple(array)
    if len(array) == 5: # filter and transmit
        return 'rgbft <%.2f, %.2f, %.2f, %.2f, %.2f>' % tuple(array)

class POVRAY(EPS):
    default_settings = {
        # x, y is the image plane, z is *out* of the screen
        'display'      : True, # Display while rendering
        'pause'        : True, # Pause when done rendering (only if display)
        'transparent'  : True, # Transparent background
        'canvas_width' : None, # Width of canvas in pixels
        'canvas_height': None, # Height of canvas in pixels 
        'camera_dist'  : 50.,  # Distance from camera to front atom
        'image_plane'  : None, # Distance from front atom to image plane
        'camera_type'  : 'orthographic', # perspective, ultra_wide_angle
        'point_lights' : [],             # [[loc1, color1], [loc2, color2],...]
        'area_light'   : [(2., 3., 40.), # location
                          'White',       # color
                          .7, .7, 3, 3], # width, height, Nlamps_x, Nlamps_y
        'background'   : 'White',        # color
        'textures'     : None, # Length of atoms list of texture names
        }

    def __init__(self, atoms, scale=1.0, **parameters):
        for k, v in self.default_settings.items():
            setattr(self, k, parameters.pop(k, v))
        EPS.__init__(self, atoms, scale=scale, **parameters)

    def cell_to_lines(self, A):
        return npy.empty((0, 3)), None, None

    def write(self, filename, **settings):
        # Determine canvas width and height
        ratio = float(self.w) / self.h
        if self.canvas_width is None:
            if self.canvas_height is None:
                self.canvas_width = min(self.w * 15, 640)
            else:
                self.canvas_width = self.canvas_height * ratio
        elif self.canvas_height is not None:
            raise RuntimeError, "Can't set *both* width and height!"

        # Distance to image plane from camera
        if self.image_plane is None:
            if self.camera_type == 'orthographic':
                self.image_plane = 1 - self.camera_dist
            else:
                self.image_plane = 0
        self.image_plane += self.camera_dist

        # Produce the .ini file
        if filename.endswith('.pov'):
            ini = open(filename[:-4] + '.ini', 'w').write
        else:
            ini = open(filename + '.ini', 'w').write
        ini('Input_File_Name=%s\n' % filename)
        ini('Output_to_File=True\n')
        ini('Output_File_Type=N\n')
        ini('Output_Alpha=%s\n' % self.transparent)
        ini('; if you adjust Height, and width, you must preserve the ratio\n')
        ini('; Width / Height = %s\n' % repr(ratio))
        ini('Width=%s\n' % self.canvas_width)
        ini('Height=%s\n' % (self.canvas_width / ratio))
        ini('Antialias=True\n')
        ini('Antialias_Threshold=0.1\n')
        ini('Display=%s\n' % self.display)
        ini('Pause_When_Done=%s\n' % self.pause)
        ini('Verbose=False\n')
        del ini

        # Produce the .pov file
        w = open(filename, 'w').write
        w('#include "colors.inc"\n')
        w('#include "finish.inc"\n')
        w('\n')
        w('global_settings {assumed_gamma 1 max_trace_level 6}\n')
        w('background {%s}\n' % pc(self.background))
        w('camera {%s\n' % self.camera_type)
        w('  right -%.2f*x up %.2f*y\n' % (self.w, self.h))
        w('  direction %.2f*z\n' % self.image_plane)
        w('  location <0,0,%.2f> look_at <0,0,0>}\n' % self.camera_dist)
        for loc, rgb in self.point_lights:
            w('light_source {%s %s}\n' % (pa(loc), pc(rgb)))

        if self.area_light is not None:
            loc, color, width, height, nx, ny = self.area_light
            w('light_source {%s %s\n' % (pa(loc), pc(color)))
            w('  area_light <%.2f, 0, 0>, <0, %.2f, 0>, %i, %i\n' % (
                width, height, nx, ny))
            w('  adaptive 1 jitter}\n')

        w('\n')
        w('#declare simple = finish {phong 0.7}\n')
        w('#declare vmd = finish {'
          'ambient .0 '
          'diffuse .65 '
          'phong 0.1 '
          'phong_size 40. '
          'specular 0.500 }\n')
        w('#declare jmol = finish {'
          'ambient .2 '
          'diffuse .6 '
          'specular 1 '
          'roughness .001 '
          'metallic}\n')
        w('#declare ase2 = finish {'
          'ambient 0.05 '
          'brilliance 3 '
          'diffuse 0.6 '
          'metallic '
          'specular 0.70 '
          'roughness 0.04 '
          'reflection 0.15}\n')
        w('#declare ase3 = finish {'
          'ambient .15 '
          'brilliance 2 '
          'diffuse .6 '
          'metallic '
          'specular 1. '
          'roughness .001 '
          'reflection .0}\n')
        w('#declare glass = finish {' # Use filter 0.7
          'ambient .05 '
          'diffuse .3 '
          'specular 1. '
          'roughness .001}\n')
        w('\n')
        w('#macro atom(LOC, R, COL, FIN)\n')
        w('  sphere{LOC, R texture{pigment{COL} finish{FIN}}}\n')
        w('#end\n')
        w('\n')
        
        z0 = self.X[:, 2].max()
        self.X -= (self.w / 2, self.h / 2, z0)

        # Draw unit cell
        if self.C is not None:
            self.C -= (self.w / 2, self.h / 2, z0)
            self.C.shape = (2, 2, 2, 3)
            for c in range(3):
                for j in ([0, 0], [1, 0], [1, 1], [0, 1]):
                    w('cylinder {')
                    for i in range(2):
                        j.insert(c, i)
                        w(pa(self.C[tuple(j)]) + ', ')
                        del j[c]
                    w('%0.3f pigment {Black}}\n' % 0.05)

        # Draw atoms
        a = 0
        for loc, dia, color in zip(self.X, self.d, self.colors):
            tex = 'ase3'
            if self.textures is not None:
                tex = self.textures[a]
            w('atom(%s, %.2f, %s, %s) // #%i \n' % (
                pa(loc), dia / 2., pc(color), tex, a))
            a += 1


def write_pov(filename, atoms, run_povray=False, **parameters):
    if isinstance(atoms, list):
        assert len(atoms) == 1
        atoms = atoms[0]
    assert 'scale' not in parameters
    POVRAY(atoms, **parameters).write(filename)
    if run_povray:
        errcode = os.system('povray %s.ini 2> /dev/null' % filename[:-4])
        if errcode != 0:
            raise OSError('Povray failed with error code %d' % errcode)
