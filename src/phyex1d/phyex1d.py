"""
Phyex1d is the main object of the phyex1d package; it deals with the execution
"""

import netCDF4
from pppy import PPPYComp
from .physics import PhysicsAromeTQ, PhysicsAromeThetaR
from .grid import Grid


class Phyex1d(PPPYComp):
    """
    Deals with the execution of the 1D model
    """
    def __init__(self, inputfile, experiments, output_dir, comp_name=''):
        """
        :param inputfile: netcdf file describing the case
        :param experiments: list of dictionaries describing the experiments
                            allowed keys are:
                              - grid: Grid name or file name containing a grid description
                                      (with '.grid' extension).
                              - dt: Timestep (s)
                              - name: experiment name
        :param timestep: timestep legth (s)
        :param comp_name: comparison name
        """
        def string2filename(s):
            return s.replace('"', '').replace("'", "").replace(' ', '_')

        self.inputfile = inputfile
        schemes = []
        for exp in experiments:
            cls = {'PhysicsAromeTQ': PhysicsAromeTQ,
                   'PhysicsAromeThetaR': PhysicsAromeThetaR,
                  }[exp['class'] if 'class' in exp else 'PhysicsAromeTQ']
            schemes.append(cls(dt=float(exp.get('dt', 1.)),
                               method='step-by-step',
                               name=exp.get('name', 'Experiment'),
                               tag=string2filename(exp.get('name', 'Experiment')),
                               inputfile=inputfile,
                               grid=Grid(exp.get('grid', 'L90arome')),
                               pyphyex=exp.get('pyphyex', None),
                               namel=exp.get('namel', 'default'),
                               dx=exp.get('dx', 1300),
                               dy=exp.get('dy', 1300)))
        super().__init__(schemes=schemes,
                         output_dir=output_dir,
                         duration=self.get_duration(),
                         init_state={},
                         name=comp_name,
                         tag=string2filename(comp_name))

    def get_duration(self):
        """
        Read the simulation duration in the netcdf file describing the case
        """
        with netCDF4.Dataset(self.inputfile, 'r') as nc:
            duration = nc['time'][-1] - nc['time'][0]

        return duration
