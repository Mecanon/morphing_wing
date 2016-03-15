# -*- coding: utf-8 -*-
"""
Created on Tue Mar 15 14:53:27 2016

@author: Pedro Leal
"""
import sqlitedict
from pprint import pprint

from openmdao.api import IndepVarComp, Group, Problem, ScipyOptimizer, ExecComp, SqliteRecorder, Component
#from openmdao.test.paraboloid import Paraboloid
from openmdao.drivers.fullfactorial_driver import FullFactorialDriver

class Paraboloid(Component):
    """ Evaluates the equation f(x,y) = (x-3)^2 + xy + (y+4)^2 - 3 """

    def __init__(self):
        super(Paraboloid, self).__init__()

        self.add_param('x', val=0.0)
        self.add_param('y', val=0.0)

        self.add_output('f_xy', val=0.0)

    def solve_nonlinear(self, params, unknowns, resids):
        """f(x,y) = (x-3)^2 + xy + (y+4)^2 - 3
        """

        x = params['x']
        y = params['y']

        unknowns['f_xy'] = (x-3.0)**2 + x*y + (y+4.0)**2 - 3.0

    def linearize(self, params, unknowns, resids):
        """ Jacobian for our paraboloid."""

        x = params['x']
        y = params['y']
        J = {}

        J['f_xy', 'x'] = 2.0*x - 6.0 + y
        J['f_xy', 'y'] = 2.0*y + 8.0 + x
        return J
        
top = Problem()
root = top.root = Group()

root.add('p1', IndepVarComp('x', 50.0), promotes=['x'])
root.add('p2', IndepVarComp('y', 50.0), promotes=['y'])
root.add('comp', Paraboloid(), promotes=['x', 'y', 'f_xy'])

top.driver = FullFactorialDriver(num_levels=4)
top.driver.add_desvar('x', lower=-50.0, upper=50.0)
top.driver.add_desvar('y', lower=-50.0, upper=50.0)

top.driver.add_objective('f_xy')

recorder = SqliteRecorder('paraboloid')
recorder.options['record_params'] = True
recorder.options['record_unknowns'] = True
recorder.options['record_resids'] = False
recorder.options['record_metadata'] = False
top.driver.add_recorder(recorder)

top.setup()
top.run()
top.cleanup()

db = sqlitedict.SqliteDict( 'paraboloid', 'openmdao' )

data = db['rank0:Driver/1']

#Collect data inside recorded file
n = len(db.keys()) #number of recorded iterations

final_data = {'x':[], 'y':[], 'f_xy':[]}
for i in range(n):
    iteration_name = 'rank0:Driver/' + str(i)
    data = db[iteration_name]
    data = data['Unknowns']
    final_data['x'].append(data['x'])
    final_data['y'].append(data['y'])
    final_data['f_xy'].append(data['f_xy'])