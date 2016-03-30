# -*- coding: utf-8 -*-
"""
Created on Fri Jul 24 18:36:06 2015

@author: Pedro Leal
"""
import pickle
import time
import os, sys
import matplotlib.pyplot as plt

import xfoil_module as xf

#adding path to static model
lib_path = os.path.abspath(os.path.join('..'))
sys.path.append(lib_path)

import static_model as model 

in_Abaqus = False

class DOE:
    """Create a Design of Experiences Environment."""
    def __init__(self, levels=5, driver='Taguchi'):
        self.levels = levels
        self.driver = driver
        # All variable will be defined through the add_variable method
        self.variables = []
        # For the influence method, we need a list of the names of all variables
        self.variables_names = []

    def add_variable(self, name, lower, upper, levels=None, type=float):
        """Add variables to the DOE problem. """
        if levels is None:
            levels = self.levels

        try:
            self.variables.append({'upper': upper, 'lower': lower, 'name': name,
                                   'levels': levels, 'type': type})
            self.variables_names.append(name)
        except:
            print 'Forgot to define upper, lower or name'

    def define_points(self):
        """
        Method to define the points to be evaluated based on the results from the
        distribution given by the array method and the bound defined by the
        add_variable method"""
        self.n_var = 0
        self.n_var_2 = 0

        for variable in self.variables:
            if variable['levels'] == self.levels:
                self.n_var += 1
            elif variable['levels'] == 2:
                self.n_var_2 += 1
            else:
                raise Exception('A variable has a number of levels that is ' +
                                'not the default or 2')
        if self.driver == 'Taguchi':
            self.Taguchi()
        elif self.driver == 'Full Factorial':
            self.FullFactorial()
        #TODO: Define a FullFactorial Driver or something like that
        self.domain = {}

        for j in range(self.n_var+self.n_var_2):
            upper = self.variables[j]['upper']
            lower = self.variables[j]['lower']
            levels = self.variables[j]['levels']
            type = self.variables[j]['type']
            
            dummy = []
            for i in range(self.runs):
                scale = self.array[i][j]
                if type == int and (scale*(upper-lower) % (levels-1.) != 0):
                    raise Exception('The bounds of the defined integer are not '+
                                    'compatible with the number of levels.')
                else:
                    dummy.append(lower + scale*(upper-lower) / (levels-1.))
            self.domain[self.variables[j]['name']] = dummy
        
    def run(self, function, inputs = None, cte_input=None,
            dependent_variables = None, parameters = None):
        """Runs and saves the results for the configurations obtained in define_points
        method.
        - inputs: if defined will define the format with which the inputs 
                  will be passed to the main function.
        - cte_input : if defined, is a dictionary containing the constant
          inputs.
        - dependent_variables: if defined, it creates a relationship between
          different variables such as {'t_spar':'t_rib'}
        - parameters: some main functions need extra parameters. this is a list
        """
        
        DataFile = open('all_data.txt','w')
        key_list = ['xs-', 'ys-', 'xs+', 'ys+', 'xl-', 'yl-', 'xl+', 'yl+']
        output_list = ['delta_xi', 'theta', 'T', 'k']
        for key in ['i'] + key_list + output_list:
            DataFile.write(key + '\t')
        DataFile.write('\n')
        DataFile.close()
        
        def set_input(self, run):
            output = {}
            for key in self.domain:
                output[key] = self.domain[key][run]
            return output
        
        for i in range(self.runs):
            input = set_input(self,i)
            # If there is a constant input, add it to input dictionary
            if cte_input != None:
                input.update(cte_input)
            if dependent_variables != None:
                for key_dependent in dependent_variables:
                    key_independent = dependent_variables[key_dependent]
                    input.update({key_dependent : input[key_independent]})
            
            DataFile = open('all_data.txt','a')
            DataFile.write(str(i))
            for key in key_list:
                DataFile.write( '\t' + str(input[key]))
            DataFile.close()
            
            if inputs != None:
                new_input = {}
                for key_input in inputs:
                    if type(inputs[key_input]) == dict:
                        new_input[key_input] = {}
                        for key_var in inputs[key_input]:
                            new_input[key_input][key_var] = input[inputs[key_input][key_var]]
                    else:
                       new_input[key_input]= input[inputs[key_input]] 
                input = new_input
            print input

            if parameters == None:
                result = function(input)
            else:
                result = function(input, parameters = parameters)
            
            if i == 0:
                # We will save the name of the putputs for plotting and etc
                self.output_names = [key for key in result]
                self.output = {}
                for key in self.output_names:
                    self.output[key] = []
            for key in self.output_names:
                self.output[key].append(result[key])
            
            DataFile = open('all_data.txt','a')
            for output in output_list:
                DataFile.write( '\t' + str(self.output[output][-1]))    
            DataFile.write('\n')
            DataFile.close()
    def find_influences(self, not_zero=False):
        """ Calculate average influence of each variable over the
        objective functions. If refinement_criteria is defined, certain points
        will be eliminated. Works for Taguchi, Full Factorial and probably
        anything else.
        """
        self.influences = {}
        
        # For refinement reasons during plotting, it is relevant to
        # know which ones have zeros
        self.equal_to_zero = {key:[False]*(self.n_var+self.n_var_2)*self.levels for
                         key in self.output_names}
        for output_name in self.output_names:
            Y = self.output[output_name]
            # List of outputs
            self.influences[output_name] = []
            
            for var in self.variables_names:
                X = self.domain[var]
                # Eliminate repetitions by transforming the list in to a set and 
                # sort them. X_set will be used for counting
                unique_X = sorted(set(X))
                # For each unique X, the average of its values will be calculated
                
                for j in range(len(unique_X)):
                    indices = [i for i, x in enumerate(X) if x == unique_X[j]]
                    # Filter option
                    if not_zero:
                        # Evaluate if any of the values of output is
                        # zero
                        for key, item in self.output.items():
                            if 0 in item:
                                # Eliminate it from the list of indices
                                for i in indices:
                                    if self.output[key][i] == 0:
                                        indices.remove(i)
                                        self.equal_to_zero[key][j] = True
                    # Count number of times the variable repeats itself
                    count = len(indices)
                    # Create an empyt slot in Y_DOE list to add all Ys
                    self.influences[output_name].append(0)
                    # Average of all points with same X value
                    dummy = 0
                    for index in indices:
                        dummy += Y[index]/count
                        # Add to the last term of Y_DOE (sum of all)
                        self.influences[output_name][-1] += Y[index]/count

    if not in_Abaqus:
        def plot(self, shadow = [], xlabel = None, ylabel = None):
            import matplotlib.pyplot as plt
            
            def list_to_string(self, separator=', '):
                """Summ all the elements of a list of strings in to a string"""
                resultant_string = ''
                for component in self.variables_names:
                    resultant_string += component + separator
                # Remove the last separator.
                resultant_string = resultant_string[:-len(separator)]
                return resultant_string
                
            def create_ticks(self):
                # In japanese mora is the length of each sylab, here it is the length of e
                if self.levels == 2:
                    mora = ['-', '+']
                elif self.levels == 3:
                    mora = ['-', 'o', '+']
                elif self.levels == 4:
                    mora = ['-', '-o', 'o+', '+']
                elif self.levels == 5:
                    mora = ['-', '-o', 'o', 'o+', '+']
                else:
                    raise Exception('n_range to high, max is 5!')
                
                # Replicate standard for all variables
                return (self.n_var_2)*['-', '+'] + (self.n_var)*mora
            
            def subtick_distance(self, border_spacing):
                """Function togenerate the distances of the second x axis
                using figtext"""
                
                # normalizing values forimage to be >0 and <1 
                norm = (2*border_spacing + self.levels*self.n_var - 1)

                # Initial proportional distance
                x0 = border_spacing/norm
                distances = []
                for i in range(len(self.variables_names)):
                    current = x0 + i*(self.levels - 1)/norm
                    print 'x0', current
                    if self.levels % 2 == 0: # if even
                        if i==0:
                            current += (self.levels)/2./norm
                        else:
                            current += (self.levels + 1)/2./norm
                    else: # if odd
                        if i == 0:
                            current += (self.levels/2 )/norm
                        else:
                            current += (self.levels/2 +1)/norm
                    print current
                    distances.append(current)
                return distances
            
            # IF the user wants to add pretty names, if not just do with the
            # variable names
            if xlabel == None:
                xlabel = self.variables_names
            if ylabel == None:
                ylabel = self.output_names
            ticks = create_ticks(self)
            border_spacing = .2
            for output in self.output_names:
                Y = self.influences[output]
                plt.subplot(100*len(self.output_names) + 11 + 
                            self.output_names.index(output))
         
                # Creating dummy values for horizontal axis
                xi = range((self.n_var+self.n_var_2) * self.levels)
                # Define ticks
                plt.xticks(xi, ticks)
#                plt.fill_between(xi, min(Y) - 0.05*(max(Y)-min(Y)),
#                                 max(Y) + 0.05*(max(Y)-min(Y)),
#                                 where = self.equal_to_zero[output],color = '0.75')
                for i in range(self.n_var+self.n_var_2):
                    plt.plot(xi[i*self.levels : (i+1) * self.levels],
                             Y[i*self.levels : (i+1) * self.levels],
                             '-o')
    #                if output in shadow:
    #                    plt.plot(xi[i*self.levels : (i+1) * self.levels],
    #                             Y[i*self.levels : (i+1) * self.levels],
    #                             '--',color=plt.getp(line, 'linewidth'))
                plt.ylabel(ylabel[self.output_names.index(output)])
#                plt.xlabel("Design Variables ("+list_to_string(self)+")")
            
                if self.output_names.index(output) == 0:
                    plt.title("Design of Experiment: %i level %s" %
                              (self.levels, self.driver))

                plt.xlim([-border_spacing, max(xi) + border_spacing])
                plt.ylim(min(Y) - 0.05*(max(Y)-min(Y)), 
                         max(Y) + 0.05*(max(Y)-min(Y)))
                plt.grid()
                plt.show()

            # Create the second x axis
            distances = subtick_distance(self, border_spacing)
            print distances
            for i in range(len(distances)):
                plt.annotate(xlabel[i], xy =(distances[i], 0),
                             xytext = (0, -25), xycoords='axes fraction',
                             textcoords='offset points', horizontalalignment='center',
                             verticalalignment='center')

        def plot_domain(self, Xaxis, Yaxis):
            """Plots all the points in a 2D plot for the definided Xaxis and 
            Yaxis
            
            param: Xaxis: string containing key for x axis
            param: Yaxis: string containing key for y axis"""
            
            plt.scatter(self.output[Xaxis],self.output[Yaxis])
            plt.xlabel(Xaxis)
            plt.ylabel(Yaxis)
            
        def load(self, filename, variables_names, outputs_names, header=None):
            """ Load data from text file with results of DOE.
                TODO: NEED TO INCLUDE POSSIBILITY FOR TWO LEVEL VARIABLE
                - input:
                    - header: If not specified, the first line in the file
                      will be considered to be the header.
            """
            if self.variables_names == []:
                if header == None:
                    Data = xf.output_reader(filename=filename)
                else:
                    Data = xf.output_reader(filename=filename, header=header)
                if True==True:
                    self.output_names = outputs_names
                    self.variables_names = variables_names
                    self.n_var = len(variables_names)
                    self.n_var_2 = 0
                    self.output = {key:Data[key] for key in self.output_names}
                    self.domain = {key:Data[key] for key in self.variables_names}
    #            except:
    #                raise Exception('Something wrong with variables_names and '+ 
    #                                'outputs_names.')
            else:
                raise Exception('Cannot atribute variables and load data at the' +
                                'same object.')
#    def save(self):
        
    def Taguchi(self):
        """ Find the necessary Taguchi array."""
        self.runs = 50
        Taguchi_L50 = [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                       [0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                       [0, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
                       [0, 0, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3],
                       [0, 0, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4],
                       [0, 1, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4],
                       [0, 1, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0],
                       [0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1],
                       [0, 1, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2],
                       [0, 1, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3],
                       [0, 2, 0, 2, 4, 1, 3, 3, 0, 2, 4, 1],
                       [0, 2, 1, 3, 0, 2, 4, 4, 1, 3, 0, 2],
                       [0, 2, 2, 4, 1, 3, 0, 0, 2, 4, 1, 3],
                       [0, 2, 3, 0, 2, 4, 1, 1, 3, 0, 2, 4],
                       [0, 2, 4, 1, 3, 0, 2, 2, 4, 1, 3, 0],
                       [0, 3, 0, 3, 1, 4, 2, 4, 2, 0, 3, 1],
                       [0, 3, 1, 4, 2, 0, 3, 0, 3, 1, 4, 2],
                       [0, 3, 2, 0, 3, 1, 4, 1, 4, 2, 0, 3],
                       [0, 3, 3, 1, 4, 2, 0, 2, 0, 3, 1, 4],
                       [0, 3, 4, 2, 0, 3, 1, 3, 1, 4, 2, 0],
                       [0, 4, 0, 4, 3, 2, 1, 3, 2, 1, 0, 4],
                       [0, 4, 1, 0, 4, 3, 2, 4, 3, 2, 1, 0],
                       [0, 4, 2, 1, 0, 4, 3, 0, 4, 3, 2, 1],
                       [0, 4, 3, 2, 1, 0, 4, 1, 0, 4, 3, 2],
                       [0, 4, 4, 3, 2, 1, 0, 2, 1, 0, 4, 3],
                       [1, 0, 0, 0, 3, 4, 3, 2, 1, 4, 1, 2],
                       [1, 0, 1, 1, 4, 0, 4, 3, 2, 0, 2, 3],
                       [1, 0, 2, 2, 0, 1, 0, 4, 3, 1, 3, 4],
                       [1, 0, 3, 3, 1, 2, 1, 0, 4, 2, 4, 0],
                       [1, 0, 4, 4, 2, 3, 2, 1, 0, 3, 0, 1],
                       [1, 1, 0, 1, 0, 2, 2, 1, 3, 4, 4, 3],
                       [1, 1, 1, 2, 1, 3, 3, 2, 4, 0, 0, 4],
                       [1, 1, 2, 3, 2, 4, 4, 3, 0, 1, 1, 0],
                       [1, 1, 3, 4, 3, 0, 0, 4, 1, 2, 2, 1],
                       [1, 1, 4, 0, 4, 1, 1, 0, 2, 3, 3, 2],
                       [1, 2, 0, 2, 2, 0, 1, 4, 4, 3, 1, 3],
                       [1, 2, 1, 3, 3, 1, 2, 0, 0, 4, 2, 4],
                       [1, 2, 2, 4, 4, 2, 3, 1, 1, 0, 3, 0],
                       [1, 2, 3, 0, 0, 3, 4, 2, 2, 1, 4, 1],
                       [1, 2, 4, 1, 1, 4, 0, 3, 3, 2, 0, 2],
                       [1, 3, 0, 3, 4, 3, 0, 1, 4, 1, 2, 2],
                       [1, 3, 1, 4, 0, 4, 1, 2, 0, 2, 3, 3],
                       [1, 3, 2, 0, 1, 0, 2, 3, 1, 3, 4, 4],
                       [1, 3, 3, 1, 2, 1, 3, 4, 2, 4, 0, 0],
                       [1, 3, 4, 2, 3, 2, 4, 0, 3, 0, 1, 1],
                       [1, 4, 0, 4, 1, 1, 4, 2, 3, 3, 2, 0],
                       [1, 4, 1, 0, 2, 2, 0, 3, 4, 4, 3, 1],
                       [1, 4, 2, 1, 3, 3, 1, 4, 0, 0, 4, 2],
                       [1, 4, 3, 2, 4, 4, 2, 0, 1, 1, 0, 3],
                       [1, 4, 4, 3, 0, 0, 3, 1, 2, 2, 1, 4]
                      ]
        # Initialize the Taguchi array.
        self.array = self.runs*[[]]
        # The reange of easch array was defined in:
        # https://controls.engin.umich.edu/wiki/index.php/
        # Design_of_experiments_via_taguchi_methods:_orthogonal_arrays
        if (self.n_var >= 7 and self.n_var <= 12) and self.n_var_2 <= 1:
            # Since the first column is for two level variables, we can ignore it.
            for i in range(self.runs):
                self.array[i] = Taguchi_L50[i][1-self.n_var_2 : self.n_var+1]
                
    def FullFactorial(self):
        """Define array for Full Factorial for a given number of
        levels.
        """
        def product(*args, **kwds):
            """ Returns all the possible combinations beween two lists
            or between itself.
            
            >>> print product('ABCD', 'xy')
            >>> Ax Ay Bx By Cx Cy Dx Dy
            
            >>> print product(range(2), repeat=3)
            >>>000 001 010 011 100 101 110 111
            
            Source: itertools
            """
            pools = map(tuple, args) * kwds.get('repeat', 1)
            result = [[]]
            for pool in pools:
                result = [x+[y] for x in result for y in pool]
            for prod in result:
                yield tuple(prod)

        self.array = []

        possibilities = [i for i in range(self.levels)]

        for subset in product(possibilities, repeat = self.n_var):
            self.array.append(subset)

        self.runs = len(self.array)

    def find_nadir_utopic(self, not_zero=True):
        """Find the minimum and maximum, nadir and utopic, for each
        output variable.
        
        This function is quite relevant for normalizing the objective
        function for the optimization.
        
        param: not_zero: filters the zero values out.
        
        returns: attributes utopic and nadir dictionaries for the output
                 variables, each containing a float value
                 
        sources: http://stackoverflow.com/questions/16122362/python-matplotlib-how-to-put-text-in-the-corner-of-equal-aspect-figure
        """
        # First verify if there are any zeros, if true, get them out
        equal_to_zero = {}
        for key in self.output_names:
            equal_to_zero[key] = [False]*len(self.output[key])
        
            if not_zero:
                for i in range(len(self.output[key])):
                    for key2 in self.output_names:
                        if self.output[key2][i] == 0:
                            equal_to_zero[key][i] = True
                        elif equal_to_zero[key][i] != True:
                            equal_to_zero[key][i] = False            
        
        # Now we can find the nadir and the utopic points
        self.nadir = {}
        self.utopic = {}                 
        for key in self.output_names:
            self.nadir[key] = 0
            self.utopic[key] = 99999999999999999.
            
            for i in range(len(self.output[key])):
                if equal_to_zero[key][i] != True and self.output[key][i] < self.utopic[key]:
                    self.utopic[key] = self.output[key][i]

                if equal_to_zero[key][i] != True and self.output[key][i] > self.nadir[key]:
                    self.nadir[key] = self.output[key][i]
                    
        
if __name__ == "__main__":
    import matlab.engine
    
    #Start Matlab engine
    eng = matlab.engine.start_matlab()
    #Go to directory where matlab file is
    eng.cd('..')
    eng.cd('SMA_temperature_strain_driven')
    
    chord = 0.6175
    x_hinge = 0.25
    
    safety = 0.05*chord
    
    problem = DOE(levels=2, driver='Full Factorial')
    problem.add_variable('xs-', lower = x_hinge/2. , upper = x_hinge, type=float)
    problem.add_variable('ys-', lower = -.9, upper = -0., type=float)
    problem.add_variable('xs+', lower = x_hinge + safety, upper = chord - safety, type=float)
    problem.add_variable('ys+', lower = 0., upper = .9, type=float)
    problem.add_variable('xl-', lower = x_hinge/2., upper = x_hinge, type=float)
    problem.add_variable('yl-', lower = -.9, upper = 0.9, type=float)
    problem.add_variable('xl+', lower = x_hinge + safety, upper = chord - safety, type=float)
    problem.add_variable('yl+', lower = -.9, upper = 0.9, type=float)
    problem.define_points()
    
    #inputs [sma, linear, sigma_o]
    inputs = {'sma':{'x-':'xs-', 'y-':'ys-', 'x+':'xs+', 'y+':'ys+'},
              'linear':{'x-':'xl-', 'y-':'yl-', 'x+':'xl+', 'y+':'yl+'}}
    problem.run(model.run, inputs = inputs, parameters = [eng])
    
    timestr = time.strftime('%Y%m%d')
    fileObject = open('DOE_'+ problem.driver + '_' + timestr,'wb') 
#    fileObject = open('DOE_FullFactorial_20150828','wb') 
    pickle.dump(problem, fileObject)   
    fileObject.close()   
        
#    problem.load(filename='FullFactorial.txt', 
#                 variables_names= ['Al0', 'Al1'],
#                 outputs_names = ['Weight', 'Lift', 'Drag', 'MaxMises', 'DispTip', 'EigenValue', 'Velocity'])
    problem.find_influences(not_zero=True)
    problem.find_nadir_utopic(not_zero=True)
    print 'Nadir: ', problem.nadir
    print 'Utopic: ', problem.utopic
#    problem.plot(xlabel = ['$A_{l_0}$', '$A_{l_1}$'],
#                 ylabel = ['Weight(N)', 'Lift(N)', 'Drag(N)', 'MaxMises(Pa)',
#                           'Displacement(m)','Eigenvalue', 'Velocity(m/s)'])
#    print problem.influences
