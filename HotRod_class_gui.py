# -*- coding: utf-8 -*-
"""
Created on Mon Mar 13 10:21:22 2017

@author: Saskia

Hot Rod class
"""
import pandas as pd
import numpy as np
from datetime import datetime
import math
from scipy.optimize import minimize
import time
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.artist
#radius_inner, radius_outer = 0.028, 0.047
#accuracy, margin = 0.06, 1
#td = TempData("templog_R2.csv",'hotrod_sensor_array.dat', radius_inner, radius_outer, accuracy, margin)

class HotRodException(Exception):
    pass

class HotRod:
    def __init__(self, temp_fn, sensor_fn, 
                 names = [],
                 initial = [], 
                 width = [],
                 uniform = [],
                 paramChoice = 1,
                 calChoice = "TNC",
                 radius_inner = 0.028, # m
                 radius_outer = 0.047, # m
                 acuracy = 0.06, # degrees celcius
                 margin = 1, # only used if the filter_data = Yes
                 pc = 2.5e6,
                 pcw = 4.1e6,
                 filter_data = 'No',
                 watt_threshold = 50., # wattage threshold, below this the heating element is off.
                 relay_dict = {1: -0.065, 2: -0.140, 3: -0.215},
                 ID = None):
        
        self.names = names
        self.initial = initial
        self.width = width
        self.uniform = uniform
        self.paramChoice = paramChoice
        self.calChoice = calChoice
        self.ID = ID
        self.acuracy = acuracy
        self.pc = pc
        self.pcw = pcw
        self.watt_threshold = watt_threshold
        self.temp_fn = temp_fn
        
        self.df_temp, self.df_batt = self.temp_parse(temp_fn)
        self.dur = self.pulse_duration()
        self.relay_depth = self.relaydepth(relay_dict)
        s = self.sensor_parse(sensor_fn, self.relay_depth)
        self.df_sensors, self.xyz, self.temp_obs, self.temp_error = self.xyzarray(s, radius_inner, radius_outer)
        self.sensor_check()
        
        if filter_data != 'No':
            self.xyz, self.temp_obs= self.filtered_data(margin)
        #self.Mean, self.Unif, self.Log, self.width, self.pmin, self.pmax = self.get_param_stats(self)
        
        
    def temp_parse(self, filename):
        with open(filename, 'r') as f:
            lines = f.readlines()    
        hds_line = lines[1]
        header_str = hds_line.strip().split(',')
        ci = len(header_str)
        while header_str[ci-1] == "":
            ci -= 1
        hds = header_str[2:ci]
        fmts =['%d/%m/%y %H:%M:%S','%d/%m/%Y %H:%M:%S']
        date_time, temp, info = [], [], []
        for i in range(2, len(lines)):
            data = lines[i].split(',')
            dt = data[0] + ' ' + data[1]
            
            dum = -1
            j = -1
            while dum < 0:
                j += 1
                try:
                    dd = datetime.strptime(dt, fmts[j])
                    # exit the loop on success
                    dum = 0
                except:
                    # repeat the loop on failure
                    pass  
            date_time.append(dd)

            dd = [float(k) for k in data[2:]]
            temp.append(dd[:-5])
            info.append(dd[-5:])
        temparray = np.asarray(temp)    
        infoarray = np.asarray(info) 
        
        df1 = pd.DataFrame(data = temparray, columns = hds)   
        dt = pd.to_datetime(date_time)
        df1.index = dt  
        df1a = df1.drop(df1.index[[0]])
        self.TT = df1a.values
        
        
        hds2 = ['Batt(V)','Pow','Relay','On Time','Heat Cur']
        df2 = pd.DataFrame(data = infoarray, columns = hds2)   
        df2.index = dt  
        df2['watts_calc'] = df2['Heat Cur'] * 12
        self.watts = df2['watts_calc'].values
 
        return(df1, df2)
        
    def pulse_duration(self):
        timestamp = self.df_temp.index.to_series().diff()
        dd = timestamp.dt.total_seconds()
        mean, mn, mx = dd.mean(), dd.min(), dd.max()
        if (mean == mn) and (mn == mx):
            self.timestep = mean
        else:
            raise TempDataException("Variable timesteps in temperature data: range %i - %i seconds"%(mn, mx))
            
        self.t = np.arange(1, np.shape(self.TT)[0] + 1, 1) * self.timestep
        heat_on = self.df_batt[self.df_batt['watts_calc'] > self.watt_threshold]
        return(len(heat_on) * self.timestep)
        
        
    def sensor_parse(self, sensor_fn, relay_depth):
        sensors = pd.read_csv(sensor_fn, sep = '\t',index_col = 0, header = 0)
        sensors['T sensor depth cor (mm)'] = sensors['T sensor depth (mm)'] * -1
        sensors['z'] =  sensors['T sensor depth cor (mm)']/1000 - relay_depth # subtract sensor depth
        sensors['x'] = np.nan
        sensors['y'] = np.nan

        return sensors
        
    def relaydepth(self,relay_dict):
        if 'Relay' in self.df_batt:
            nn = self.df_batt.Relay.unique()
            if len(nn) > 1:
                raise TempDataException('Multiple heat sources')
            else:
                heat_source = int(nn[0])
        else:
            raise TempDataException("Relay column not in dataframe")
        if heat_source not in relay_dict:
            raise TempDataException("Invalid heat source: {}".format(heat_source))
        relay_depth = relay_dict[heat_source]

        return relay_depth
        
    def sensor_check(self):
        temp_sensors = set(self.df_temp.columns)
        sens_sensors = set(self.df_sensors.index)
        if temp_sensors != sens_sensors:
            diff1 = sens_sensors - temp_sensors
            diff2 = temp_sensors - sens_sensors
            raise HotRodException("Sensor Mismatch: {} not in temperature sensors "\
                "and {} not in sensors data".format(diff1, diff2))
    
    
    def xyzarray(self, sensors, radius_inner, radius_outer):
        angle = np.arange(2*math.pi + math.pi/2, math.pi/2, -math.pi/4) 
        radii = [radius_outer, radius_inner] * int(len(angle)/2)
        
        x = np.cos(angle) * radii   
        y = np.sin(angle) * radii

        rod_id = ['T1','T2','T3','T4','T5','T6','T7','T8']
        coords = pd.DataFrame({'rod': rod_id, 'x': x, 'y': y})
        coords = coords.set_index('rod')
              
        for idx, row in sensors.iterrows():
            var = row['Sensor position ID'][:2]
            sensors.loc[idx, 'x'] = coords.loc[var,'x']
            sensors.loc[idx, 'y'] = coords.loc[var,'y']
        
        X = sensors[['x','y','z']].values 
        X[abs(X) < 1e-5] = 0

        self.sensor_id = sensors['Sensor position ID'].tolist()

        sensor_order = list(sensors.index.values) 
        temp_ordered = self.df_temp[sensor_order]        
        temp_ob = temp_ordered.values
        temp_obs = temp_ob[1:,:]
        temp_error = np.ones(np.shape(temp_obs)) * self.acuracy

        direction_bounds = []
        for i in range(3):
            var_mn = X[:,i].min()
            var_mx = X[:,i].max()
            direction_bounds.append([var_mn, var_mx])
        self.direction_bounds = np.asarray(direction_bounds)
        
        return sensors, X.T, temp_obs, temp_error

    def filtered_data(self, margin):
        temp_stats = self.df_temp.describe()
        temp_stats.ix['range', :] = temp_stats.ix['max', :] # - temp_stats.ix['min', :] 
        temp_summary = temp_stats.transpose()
        temp_resp = temp_summary[temp_summary['range'] >= (self.acuracy * self.margin)]
        sensors_filtered = list(temp_resp.index.values)  
                
        temp_filtered = self.df_temp[sensors_filtered] # select the data for the filtered sensors
        sensors_filtered = self.df_sensors.loc[sensors_filtered]
        
        filtered_loc = sensors_filtered[['x','y','z']].values
        filtered_loc[abs(filtered_loc) < 1e-5] = 0
        filtered_loc = filtered_loc.T
        
        self.sensor_id = sensors_filtered['Sensor position ID'].tolist()
        
        filtered_temp = temp_filtered.iloc[1:,:].values
        filtered_temp = filtered_temp.T # transpose to generate the correct cols/row orientation
        
        temp_error = np.ones(np.shape(filtered_temp)) * self.acuracy
                
        return filtered_loc, filtered_temp, temp_error


    def get_param_stats(self):
       
        self.Mean = np.asarray(self.initial)
        self.par_names = self.names

        boolean_unif = [True if d == "Uniform" else False for d in self.uniform]

        
        self.Unif = np.asarray(boolean_unif)

        self.pmin = []
        self.pmax = []        
        for i in range(len(boolean_unif)):
            if boolean_unif[i] == True:
                mn = self.Mean[i] - self.width[i]
                mx = self.Mean[i] + self.width[i]
                self.pmin.append(mn)
                self.pmax.append(mx)                
            else:
                mn = self.Mean[i] - (3*self.width[i])
                mx = self.Mean[i] + (3*self.width[i])
                self.pmin.append(mn)
                self.pmax.append(mx) 

        self.pmin[2:4] = [-1*np.pi] * 2#[-2*np.pi] * 2   
        self.pmax[2:4] = [np.pi] * 2 
           
        mn_diff =  2*np.pi / 2
        print(self.Mean)
        self.Mean[2:4] = [0.0, 0.0] 
        self.width[2:4] = 2*[np.pi]                 

        
    def max_responder(self):
        # identify the sensor with the largest response
        max_idx = np.argwhere(self.temp_obs == np.max(self.temp_obs))
        peak = max_idx[0]

        max_idx = peak[0]
        max_sen = peak[1]
        t_adj = self.t.reshape((1, len(self.t)))

        self.t_mu = np.dot(t_adj, self.temp_obs[:,max_sen]) / np.sum(self.temp_obs[:, max_sen])
        self.XYZ_max = self.xyz[:, max_sen]
        self.dist_mx = self.distance()

    
    def distance(self):
        # distance of the maximum response sensor to the heating element
        return((self.XYZ_max[0]**2 + self.XYZ_max[1]**2 + self.XYZ_max[2]**2)**0.5)
        
    
    def flow_magnitude(self):
        # Magnitude of the flow
        self.q_mag = (self.dist_mx / self.t_mu) * (self.pc/self.pcw)        
        
    def q_initial(self):
        self.max_responder()  
        self.flow_magnitude()
        self.theta = np.arctan(self.XYZ_max[1]/self.XYZ_max[0]) 
        self.phi = -np.arctan(self.XYZ_max[2])#/self.XYZ_max[0])  
        
    def angle2flux(self):
        qz = np.sin(self.phi) * self.q_mag
        qy = np.cos(self.phi) * np.sin(self.theta) * self.q_mag
        qx = np.cos(self.phi) * np.cos(self.theta) * self.q_mag     
        vector = np.asarray((qx,qy,qz))
        self.flux_xyz = np.nan_to_num(vector)
                                          
    def rotate(self):
        theta1 = self.theta 
        # Eq 8
        J = np.array([[np.cos(theta1),np.sin(theta1),0.],
                 [-np.sin(theta1),np.cos(theta1),0.],
                  [0.,0.,1.]])

        x_ = np.dot(J, self.xyz)
        
        phi = self.phi
        J_ = np.array([[np.cos(phi),0.,np.sin(phi)],
                 [0.,1.,0.],
                  [-np.sin(phi),0., np.cos(phi)]])
        _x_ = np.dot(J_,x_)
        
        qq = np.array([self.q_mag, 0, 0])
        return(_x_, qq)
    
    
    def temp_pulse_disp(self, q, U):
        qx = q[0]
        return(T)

    def temp_pulse(self, q, U):
        qx = q[0]

        if self.paramChoice == 1 :
            k2 = self.pcw/self.pc
            '''
            Eq 2: Dispersivity
            ''' 
            DispL = self.dispL
            DispT = self.dispT # Dz == Dy
            print(k2, DispL, DispT, self.pc)
            ns = np.shape(U)[1]
            nt = np.size(self.t)
            T = np.zeros((nt,ns))
#            print(U)
            for i in range(ns):
                '''
                Eq 15
                
                '''
#                print("here")
                x, y, z = U[0,i], U[1,i], U[2,i]
                front = (self.timestep/(8. * self.pc * (DispL**0.5) * DispT * (math.pi*self.t)**1.5))
#                print(front[0], DispL, DispT, self.pc)
#                print (i,x, y, z,self.timestep )
                back = np.exp(-((x-k2*qx*self.t)**2.)/(4.*DispL*self.t)-\
                              ((y**2.+z**2)/(4.*DispT*self.t))) 
#                print(back[0])
                T[:,i] = front * back
#                print(lkml)
        else:
            k1 = 4.*self.k/self.pc
            k2 = self.pcw/self.pc
            ns = np.shape(U)[1]
            nt = np.size(self.t)
            T = np.zeros((nt,ns))
            for i in range(ns):
                '''
                Eq 15 Assumes uniform dispersion
                '''
                x, y, z = U[0,i], U[1,i], U[2,i]
                front = (self.timestep/(self.pc * (4*math.pi*k1*self.t)**1.5))
                
                back = np.exp(-((x-k2*qx*self.t)**2.)/(4.*k1*self.t)-\
                              ((y**2.+z**2)/(4.*k1*self.t))) 
                
                T[:,i] = front * back
        
        return(T)
        

        
    def calculation(self):
        U, qq = self.rotate()
        T = self.temp_pulse(qq, U)
        npt = int(self.dur/(self.t[1]-self.t[0]))
        nT = np.shape(T)[0]
        Tp = np.copy(T)
        T = T*self.watts[0]
        '''
        Wattage is not constant, if it where then temp_pulse would have 
        watts*timestep (Eq 5). Given that it is not constant these 2 variables 
        (watts and timestep) have been seperated 
        '''
        for i in range(1,npt,1):
            T[i:,:] = T[i:,:] + self.watts[i] * Tp[:nT-i,:]
        return(T) 
        
    def launch(self, cal_variables):       
        # asign parameters
        if self.paramChoice == 1:
            self.q_mag_est = 10**cal_variables[0]
            self.pc = cal_variables[1]
            self.theta = cal_variables[2]
            self.phi = cal_variables[3]
            self.dispL = cal_variables[4]
            self.dispT = cal_variables[5]
        else:
            self.q_mag_est = 10**cal_variables[0]
            self.pc = cal_variables[1]
            self.theta = cal_variables[2]
            self.phi = cal_variables[3]
            self.k = cal_variables[4]
            
        T = self.calculation()
        obj = 0
        dum = (100 * (self.temp_obs - T))**2.
        obj = dum.sum()
#        print("Objective function value: ",obj)
        return(obj)

    def calibrated_run(self, cal_variables):       
        # asign parameters
        self.q_mag = 10**cal_variables[0]
        self.theta = cal_variables[1]
        self.phi = cal_variables[2]
        self.pc = cal_variables[3]
        self.dispL = cal_variables[4]
        self.dispT = cal_variables[5]

        T = self.calculation()
        return(T)        
        
  

    def set_propositions(self, proposed_val):
        proposed_val1 = np.copy(proposed_val)
        self.q_mag = 10**proposed_val1[0]

        self.theta = proposed_val[1]
        self.phi = proposed_val[2]

        self.pc = proposed_val1[3]

        self.dispL = proposed_val1[4]
        self.dispT =proposed_val1[5]


    def get_set_solve(self):
        self.temp_sim = self.calculation()
#        print("Solver End = ", time.strftime("%Y-%m-%d %H:%M:%S"))
    
    def solve(self, names):
        self.varNames = names
        self.q_initial()
        if self.paramChoice == 1:
            cal_variables = [np.log10(self.q_mag[0]),  
                            self.Mean[1],
                            self.theta, 
                            self.phi,
                            self.Mean[4],
                            self.Mean[5]]
                                    
        else:
            cal_variables = [np.log10(self.q_mag[0]),  
                            self.Mean[1],
                            self.theta, 
                            self.phi,
                            self.Mean[4]] 
        print(names)
        print(cal_variables)
#        print(self.Mean)
        self.bounds = tuple(zip(self.pmin, self.pmax))
        print(self.bounds)

        solution = minimize(self.launch, cal_variables, 
                        method = 'TNC', 
                        bounds = self.bounds, 
                        options = {"maxiter":5000})
        
        self.solved_output(solution)
        self.calibrated_sim = self.calibrated_run(solution.x)
        

    def solved_output(self, solution):
#        print(" ************ ")
        print(" ************ ")
        print("Minimize function Output")
        print("------------------------")
        print(solution)
        print(" ************ ")
        print(" ************ ")
        print("Calibration Output Values: ")
        self.calibrated_variables = solution.x
        for i in range(len(self.calibrated_variables)):
            if self.varNames[i] == 'ln(q)':
                print('Q_mag = ', 10 ** self.calibrated_variables[i])
            else:
                print(self.varNames[i],' = ', self.calibrated_variables[i])
        
        print(" ************ ")
        print(" ************ ") 
        print(" Flux: x, y, z ") 
        self.angle2flux()
        print(self.flux_xyz)
        print(" ************ ") 
        
#        fn = self.temp_fn[:-4] + "_calibrated.dat"
#        with open(fn, 'w') as f:
#            f.write("%s\t%s\t%s\t%s\n"%('var','cal','bound_min','bound_max'))
#            for i in range(len(self.calibrated_variables)):
#                f.write("%s\t%g\t%g\t%g\n"%(self.varNames[i],
#                                            self.calibrated_variables[i],
#                                            self.bounds[i][0], 
#                                            self.bounds[i][1]))

                
                
    def figures_cal(self, save_figures):      
          
        temp_max = self.calibrated_sim.max() + self.acuracy
        sets = len(self.df_sensors)/9.
        sets = math.ceil(sets)
        print(sets)
        count = 0 
        for j in range(sets):
            
            fig = plt.figure(figsize=(10, 10),facecolor='w', edgecolor='k')
            ss = 1
            
            if (count + 9) <= len(self.df_sensors):
                count_up = count + 9
            else:
                count_up = len(self.df_sensors)
                
            for i in range(count, count_up):
                ax = plt.subplot(3, 3, ss)
                
                ax.plot(td.t/60 , self.temp_obs[:, i])
                ax.plot(td.t/60 , self.calibrated_sim[:, i])
                
                ax.set_title("%s"%(self.sensor_id[i]))
                ax.set_ylim(-0.1,temp_max)
                ss += 1
                count += 1
                
            plt.tight_layout()
            
            if save_figures == "Yes":
                figname = "fig_%s_%i.png"%(file[:-4], j)
                plt.savefig(figname, dpi = 400)
                
            plt.show()     

       
    def figure_3d(self):
        
        self.tr1 = self.xyz[:,:7]    
        self.dt = np.arange(0, len(self.t))
         
        mag = np.sqrt(self.flux_xyz[0] ** 2 + self.flux_xyz[1] ** 2. + self.flux_xyz[2] ** 2.)    
        # scale the flow vector so that its length equals 0.05m (i.e., so that its visable)
        a = 0.05/mag
        self.ini_time = 60
        self.max_responder()    
        
        self.fig = plt.figure(figsize=(10, 10))
        ax = self.fig.add_subplot(111,projection='3d')
        self.fig.subplots_adjust(bottom=0.25)      
        
        # location of rod 1
        ax.plot(self.tr1[0,[0,6]], self.tr1[1,[0,6]], self.tr1[2,[0,6]] + self.relay_depth,
                'k-', 
                zdir = 'z',
                label = 'T1 sensors')
        # sensor locations
        ax.scatter(self.xyz[0,:], self.xyz[1,:], self.xyz[2,:] + self.relay_depth, 
                   zdir = 'z',  
                   s = 10, 
                   c = 'r', 
                   edgecolor = '')
        # heat source location
        ax.scatter(0, 0, self.relay_depth, 
                   zdir = 'z', 
                   s = 50, 
                   c = 'orange', 
                   edgecolor = '', 
                   label = 'heat source')
        # initial conditions
        ax.scatter([self.XYZ_max[0]], [self.XYZ_max[1]], [self.relay_depth + self.XYZ_max[2]],
                zdir = 'z', 
                c = 'b', 
                s = 100, 
                label = 'max responder')
        # calibrated vector
        ax.plot([0, self.flux_xyz[0]*a], [0, self.flux_xyz[1]*a], [self.relay_depth, self.relay_depth + self.flux_xyz[2]*a],
                zdir = 'z', 
                c = 'k',
                label = 'flow')
        # observed temperature
        self.sc = ax.scatter(self.xyz[0,:], self.xyz[1,:], self.xyz[2,:] + self.relay_depth, 
                             zdir = 'z',  
                             s = 2000 * self.temp_obs[self.ini_time, :],
                            c = 'red', 
                            edgecolor = '', 
                            label = 'sensors')
    
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        
        pos, name, min, max = 0.03, ' Time \n(x*%i sec)'%(self.timestep), 0, len(self.t)        
        time_slider_ax = plt.axes([0.1, (2*pos), 0.8, 0.03], axisbg='lightgoldenrodyellow')
        time_slider = Slider(time_slider_ax, name, min, max, valinit = self.ini_time, valfmt = '%0.0f')
   
        def update(val):
            new_val = int(time_slider.val)
            new_size_array = np.ravel(2000 * (self.temp_obs[new_val, :]))
            matplotlib.artist.setp(self.sc, sizes = new_size_array)
            self.fig.canvas.draw_idle()
            
        time_slider.on_changed(update)  
        plt.show()
 
# ***************************************************************************
if __name__ == '__main__':
    
    td = HotRod("templog_H_R2.csv",'hotrod_sensor_array_updated.dat')
    td.get_param_stats()
    td.solve()
    save_figures = "No"
    td.figures_cal(save_figures)
    td.figure_3d()

