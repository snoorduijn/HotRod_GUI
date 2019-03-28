# -*- coding: utf-8 -*-
"""
Created on Tue Mar 28 14:47:20 2017

@author: Saskia
"""

import multiprocessing as mp
import matplotlib
from matplotlib.backend_bases import key_press_handler
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.artist
from HotRod_class_gui import HotRod
from numpy import arange, sin, pi
import numpy as np
import pandas as pd
import tkinter.ttk as ttk
from tkinter import messagebox
import sys
from DREAM import DREAM

if sys.version_info[0] < 3:
    print("System Version: %i"%sys.version_info[0])
    import Tkinter as Tk

else:
    print("System Version: %i"%sys.version_info[0])
    import tkinter as Tk

print("matplotlib: ",matplotlib.__version__)
print("numpy: ",np.__version__)
print("pandas: ",pd.__version__)
print("tkinter.ttk: ",ttk.__version__)

class hotrod3Dfigure(Tk.Frame):
    
    def __init__(self, master, flds, value, lbls):
        Tk.Frame.__init__(self) # to initialise the class Tk.Frame
        self.master = master
        self.flds = flds
        self.lbls = lbls
        self.value = value
        
        self.figFrame = Tk.Frame(self.master)
        for row_index in range(10):
            self.figFrame.rowconfigure(row_index, weight=1)
            for col_index in range(10):
                self.figFrame.columnconfigure(col_index, weight=1)  
        
        self.figFrame.grid(row=3, column=0)
        
        self.figure_drawn = False
        self.setup_figure()
        
    def setup_figure(self):
#        self.f = Figure(figsize=(5, 4), dpi=100)
        self.f = Figure(dpi=100)
        # a tk.DrawingArea
        # the canvas need to be drawn first to enable rotation of the 3D plot
        self.canvas = FigureCanvasTkAgg(self.f, master=self.figFrame)
        self.canvas.get_tk_widget().configure(background='white')

#        self.canvas.get_tk_widget().grid(column = 0, row = 0)#(side=Tk.TOP, fill=Tk.BOTH, expand=1) 
        self.canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1) 
        
        toolbar = NavigationToolbar2Tk(self.canvas, self.figFrame)
        toolbar.update()

        self.slider_val0 = 1 # initial value
        self.slider_val = Tk.IntVar() # set val as an interger value
#        print(self.slider_val, self.slider_val0)
        self.time_slider = Tk.Scale(self.figFrame,
                               label = 'Time (sec*3)',
                               from_ = 1,
                               to = 600,
                               length = 600,
                               tickinterval = 50,
                               orient = 'horizontal',
                               command = self.updateScale)
        
        self.time_slider.set(self.slider_val0)

        
    def updateScale(self,val):
#        print ("scale is now %s" % (val))
        if not self.figure_drawn:
            return
        new_val = int(val)
        new_size_array = np.ravel(500 * self.td.temp_obs[new_val, :])
        matplotlib.artist.setp(self.sc, sizes = new_size_array)
        self.f.canvas.draw()

    def draw(self):

        for i, f in enumerate(self.flds):
            if i < 2:
                try:
                    f.get()
                except:
                     messagebox.showerror("Error", "Please enter file name for %s"%(self.lbls[i][:-1]))
                if f.get() is "":
                    messagebox.showerror("Error", "Please enter file name for %s"%(self.lbls[i][:-1]))

        vlbls = ["Radius: Inner (m)","Radius: Outter (m)","Accuracy (C)", "Margin (C)"]            
        for i, v in enumerate(self.value):
            try:
                v.get()
            except:
                 messagebox.showerror("Error", "Please enter value for %s"%(vlbls[i]))
                 
        sinfo_file = self.flds[0].get()
        tdata_file = self.flds[1].get()
        radius_i = self.value[0].get()
        radius_o = self.value[1].get() 
        acuracy = self.value[2].get() 
        margin = self.value[3].get()
        

        self.td = HotRod(tdata_file, 
                         sinfo_file,
                         radius_inner=radius_i, 
                         radius_outer=radius_o, 
                         acuracy=acuracy, 
                         margin=margin)

        if not self.figure_drawn:
            self.time_slider.pack()
 
        self.time_slider.set(self.slider_val0)
        
        self.f.clear()
        ax = self.f.add_subplot(111,projection='3d')
            
        self.tr1 = self.td.xyz[:,:7]    
        self.dt = np.arange(0, len(self.td.t))
         
        self.td.max_responder() 
        
        # location of rod 1
        ax.plot(self.tr1[0,[0,6]], self.tr1[1,[0,6]], self.tr1[2,[0,6]] + self.td.relay_depth,
                'k-', 
                zdir = 'z',
                label = 'T1 sensors')
        # sensor locations
        ax.scatter(self.td.xyz[0,:], self.td.xyz[1,:], self.td.xyz[2,:] + self.td.relay_depth, 
                   zdir = 'z',  
                   s = 10, 
                   c = 'r', 
                   edgecolor = '')
        # heat source location
        ax.scatter(0, 0, self.td.relay_depth, 
                   zdir = 'z', 
                   s = 50, 
                   c = 'orange', 
                   edgecolor = '', 
                   label = 'heat source')
        # initial conditions
        ax.scatter([self.td.XYZ_max[0]], [self.td.XYZ_max[1]], 
                   [self.td.relay_depth + self.td.XYZ_max[2]],
                zdir = 'z', 
                c = 'b', 
                s = 100, 
                label = 'max responder')
        
        if hasattr(self.td, 'flux_xyz'):
            mag = np.sqrt(self.td.flux_xyz[0] ** 2 + self.td.flux_xyz[1] ** 2. + self.td.flux_xyz[2] ** 2.)    
            # scale the flow vector so that its length equals 0.05m (i.e., so that its visable)
            a = 0.05/mag        
            # calibrated vector
            ax.plot([0, self.td.flux_xyz[0]*a], [0, self.td.flux_xyz[1]*a], 
                    [self.td.relay_depth, self.td.relay_depth + self.td.flux_xyz[2]*a],
                    zdir = 'z', 
                    c = 'k',
                    label = 'flow')
        
        # observed temperature
        self.sc = ax.scatter(self.td.xyz[0,:], self.td.xyz[1,:], 
                             self.td.xyz[2,:] + self.td.relay_depth, 
                             zdir = 'z',  
                             s = 500 * self.td.temp_obs[self.slider_val0, :],
                            c = 'red', 
                            edgecolor = '', 
                            label = 'sensors')
        
        ax.set_xlabel('X (m)')
        ax.set_ylabel('Y (m)')
        ax.set_zlabel('Z (m)')
            
        self.figure_drawn = True

#        if self.hr
#            self.hr.angle2flux()
    
class Application(Tk.Frame):
    
    def __init__(self, master):
        Tk.Frame.__init__(self) # to initialise the class Tk.Frame    
        self.master = master
        self.master.title(u"HotRod")

        self.dream_proc = None
        self.tab2 = None

        for row_index in range(10):
            self.master.rowconfigure(row_index, weight=1)
            for col_index in range(10):
                self.master.columnconfigure(col_index, weight=1)     
        self.master.minsize(width=600, height=600)   
        
        self.nb = ttk.Notebook(self.master)
        self.nb.grid()
        self.tab1 = ttk.Frame(self.nb)
        self.nb.add(self.tab1, text = "Console")       
         
        self._create_input_panels()

        # catch user exit request
        self.bind("<Destroy>", self.quit_gui)

    def quit_gui(self, e):
        # Kill the dream process if necessary
        if self.dream_proc is not None:
            self.dream_proc.terminate()
        self.destroy()
        
    def _create_input_panels(self):
        self.frame1 = Tk.Frame(self.tab1, background='green')
        self.frame1.grid(row=0, column=0, columnspan=7, 
                         rowspan=10, sticky='nwse')
        
        for row_index in range(10):
            self.frame1.rowconfigure(row_index, weight=1)
            for col_index in range(7):
                self.frame1.columnconfigure(col_index, weight=1)
      
        self._create_input_pane(self.frame1)
        self.lf.grid(column=0, row=0, 
                       columnspan=7, rowspan = 1,
                       sticky="nsew")

        #######################################################################
        self.frame2 = Tk.Frame(self.tab1, background='red')
        self.frame2.grid(row=0, column=7, columnspan=3, rowspan=10, sticky='nwse')
        
        for row_index in range(10):
            self.frame2.rowconfigure(row_index, weight=1)
            for col_index in range(3):
                self.frame2.columnconfigure(col_index, weight=1)

        
        self._create_setup_pane(self.frame2)
        self.lfSetup.grid(column=0, row=0, rowspan=1,  columnspan=3, sticky="nsew")
      
        self._create_calibration_panes(self.frame2)   
        self.lfInput.grid(column=0, row=1, columnspan=3, sticky="nsew")
        self.lfType.grid(column=0, row=2, columnspan=3, sticky="nsew")
        self._param_option_panes()
        self._calibration_option_panes()
        #######################################################################

        self.hr3d = hotrod3Dfigure(self.frame1, self.flds, self.setup, self.lbls)
        self.button = Tk.Button(self.frame1, text="Graph Data", 
                                    command = self.hr3d.draw)

        self.button.grid()
        
        self.solveButton = ttk.Button(self.frame2, 
                                      text = 'Calibrate', 
                                      command = self._solve_frame)  
        self.solveButton.grid(row=3, column = 1)  
        
    def _create_input_pane(self, parent):

        self.lf = Tk.LabelFrame(parent, text='File Information', labelanchor="n")
        self.lbls = ['Sensor Info:','Temp. Data:','Output Location:']
        default_vals = ['./inputs/hotrod_sensor_array_updated.dat',"./inputs/templog_H_R2.csv", '']
        self.flds = []

        for i in range(3):
            f = Tk.Frame(self.lf)
            f.grid(row=i+1, column=0)
             
            lbl = Tk.Label(f,
                           justify = Tk.LEFT,
                           text=self.lbls[i],
                            width=15)
            e = Tk.Entry(f, width=85)
            e.insert(0,default_vals[i])
            self.flds.append(e) # save reference to field
             

            lbl.grid(row=i+1, column=0)
            e.grid(row=i+1, column=1, sticky='nsew')
             
        self.flds[0].focus_set()

    def _create_setup_pane(self, parent):
        self.lfSetup = Tk.LabelFrame(parent, 
                           text='Setup', 
                           labelanchor="n")
        
        degree = u"\u00b0"
        self.vlbls = ["Radius (m)","Accuracy ("+degree+"C)", "Margin ("+degree+"C)"]
        default_vals = ['0.028', '0.047','0.06', '1']
        
        f = Tk.Frame(self.lfSetup)
        f.grid(column=0, row=0, columnspan=5, sticky="nsew")
        
        Tk.Label(f, text="Inner").grid(row=0, column=1)
        Tk.Label(f, text="Outer").grid(row=0, column=2)
        
        self.setup = []

        c = -1
        for i in range(len(self.lbls)):
            lbl = Tk.Label(f, text=self.vlbls[i], anchor="e").grid(row = i+1, column=0)
            
            v1 = Tk.DoubleVar()
            eVal = Tk.Entry(f, width=10, text = v1).grid(row=i+1, column=1)
            v1.set(default_vals[i])
            self.setup.append(v1) # save reference to field
            
            if i == 0:
                v1 = Tk.DoubleVar()
                eVal = Tk.Entry(f, width=10, text = v1).grid(row=i+1, column=2)
                v1.set('0.047')      
                self.setup.append(v1) # save reference to field
        

        
    def _create_calibration_panes(self, parent):
        self.lfInput = Tk.LabelFrame(parent, 
                           text='Calibration: input', 
                           labelanchor="n")
        
        self.lfType = Tk.LabelFrame(parent, 
                           text='Calibration: type', 
                           labelanchor="n") 

        
    def cb_check1(self):
        if self.Checkb1.get() == 0:
            self.Checkb2.set(value = 1) 
            self.Checkb1.set(value = 0)
        else:
            self.Checkb1.set(value = 1) 
            self.Checkb2.set(value = 0) 
                      
    def cb_check2(self):
        if self.Checkb2.get() == 1:
            self.Checkb1.set(value = 0)  
        else:
            self.Checkb1.set(value = 1) 
            self.Checkb2.set(value = 0)


    def naccheck(self, entry, var):
        if var.get() == 0:
            entry.configure(state='disabled')
        else:
            entry.configure(state='normal')
            
    def _param_option_panes(self):
        
        plusMinus, sigma, rhoC = u"\u00b1", u'\u03C3', u"\u03c1"+"c" 
        textTitle = plusMinus + r" or " + sigma
        
        for row_index in range(8):
            self.lfInput.rowconfigure(row_index, weight=1)
            for col_index in range(4):
                self.lfInput.columnconfigure(col_index, weight=1) 
        
        self.Checkb1 = Tk.IntVar(value = 1)    
        b1 = Tk.Checkbutton(self.lfInput, 
                            text = 'Option 1', 
                            disabledforeground="green",
                            variable = self.Checkb1,
                            command=self.cb_check1,
                            onvalue=1,
                            offvalue=0)

        self.Checkb2 = Tk.IntVar(value = 0)   
        b2 = Tk.Checkbutton(self.lfInput, 
                            text = 'Option 2', 
                            disabledforeground="green",
                            variable = self.Checkb2,
                            command=self.cb_check2,
                            onvalue=1,
                            offvalue=0)
               
        b1.grid(row = 3, column = 0, columnspan = 2)
        b2.grid(row = 3, column = 2, columnspan = 2)    
        
        Tk.Label(self.lfInput, text="Initial Value").grid(row=0, column=1)
        Tk.Label(self.lfInput, text=textTitle).grid(row=0, column=2)
        Tk.Label(self.lfInput, text="Distribution").grid(row=0, column=3)
        
        self.lbls = ['ln(q)', rhoC]
        default_vals = ['-5', '2.75E+6']
        default_width = ['3','1.25E+06']
        distrChoice = ["","Uniform", "Normal"]
        self.ini, self.widths, self.Distr = [], [], []
        
        for i in range(len(self.lbls)):
            lbl = Tk.Label(self.lfInput, text=self.lbls[i],
                       font = "Calibri 10 italic").grid(row = i+1)
   
            v1 = Tk.DoubleVar()
            eVal = Tk.Entry(self.lfInput, width=10, text = v1).grid(row=i+1, column=1)
            v1.set(default_vals[i])
            
            v2 = Tk.DoubleVar()
            eWidth = Tk.Entry(self.lfInput, width=10, text = v2).grid(row=i+1, column=2)
            v2.set(default_width[i])
            
            self.ini.append(v1) # save reference to field
            self.widths.append(v2)
 
            distr = Tk.StringVar()    
            distr.set(distrChoice[1])
            bU1 = ttk.OptionMenu(self.lfInput, 
                                distr,
                                *distrChoice
                                )
            bU1.grid(row = i+1, column = 3)
        
            self.Distr.append(distr)        
        
         
        self.opt1 = Tk.LabelFrame(self.lfInput, text = "Option 1")   
        self.opt1.grid(row = 4, column=0, columnspan=8, rowspan = 5, sticky='nsew')
        
        for row_index in range(2):
            self.opt1.rowconfigure(row_index, weight=1)
            for col_index in range(5):
                self.opt1.columnconfigure(col_index, weight=1)        

        self.lbls1 = ['D'+u"\u2097", 'D'+u"\u209c"]
        self.lbls2 = ['k'] 
                 
        default_vals1 = ['1.0E-06','1.0E-06']
        default_width1 = ['1.0E-06','1.0E-06']

        default_vals2 = ['2.5']
        default_width2 = ['1.25']  

               
        self.ini1, self.width1, self.opt1Distr = [], [], []

        for i in range(len(self.lbls1)):
            lbl = Tk.Label(self.opt1, text=self.lbls1[i],
                       font = "Calibri 10 italic").grid(row = i+2, column=0, 
                                                        columnspan=2, padx=8)
   
            v1 = Tk.DoubleVar()
            eVal = Tk.Entry(self.opt1, width=10, text = v1).grid(row=i+2, column=2)
            v1.set(default_vals1[i])
            
            v2 = Tk.DoubleVar()
            eWidth = Tk.Entry(self.opt1, width=10, text = v2).grid(row=i+2, column=3)
            v2.set(default_width1[i])
            
            self.ini1.append(v1) # save reference to field
            self.width1.append(v2)
 
            distr = Tk.StringVar()    
            distr.set(distrChoice[1])
            bU1 = ttk.OptionMenu(self.opt1, 
                                distr,
                                *distrChoice
                                )
            bU1.grid(row = i+2, column = 4)
        
            self.opt1Distr.append(distr)

        self.opt2 = Tk.LabelFrame(self.lfInput, text = "Option 2")   
        self.opt2.grid(row = 10, column=0, columnspan=8, rowspan = 5, sticky='ew')
        
        for row_index in range(2):
            self.opt2.rowconfigure(row_index, weight=1)
            for col_index in range(5):
                self.opt2.columnconfigure(col_index, weight=1)       
                
 
        lbl = Tk.Label(self.opt2, text=self.lbls2[0],
                       font = "Calibri 10 italic").grid(row = i+2, column=0, 
                                                        columnspan=2, padx=10)
        
        self.ini2 = Tk.DoubleVar()
        eVal = Tk.Entry(self.opt2, width=10, text = self.ini2).grid(row=i+2, column=2)
        self.ini2.set(default_vals2[0])
        
        self.width2 = Tk.DoubleVar()
        eWidth = Tk.Entry(self.opt2, width=10, text = self.width2).grid(row=i+2, column=3)
        self.width2.set(default_width2[0])

        self.distr = Tk.StringVar()    
        self.distr.set(distrChoice[1])
        bU1 = ttk.OptionMenu(self.opt2, 
                            self.distr,
                            *distrChoice
                            )
        bU1.grid(row = i+2, column = 4)

            
    def _calibration_option_panes(self):
       
        for column in range(5):
            self.lfType.columnconfigure(column, weight=1)
            self.lfType.rowconfigure(column, weight=1)    

        solverChoice = ["","TNC", "DREAM"]
        Tk.Label(self.lfType, text="Solver").grid(row=0, column=0)         
        self.solv = Tk.StringVar()    
        self.solv.set(solverChoice[2])
        b1 = ttk.OptionMenu(self.lfType, 
                            self.solv,
                            *solverChoice
                            )
        b1.grid(row = 0, column = 1)
        
        self.DREAMlbls = ["# chains", "Burn length", "# pairs"]
        DREAMdefault = ["10", "300", "1"]

        self.DREAMparam = []

        for i in range(len(self.DREAMlbls)):
            
            lbl = Tk.Label(self.lfType, text=self.DREAMlbls[i],
                           font = "Calibri 10").grid(row = i+1, column=0)
            
            v1 = Tk.IntVar()
            eVal = Tk.Entry(self.lfType, width=10, text = v1).grid(row=i+1, column=1)
            v1.set(DREAMdefault[i])
                       
            self.DREAMparam.append(v1) # save reference to field      
            
     

    def check_input(self, choice):
        
        plusMinus, sigma = u"\u00b1", u'\u03C3'
        allInputs = [self.ini, self.widths]
        labels=self.lbls
        for sets in allInputs:
            if sets == 0:
                text = "initial"
            else:
                text = plusMinus + r" or " + sigma
            for var in sets:
                try:
                    var.get()
                except:
                     messagebox.showerror("Missing Value", "Please enter %s value for %s"%(text, self.lbls[i][:-1]))
                if var.get() is "":
                    messagebox.showerror("Missing value", "Please enter value for %s"%(self.lbls[i][:-1]))
        
        if choice == 1:
            allInputs = [self.ini1, self.width1]
            labels=self.lbls1
        else:
            allInputs = [self.ini2, self.width2]
            labels=self.lbls2
            
#        print(labels)
        for sets in allInputs:
            if sets == 0:
                text = "initial"
            else:
                text = plusMinus + r" or " + sigma
            for var in sets:
                try:
                    var.get()
                except:
                     messagebox.showerror("Missing Value", "Please enter %s value for %s"%(text, self.lbls[i][:-1]))
                if var.get() is "":
                    messagebox.showerror("Missing value", "Please enter value for %s"%(self.lbls[i][:-1]))
        
    def _solve_frame(self):
        # add tab for solver output
        if self.tab2 is None:
            self.tab2 = ttk.Frame(self.nb)
            self.nb.add(self.tab2, text = "Ouput")
        # run the calibration code
        self._solve()
        
    def _solve(self):
        
        choice = self.Checkb1.get()
        if choice == 0:
            choice = 2
        self.check_input(choice)
        angles = ['theta','phi']
        
        self.names, self.avg, self.width, self.distr = [], [], [], []
        for i, ii in enumerate(self.ini):
            self.avg.append(float(ii.get()))
            self.distr.append(self.Distr[i].get())
            self.width.append(float(self.widths[i].get()))        
        # add the optional parameters to the list    
        if choice == 1:
            print("Option 1")
            self.names = self.lbls + angles + self.lbls1
            for i, ii in enumerate(self.ini1):
                self.avg.append(float(ii.get()))
                self.distr.append(self.opt1Distr[i].get())
                self.width.append(float(self.width1[i].get()))
        else:
            print("Option 2")
            self.names = self.lbls + angles + self.lbls2
            for i, ii in enumerate(self.ini2):
                self.avg.append(float(ii.get()))
                self.distr.append(self.opt2Distr[i].get())
                self.width.append(float(self.width2[i].get()))


        self.avg = self.avg[:2] +[0,0] +self.avg[2:]
        self.width = self.width[:2] +[0,0]+self.width[2:]
        self.distr = self.distr[:2] + ["Uniform", "Uniform"] + self.distr[2:]
        


        calType = self.solv.get()
        if calType == "DREAM":
            self.DREAMprm = []
            for dm in self.DREAMparam:
                self.DREAMprm.append(int(dm.get()))

        
        sinfo_file = self.flds[0].get()
        tdata_file = self.flds[1].get()
        radius_i = self.setup[0].get()
        radius_o = self.setup[1].get() 
        acuracy = self.setup[2].get() 
        margin = self.setup[3].get()
        
        self.hr = HotRod(tdata_file, 
                         sinfo_file,
                         self.names,
                         self.avg,
                         self.width,
                         self.distr,
                         paramChoice = choice,
                         calChoice = calType,
                         radius_inner=radius_i, 
                         radius_outer=radius_o, 
                         acuracy=acuracy, 
                         margin=margin)

        
        self.hr.get_param_stats()
        self.hr.paramChoice = choice
        self.hr.calChoice = calType    
        print(self.hr.calChoice) 
        if  calType == 'TNC':                   
            #self.hr.solve(self.names)
            print("TNC is disabled, use DREAM instead")

        def dream_worker():
            # set up markov chain; play around with burn in and chains to get better fit
            nchains = self.DREAMprm[0]
            nburn = self.DREAMprm[1]
            npairs = self.DREAMprm[2]
            npars = len(self.hr.par_names)
            self.D = DREAM(nchains, nburn = nburn, npairs = 1)
            self.D.sampler(self.hr)

        TICK_INTERVAL = 100 # ms

        def wait_for_dream():
            if self.dream_proc is not None and self.dream_proc.is_alive():
                self.master.after(TICK_INTERVAL, wait_for_dream)
            else:
                self.dream_proc.join()
                self.dream_proc = None

        if  calType == 'DREAM':
            if self.dream_proc is not None:
                print("DREAM already running")
            else:
                self.dream_proc = mp.Process(target=dream_worker)
                self.dream_proc.start()
                self.master.after(TICK_INTERVAL, wait_for_dream)
#------------------------------------------------------------------------------    
            
def main():
    
    root = Tk.Tk()    
    root.rowconfigure(0, weight=1)
    root.columnconfigure(0, weight=1)
    root.resizable(0,0)
    
    #hr = hotrod3Dfigure(root)        
    hotrod = Application(root)
    root.mainloop()
#------------------------------------------------------------------------------    
#------------------------------------------------------------------------------ 

if __name__ == '__main__':
    main()

    
