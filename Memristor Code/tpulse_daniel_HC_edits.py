from glob import glob
import tkinter as tk
from tkinter import filedialog, ttk, font
import matplotlib.backends.backend_tkagg as mplTk
import matplotlib.pyplot as plt
from matplotlib import rc, cm
from datetime import datetime
import os
from instr_inter import *
from time import sleep
from time import perf_counter as timer
import time
import numpy as np
from scan import meas_current
# from fast_current_scan import meas_current_
from tpulse_current_scan import meas_current_

def select_dir(*args):
   path= filedialog.askdirectory(title="Select a Folder")
   save_dir_val.set(path)

# Helper function to parse a comma seperated input string
def parse_csl(str):
    sub_strs = str.split(',')
    return np.array([float(sub.strip()) for sub in sub_strs])

def replot():  
    global meas_bias_v
    ax_meas.cla()
    ax_meas.set_xlabel('Measurement #')
    ax_meas.set_ylabel('Current (A)')
    if meas_bias_v is not None:
        ax_meas.set_title(f'Measurement Bias = {meas_bias_v:.3}V')
    ax_in.cla()
    ax_in.set_xlabel('Measurement #')
    ax_in.set_ylabel('Pulse Voltage (V)')

    meas_local_idx = 0
    #pulse_local_idx = 0

    # Just making sure the currents and pulse_hist are not None
    if currents is not None and pulse_hist is not None:
        for ii in range(min(currents.shape[0], pulse_hist.shape[0])):
            ax_meas.axvline(x=meas_local_idx-.5, c='k')
            v = pulse_hist[ii,0]
            for jj in range(currents.shape[1]):
                if currents[ii, jj] != np.NAN:
                    if v > 0:
                        ax_meas.plot(meas_local_idx, currents[ii, jj], 'ro')
                    else:
                        ax_meas.plot(meas_local_idx, currents[ii, jj], 'bo')
                    # Plot Pulses
                    if jj == 0:
                        plot_square(meas_local_idx, v)
                    else:
                        plot_zero(meas_local_idx, v)
                    meas_local_idx += 1

    fig_meas.tight_layout()
    fig_in.tight_layout()
    canvas_in.draw()
    canvas_meas.draw()
    canvas_in.flush_events()
    canvas_meas.flush_events()
    
def scale_change(*args):
    global plot_yscale_clicked
    try:
        val = plot_scale_clicked.get()
    except:
        return

    if val == 'symlog':
        plot_linthresh_meas.grid(column=2, row=6, sticky=tk.E, padx=5, pady=5)
        plot_linthresh_label.grid(column=0, row=6, sticky=tk.W, padx=5, pady=5)
        try:
            linthresh = plot_linthresh_val.get()
        except:
            return
        ax_meas.set_yscale(val, linthresh=linthresh)
    else:
        replot()
        plot_linthresh_meas.grid_remove()
        plot_linthresh_label.grid_remove()
        ax_meas.set_yscale(val)

    fig_meas.tight_layout()
    canvas_meas.draw()

def plot_square(idx, v):
    if v > 0:
        c = 'r'
    else:
        c = 'b'
    x = [idx-1, idx-.75, idx-.75, idx-.25, idx-.25, idx]
    y = [0, 0, v, v, 0, 0]
    ax_in.plot(x, y, c=c)

# plots the pulse graph
def plot_zero(idx, v):
    if v > 0:
        c = 'r'
    else:
        c = 'b'
    x = [idx-1, idx]
    y = [0, 0]
    ax_in.plot(x, y, c=c)

def save_data():
    lna_gain = lna_gain_mode_clicked.get().replace(' ', '_')
    lna_filter = lna_filter_clicked.get().replace(' ', '_')
    if lna_filter != 'None':
        filter_str = f'{lna_gain}_LP{lna_filter}{lna_filter_clicked.get()}Hz'
    else:
        filter_str = f'{lna_gain}_None'
    device_name = device_name_val.get()
    device_idx = device_idx_val.get()
    dir_name = save_dir_val.get()
    integ_time=dmm_integ_clicked.get()
    title = f'{device_name}_{device_idx}_{filter_str}_Integ{integ_time}'
    print('stop')
    # Create dir if it does not exist
    if not os.path.exists(dir_name):
        os.makedirs(dir_name, exist_ok=True)
    fig_meas.savefig(os.path.join(dir_name, f'{title}.png'))
    fig_in.savefig(os.path.join(dir_name, f'{title}_pulse_figure.png'))

    # Save the raw data to a directory
    headers = ['pulse_v', 'pulse_width', 'num_applied', 'meas_v'] + [f'i_{ii}' for ii in range(currents.shape[1])] # it is getting the number of measurements
    header_str = ','.join(headers)
    print("Final Pulse History: \n", pulse_hist)
    print("Unedited Currents Array\n: ", currents)
    print("Saved Currents Array\n: ", currents[-pulse_hist.shape[0]:])
    save_arr = np.hstack((pulse_hist, 
                          meas_bias_v*np.ones((pulse_hist.shape[0],1)),
                          currents[-pulse_hist.shape[0]:]))
    np.savetxt(os.path.join(dir_name, f'{title}.csv'), save_arr, 
                delimiter=",",header=header_str)

# This does the sweeping
def pulse_train(v, pulse_on_period, burst_period, 
                pulse_number, trig_freq, curr_min, curr_max):
    
    print(f"Pulse Train: {v}V")
    # Set up
    relay.switch_relay(relay.NONE)
    PGen.config_pulse(pulse_voltage=v, pulse_width='in')
    # Set the LNA bias
    LNA.set_bias(round(meas_bias_v*1000))
    sleep(0.5)
    num_meas = dmm_num_val.get()
    global pulse_idx, meas_idx, currents, e_stop

    # Set up burst
    DGen.config_burst(pulse_on_period, burst_period, pulse_number, 
                      trig_freq, 2, bnc=1, channel_start=2, channel_end=3)
    DGen.config_burst(pulse_on_period, burst_period, pulse_number, 
                      trig_freq, 2, bnc=2, channel_start=4, channel_end=5)
    PGen.config_pulse(pulse_voltage=v, pulse_width='in')
    sleep(0.5)
    # sens = lna_sens_clicked.get() # <! in pulse_train
    sens = LNA.sens

    # Start with pulse
    relay.switch_relay(relay.DGEN)
    DGen.trig_burst() # Triggers burst
    relay.switch_relay(relay.LNA)
    
    plot_square(meas_idx, v)
    ax_meas.axvline(x=meas_idx-.5, c='k')

    # plots the max curr and min curr
    ax_meas.plot([ii for ii in range(meas_idx, meas_idx+num_meas+1)],
                 [curr_min for jj in range(0, num_meas+1)], 'g')
    ax_meas.plot([ii for ii in range(meas_idx, meas_idx+num_meas+1)],
                [curr_max for jj in range(0, num_meas+1)], 'g')
    
    # Measurements after each burst, does plotting
    for kk in range(num_meas):

        # This is the actual measuring portion below, meas_current is in scan.py file, why does it return sens?
        currents[pulse_idx, kk], sens = meas_current(sens, LNA, DMM,
                                                    LNA_gui_var=curr_var, 
                                                    Sens_gui_var=sens_var)

        # all plotting
        if v > 0:
            ax_meas.plot(meas_idx, currents[pulse_idx, kk], 'ro')
        else:
            ax_meas.plot(meas_idx, currents[pulse_idx, kk], 'bo')
        if kk > 0:
            plot_zero(meas_idx, v)

        # ax_meas.plot(meas_idx, curr_min, 'go')
        # ax_meas.plot(meas_idx, curr_max, 'go')

        # Pause to allow plots and progress bar to update
        fig_in.tight_layout()
        canvas_in.draw()
        canvas_in.flush_events()
        fig_meas.tight_layout()
        canvas_meas.draw()
        canvas_meas.flush_events()
        sleep(0.02)

        meas_idx += 1
    pulse_idx += 1
    new_current_row = np.ones((1,currents.shape[1]))
    new_current_row[:] = np.nan
    currents = np.vstack((currents, new_current_row))
    LNA.set_bias(0)

#-----------------------------------------------------------------------------------
# Attempt at adding a de-noise code
#-----------------------------------------------------------------------------------
def denoise(v, pulse_on_period, burst_period, 
                pulse_number, trig_freq):
    
    print(f"Denoise: {v}V")
    # Set up
    relay.switch_relay(relay.NONE)
    PGen.config_pulse(pulse_voltage=v, pulse_width='in')
    # Set the LNA bias
    LNA.set_bias(round(meas_bias_v*1000))
    sleep(0.5)
    global pulse_idx, meas_idx, currents, e_stop

    # Set up burst
    DGen.config_burst(pulse_on_period, burst_period, pulse_number, 
                      trig_freq, 2, bnc=1, channel_start=2, channel_end=3)
    DGen.config_burst(pulse_on_period, burst_period, pulse_number, 
                      trig_freq, 2, bnc=2, channel_start=4, channel_end=5)
    PGen.config_pulse(pulse_voltage=v, pulse_width='in')
    sleep(0.5)

    # Start with pulse
    relay.switch_relay(relay.DGEN)
    DGen.trig_burst() # Triggers burst
    relay.switch_relay(relay.LNA)
  
    LNA.set_bias(0)
#-----------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------


def measure_loop(first_v, res_min, res_max, max_tries, pulse_on_period, 
                 burst_period, pulse_num, trig_freq):
    global currents, pulse_idx, meas_idx, pulse_hist, meas_bias_v, e_stop
    curr_min = meas_bias_v/res_max
    curr_max = meas_bias_v/res_min
    sens = lna_sens_clicked.get() # <! in measure_loop
    relay.switch_relay(relay.LNA)
    # First measure the current to see what's up with the device
    LNA.set_bias(round(meas_bias_v*1000))
    curr, sens = meas_current(sens, LNA, DMM,
                            LNA_gui_var=curr_var, 
                            Sens_gui_var=sens_var)
    res = abs(meas_bias_v/curr)
    attempt = 0
    amplitude = first_v
    add = 0.5
    pulse_config_arr = np.ones((1, 3))
    # If resistance is not within range
    while ((res < res_min or res > res_max) and attempt < max_tries):
        if e_stop:
            relay.switch_relay(relay.NONE)
            LNA.set_bias(0)
            print("Emergency Stopped!")
            break
        print("Resistance: ", res)
        if (attempt == 0):
            pulse_train(amplitude, pulse_on_period, burst_period, pulse_num, trig_freq,
                        curr_min, curr_max)\

        else:
            # if the resistance is close to the resistance range, decrease the increment
            if (1 > res_max/res > 0.5 or 1 < res_min/res < 2):
                add = 0.20
            else:
                add = 0.50

            # go through the measurement loop
            if (res < res_min):
                # if the resistance is too low but the pulse is positive, must make pulse negative
                if (amplitude > 0):
                    amplitude = -add
                    pulse_train(amplitude, pulse_on_period, burst_period, pulse_num, trig_freq,
                                curr_min, curr_max)
                # if previous pulse changed the current enough, repeat
                elif (currents.shape[0] >= 3 and abs(currents[-2,-1]/currents[-3,-1]) < 4/7) or (amplitude >= 10 or amplitude <= -8):
                    pulse_train(amplitude, pulse_on_period, burst_period, pulse_num, trig_freq,
                                curr_min, curr_max)
                else:
                    amplitude -= add
                    pulse_train(amplitude, pulse_on_period, burst_period, pulse_num, trig_freq,
                                curr_min, curr_max)
            elif (res > res_max):
                # if the resistance is too high but the pulse is negative, must make pulse positive
                if (amplitude < 0):
                    amplitude = add
                    pulse_train(amplitude, pulse_on_period, burst_period, pulse_num, trig_freq,
                                curr_min, curr_max)
                elif (currents.shape[0] >= 3 and abs(currents[-2,-1]/currents[-3,-1]) > 1.75) or (amplitude >= 10 or amplitude <= -8):
                    pulse_train(amplitude, pulse_on_period, burst_period, pulse_num, trig_freq,
                                curr_min, curr_max)
                else:
                    amplitude += add
                    pulse_train(amplitude, pulse_on_period, burst_period, pulse_num, trig_freq,
                                curr_min, curr_max)
        
        # Add array to pulse history
        pulse_config_arr[0, 0] = amplitude
        pulse_config_arr[0, 1] = pulse_on_period
        pulse_config_arr[0, 2] = pulse_num
        res = abs(meas_bias_v/currents[-2, -1])
        if attempt == 0:
            print("Should only print once...")
            pulse_hist = pulse_config_arr
        else:
            print("config arr: ", pulse_config_arr)
            print("small pulse hist: ", pulse_hist)
            pulse_hist = np.vstack((pulse_hist, pulse_config_arr))
        attempt += 1
        print("Number of attempts: ", attempt)
        print("Pulse History Update: \n", pulse_hist)
    
    # last row is an extra row of nan        
    currents = currents[0:-1]

def meas_something(DMM):
    stuff = DMM.read()
    return stuff


#followup_measurement
def follow_loop(res_min=1.02e7, res_max=1.04e7):
    global currents, pulse_idx, meas_idx, pulse_hist, meas_bias_v, e_stop
    #e_stop = False
    dir_name = save_dir_val.get()
    sens = lna_sens_clicked.get() # <! in follow_loop
    if sens < 9:
        orig_sens = sens
        sens = 9

    lna_gain = lna_gain_mode_clicked.get().replace(' ', '_')
    lna_filter = lna_filter_clicked.get().replace(' ', '_')
    if lna_filter != 'None':
        filter_str = f'{lna_gain}_LP{lna_filter}{lna_filter_clicked.get()}Hz'
    else:
        filter_str = f'{lna_gain}_None'
    device_name = device_name_val.get()
    device_idx = device_idx_val.get()

    relay.switch_relay(relay.LNA)

    # get the added_time, sample_count, and trig_delay from interface
    added_time = follow_time_val.get()
    sample_count = follow_num_val.get()
    start_delay = start_delay_val.get()
    auto_range_status = auto_range_on.get()
    trig_delay = 1

    meas_bias_v = lna_bias_clicked.get()
    print("LNA bias:", meas_bias_v)
    LNA.set_bias(round(meas_bias_v*1000))

    # First measure the current to see what's up with the device
    _, sens = meas_current_(sens, LNA, DMM,
                            LNA_gui_var=curr_var, 
                            Sens_gui_var=sens_var)
    print('sens: ', sens)

    # fetchs the integ time associated with the returned sens
    # integ_time = DMM.integ_times_table[sens] 
    integ_time = .2
    print("Integration time:", integ_time)
    DMM.set_integ_time(integ_time)

    # calculates the data interval based on integration time and follow time
    integ_time_real = integ_time / 60
    data_interval = integ_time_real*1.5 + added_time
    print("Sample time:", data_interval)

    curr_str, elapsed_time = DMM.timer_sample(sample_count=sample_count, start_delay=start_delay, 
                                source_timer=data_interval, trig_delay=trig_delay,
                                auto_range_off=(not auto_range_status))
    # curr_str, elapsed_time = DMM.delay_sample(sample_count=sample_count, start_delay=added_time)
    print("Elapsed time:", elapsed_time)

    # processing data
    curr_str = curr_str[0:-1] # removes trailing \n (newline character)
    curr_float = np.array([float(c) for c in curr_str.split(',')])
    print("Sensitivity multiplier:", LNA.sens_table[sens])
    curr_vals = curr_float * LNA.sens_table[sens]
    res_vals = abs(meas_bias_v/curr_vals)
    time_data = np.array([trig_delay + data_interval*i 
                          for i in range(0, sample_count)])
    
    DMM.write('*RST') # maybe this will fix it?
    
    # Create dir if it does not exist
    if not os.path.exists(dir_name):
        os.makedirs(dir_name, exist_ok=True)

    # Save the raw data to a directory
    title = f'{device_name}_{device_idx}_{filter_str}_Integ{integ_time}_retentiondata'
    headers = ['resistance (ohms)', 'current (A)', 'time (s)', 
               'res min', 'res_max'] 
    header_str = ','.join(headers)
    # print("Time Data: \n", time_data)
    # print("Resistance value\n: ", res_vals)
    save_arr = np.transpose(
        np.vstack((res_vals, curr_vals, time_data, 
                   res_min*np.ones(sample_count), 
                   res_max*np.ones(sample_count))))
    np.savetxt(os.path.join(dir_name, f'{title}.csv'), save_arr, 
                delimiter=",",header=header_str)

    # Save some useful information about the scan to a file
    info_file = open(os.path.join(dir_name, f'{title}_info.txt'), "w")

    info_file.write("Useful scan information:\n\n")
    info_file.write(f"Original sensitivity: {orig_sens}\n")
    info_file.write(f"Sensitivity: {sens}\n")
    info_file.write(f"Integration time (*PLC): {integ_time}\n")
    info_file.write(f"Integ time real (s): {integ_time_real}\n")
    info_file.write(f"Sample time (s): {data_interval}\n")
    info_file.write(f"Elapsed time (s): {elapsed_time}\n")
    info_file.write(f"Sensitivity multiplier: {LNA.sens_table[sens]}\n")
    info_file.write(f"LNA measurement bias (V): {meas_bias_v}\n")
    info_file.write(f"Auto range on: {auto_range_status}\n")
    info_file.write(f"Start delay (s): {start_delay}\n")

    info_file.close()



# Biggest loop
def dynamic_scan(*args):
    number_looped = len(parse_csl(res_range_val.get()))
    loop_counter = 0
    global e_stop
    e_stop = False
    while loop_counter < number_looped/2 and not e_stop:
        print(loop_counter)
        print(e_stop)
        # Get pulse configs from GUI
        first_v = pulse_v_val.get()
        pulse_on_period = float(pulse_period_val.get())
        pulse_off_period = pulse_on_period
        pulse_per_burst = int(pulse_per_burst_val.get())
        pulse_number = int(pulse_per_burst_val.get())
        burst_period = (pulse_on_period+pulse_off_period)
        trig_period = pulse_per_burst*(burst_period)*2
        trig_freq = 1/trig_period
        res_range = parse_csl(res_range_val.get())
        max_tries = max_tries_val.get()
        res_min = res_range[2*loop_counter]
        res_max = res_range[2*loop_counter + 1]
        followup = follow_up.get()
        assert res_max > res_min, "1st argument must be less than 2nd argument"

        # Cease operating GUI
        save_dir_btn['state'] = tk.DISABLED
        pb.grid(column=0, row=2, columnspan=3, sticky='EW', padx=5, pady=5)
        stop_btn.grid(column=0, row=1, columnspan=3, sticky='EW', padx=5, pady=5)
        sens_label.grid(column=0, row=3, sticky='EW', padx=5, pady=5)
        curr_label.grid(column=2, row=3, sticky='EW', padx=5, pady=5)
        sleep(0.2)

        # Pull parameters from GUI to initialize instruments
        DMM.set_integ_time(dmm_integ_clicked.get())
        lna_gain = lna_gain_mode_clicked.get().replace(' ', '_')
        lna_filter = lna_filter_clicked.get().replace(' ', '_')
        LNA.set_filter(lna_gain, lna_filter, lna_filter_freq_clicked.get())
        meas_per_volt = dmm_num_val.get()

        global currents, pulse_idx, meas_idx, pulse_hist, meas_bias_v, single_var
        
        if single_var == True:
            first_v = s_pulse_v_val.get()
            max_tries = 1

        # Get and Set bias for the LNA
        new_bias_v = lna_bias_clicked.get()
        is_new_bias = meas_bias_v is None or new_bias_v != meas_bias_v


        # Create/Adjust currents array
        if hold_on.get() and currents is not None and not is_new_bias:
            if currents.shape[1] < meas_per_volt:       # checks if new meas_per_volt is greater
                # First add additional columns to the original array, hstack extends an np array in horizontal dimension
                filler_arr = np.empty(currents.shape[0], meas_per_volt-currents.shape[1])
                filler_arr[:] = np.nan
                currents = np.hstack(currents, filler_arr)
            # Array: bursts vs measurements
            sub_arr = np.empty((1, currents.shape[1]))
            sub_arr[:] = np.nan
            
            # extend empty cells in the array to prepare for next measurements in the pulse_train
            currents = np.vstack((currents, sub_arr))
        else:
            meas_bias_v = new_bias_v
            pulse_idx = 0
            meas_idx = 0
            ax_meas.cla()
            ax_in.cla()

            ax_meas.set_xlabel('Measurement #')
            ax_meas.set_ylabel('Current (A)')
            ax_meas.set_title(f'Measurement Bias = {meas_bias_v:.3}V')
            ax_meas.relim()
            ax_meas.autoscale_view()
            fig_meas.tight_layout()

            ax_in.set_xlabel('Measurement #')
            ax_in.set_ylabel('Pulse Voltage (V)')
            ax_in.relim()
            ax_in.autoscale_view()
            fig_in.tight_layout()

            # THIS IS WHERE currents IS INITIALIZED
            currents = np.empty((1, meas_per_volt))
            currents[:] = np.nan

        # Measure the current to see what's up with the device
        measure_loop(first_v, res_min, res_max, max_tries, pulse_on_period,
                    burst_period, pulse_number, trig_freq)
        
        #daniel edits -- feel free to disregard --------------------------
        # def denoise(v, pulse_on_period, burst_period, 
        #                 pulse_number, trig_freq):

        dn_v = 0.3
        dn_pulse_on_period = float(500e-9)
        dn_pulse_off_period = dn_pulse_on_period
        dn_pulse_per_burst = int(1000)
        dn_pulse_number = int(1000)
        dn_burst_period = (dn_pulse_on_period+dn_pulse_off_period)
        dn_trig_period = dn_pulse_per_burst*(dn_burst_period)*2
        dn_trig_freq = 1/dn_trig_period

        #can add a button later for the denoising
        #if denoise_check == True and not e_stop:
        #denoise(dn_v, dn_pulse_on_period, dn_burst_period, dn_pulse_number, dn_trig_freq)
        #-----------------------------------------------------------------
        
        if followup == True and not e_stop:
            follow_loop(res_min, res_max)
        
        save_data()

        # Increment index
        device_idx_val.set(device_idx_val.get()+1)
        save_dir_btn['state'] = tk.NORMAL
        pb.grid_remove()
        sens_label.grid_remove()
        curr_label.grid_remove()
        stop_btn.grid_remove()

        res = abs(meas_bias_v/currents[-2, -1])
        # print(currents[-2,-1])
        # print(res)
        if (res > res_min and res < res_max) or single_var==True:
            print("done!")
        else:
            print(f"Completed all {max_tries} tries, but missed the resistance range of {res_min} to {res_max} Ohms. Final resistance value: {res}")
        
        loop_counter += 1

def dynamic_scan_call():
    global single_var
    single_var = False
    dynamic_scan()

def single_pulse():
    global single_var
    single_var = True
    dynamic_scan()
    
# Intialize Instruments
relay = relay_inter.relay_inter()
relay.switch_relay(relay.NONE) # Close relay during setup
DGen = DG645.DG645()
DMM = A34410A.A34410A()
LNA = SR570.SR570()
PGen = AV1010B.AV1010B()
PGen.set_trigger('EXT')

# Initial global variables
currents = None
meas_idx = None
pulse_idx = None
pulse_hist = None
meas_bias_v = None
e_stop = False
single_var = False

root = tk.Tk()
root.state('zoomed')
root.title('Target Pulse GUI')

# Dashboard
dashboard_frame = tk.Frame(root)
dashboard_frame.grid(column=0, row=0, rowspan=4, sticky='ns')

## Device options
device_frame = tk.Frame(dashboard_frame, bd=2, width=450)
device_frame.grid(column=0, row=0,sticky=tk.W+tk.E)
device_frame.grid_columnconfigure(0, weight=1)
# Save Dir
save_dir_label = tk.Label(device_frame, text='Save Folder')
save_dir_label.grid(column=0, row=1, sticky=tk.W, padx=5, pady=5)
save_dir_btn = tk.Button(device_frame, text='Browse', command=select_dir)
save_dir_btn.grid(column=1, row=1, sticky=tk.E, padx=5, pady=5)
save_dir_val = tk.StringVar()
save_dir_entry = tk.Entry(device_frame, text=save_dir_val, width=40)
save_dir_entry.grid(column=0, row=2, columnspan=2, sticky=tk.W, padx=5, pady=5)
now = datetime.now()
date = now.strftime('%Y-%m-%d')
save_dir_val.set(os.path.join('.', 'Measurements','Pulse', date))
# Device Name
device_name_label = tk.Label(device_frame, text='Device Name')
device_name_label.grid(column=0, row=3, sticky=tk.W, padx=5, pady=5)
device_name_val = tk.StringVar()
device_name_meas = tk.Entry(device_frame, text=device_name_val, width=20)
device_name_meas.grid(column=1, row=3, sticky=tk.E, padx=5, pady=5)
device_name_val.set('FIB3_K3_1')

# Meas Idx
device_idx_label = tk.Label(device_frame, text='Index')
device_idx_label.grid(column=0, row=4, sticky=tk.W, padx=5, pady=5)
device_idx_val = tk.IntVar()
device_idx_meas = tk.Entry(device_frame, text=device_idx_val, width=20)
device_idx_meas.grid(column=1, row=4, sticky=tk.E, padx=5, pady=5)
device_idx_val.set(0)

## LNA Options
lna_frame = tk.Frame(dashboard_frame, bd=2, width=375)
lna_frame.grid(column=0, row=1, sticky=tk.W+tk.E)
lna_frame.grid_columnconfigure(0, weight=1)
# Label for the subframe
tk.Label(lna_frame, text='LNA', font=('Arial',12)).grid(column=0, row=0, sticky=tk.W, padx=5, pady=5)
# Gain Mode
lna_gain_mode_label = tk.Label(lna_frame, text='Gain Mode')
lna_gain_mode_label.grid(column=0, row=1, sticky=tk.W, padx=5, pady=5)
lna_gain_mode_options = [m.replace('_', ' ') for m in list(SR570.SR570.gain_modes.keys())]
lna_gain_mode_clicked = tk.StringVar()
lna_gain_mode_clicked.set(lna_gain_mode_options[2])
lna_gain_mode_menu = tk.OptionMenu(lna_frame, lna_gain_mode_clicked ,*lna_gain_mode_options)
lna_gain_mode_menu.config(width=len(max(lna_gain_mode_options)))
lna_gain_mode_menu.grid(column=1, row=1, sticky=tk.E, padx=5, pady=5)
# Filter Type
lna_filter_label = tk.Label(lna_frame, text='LP Cutoff')
lna_filter_label.grid(column=0, row=2, sticky=tk.W, padx=5, pady=5)
lna_filter_options = list(SR570.SR570.filters.keys())
lna_filter_options[0] = 'None' # Convert None to a string
lna_filter_clicked = tk.StringVar()
lna_filter_clicked.set(lna_filter_options[1])
lna_filter_menu = tk.OptionMenu(lna_frame, lna_filter_clicked ,*lna_filter_options)
lna_filter_menu.grid(column=1, row=2, sticky=tk.E, padx=5, pady=5)
# Filter cutoff frequency
lna_filter_freq_label = tk.Label(lna_frame, text='LP Cutoff Frequency (Hz)')
lna_filter_freq_label.grid(column=0, row=3, sticky=tk.W, padx=5, pady=5)
lna_filter_freq_options = [f'{f}' for f in SR570.SR570.filter_freqs]
lna_filter_freq_clicked = tk.DoubleVar()
lna_filter_freq_clicked.set(lna_filter_freq_options[3])
lna_filter_freq_menu = tk.OptionMenu(lna_frame, lna_filter_freq_clicked ,*lna_filter_freq_options)
lna_filter_freq_menu.grid(column=1, row=3, sticky=tk.E, padx=5, pady=5)
# Initial Sens
lna_sens_label = tk.Label(lna_frame, text='Initial Sensitivity')
lna_sens_label.grid(column=0, row=4, sticky=tk.W, padx=5, pady=5)
lna_sens_options = list(range(len(SR570.SR570.sens_table)))
lna_sens_clicked = tk.IntVar()
lna_sens_clicked.set(lna_sens_options[0])
lna_sens_menu = tk.OptionMenu(lna_frame, lna_sens_clicked ,*lna_sens_options)
lna_sens_menu.grid(column=1, row=4, sticky=tk.E, padx=5, pady=5)
# LNA bias during measurement
lna_bias_label = tk.Label(lna_frame, text='LNA Measurement Bias (V)')
lna_bias_label.grid(column=0, row=5, sticky=tk.W, padx=5, pady=5)
lna_bias_clicked = tk.DoubleVar()
lna_bias_clicked.set(-0.1)
lna_bias_menu =  tk.Entry(lna_frame, text=lna_bias_clicked, width=10)
lna_bias_menu.grid(column=1, row=5, sticky=tk.E, padx=5, pady=5)
# Offset correction
lna_offset = tk.BooleanVar()
lna_offset_check = tk.Checkbutton(lna_frame, text='Offset Correction', variable=lna_offset)
lna_offset_check.grid(column=0, row=6, columnspan=2, padx=5, pady=5)

## DMM Options
dmm_frame = tk.Frame(dashboard_frame, bd=2, width=375)
dmm_frame.grid(column=0, row=2, sticky=tk.W+tk.E)
dmm_frame.grid_columnconfigure(0, weight=1)
# Label for the subframe
tk.Label(dmm_frame, text='DMM', font=('Arial',12)).grid(column=0, row=0, sticky=tk.W, padx=5, pady=5)
# Integration time
dmm_integ_label = tk.Label(dmm_frame, text='Integration Time (* 1/60 s)')
dmm_integ_label.grid(column=0, row=1, sticky=tk.W, padx=5, pady=5)
dmm_integ_options = A34410A.A34410A.integ_times
dmm_integ_clicked = tk.DoubleVar()
dmm_integ_clicked.set(dmm_integ_options[4])
dmm_integ_menu = tk.OptionMenu(dmm_frame, dmm_integ_clicked ,*dmm_integ_options)
dmm_integ_menu.grid(column=1, row=1, sticky=tk.E, padx=5, pady=5)
# Measurements per voltage
dmm_num_label = tk.Label(dmm_frame, text='Measurements per voltage')
dmm_num_label.grid(column=0, row=2, sticky=tk.W, padx=5, pady=5)
dmm_num_val = tk.IntVar()
dmm_num_meas = tk.Entry(dmm_frame, text=dmm_num_val, width=5)
dmm_num_meas.grid(column=1, row=2, sticky=tk.E, padx=5, pady=5)
dmm_num_val.set(5)
# Initial measurements to discard per voltage
dmm_discard_label = tk.Label(dmm_frame, text='Initial discarded')
dmm_discard_label.grid(column=0, row=3, sticky=tk.W, padx=5, pady=5)
dmm_discard_val = tk.IntVar()
dmm_discard_meas = tk.Entry(dmm_frame, text=dmm_discard_val, width=5)
dmm_discard_meas.grid(column=1, row=3, sticky=tk.E, padx=5, pady=5)
dmm_discard_val.set(0)

## Pulse Parameters
pulse_frame = tk.Frame(dashboard_frame, bd=2, width=375)
pulse_frame.grid(column=1, row=0, sticky=tk.W+tk.E)
pulse_frame.grid_columnconfigure(0, weight=1)
# Label for the subframe
tk.Label(pulse_frame, text='Pulse Measurement', font=('Arial',12)).grid(column=0, row=0, columnspan=2, sticky=tk.W, padx=5, pady=5)
# Starting Pulse Voltage
pulse_v_label = tk.Label(pulse_frame, text='Starting Pulse Amplitude (V)')
pulse_v_label.grid(column=0, row=1, sticky=tk.W, padx=5, pady=5)
pulse_v_val = tk.DoubleVar()
pulse_v_meas = tk.Entry(pulse_frame, text=pulse_v_val, width=25)
pulse_v_meas.grid(column=1, row=1, sticky=tk.E+tk.W, padx=5, pady=5)
pulse_v_meas.grid_columnconfigure(1, weight=1)
pulse_v_val.set(2.0)
# Pulse period
pulse_period_label = tk.Label(pulse_frame, text='Pulse On Period (s)')
pulse_period_label.grid(column=0, row=2, sticky=tk.W, padx=5, pady=5)
pulse_period_val = tk.StringVar()
pulse_period_meas = tk.Entry(pulse_frame, text=pulse_period_val, width=25)
pulse_period_meas.grid(column=1, row=2, sticky=tk.E+tk.W, padx=5, pady=5)
pulse_period_meas.grid_columnconfigure(1, weight=1)
pulse_period_val.set('500e-9')
# Pulse per Burst
pulse_per_burst_label = tk.Label(pulse_frame, text='Pulses per burst')
pulse_per_burst_label.grid(column=0, row=3, sticky=tk.W, padx=5, pady=5)
pulse_per_burst_val = tk.StringVar()
pulse_per_burst_meas = tk.Entry(pulse_frame, text=pulse_per_burst_val, width=25)
pulse_per_burst_meas.grid(column=1, row=3, sticky=tk.E+tk.W, padx=5, pady=5)
pulse_per_burst_meas.grid_columnconfigure(1, weight=1)
pulse_per_burst_val.set('1000')
# Single Pulse Voltage
s_pulse_v_label = tk.Label(pulse_frame, text='Single Pulse Amplitude (V)')
s_pulse_v_label.grid(column=0, row=4, sticky=tk.W, padx=5, pady=5)
s_pulse_v_val = tk.DoubleVar()
s_pulse_v_meas = tk.Entry(pulse_frame, text=s_pulse_v_val, width=25)
s_pulse_v_meas.grid(column=1, row=4, sticky=tk.E+tk.W, padx=5, pady=5)
s_pulse_v_meas.grid_columnconfigure(1, weight=1)
s_pulse_v_val.set('5')


# Dynamic Parameters
dynamic_frame = tk.Frame(dashboard_frame, bd=2, width=375)
dynamic_frame.grid(column=1, row=1, sticky=tk.W+tk.E)
dynamic_frame.grid_columnconfigure(0, weight=1)
# Label for the subframe
tk.Label(dynamic_frame, text="Dynamic Targets", font=("Arial",12)).grid(column=0, row=0, columnspan=2, sticky=tk.W, padx=5, pady=5)
# Resistance Limit
res_range_label = tk.Label(dynamic_frame, text = 'Resistance (min, max)')
res_range_label.grid(column=0, row=1, sticky=tk.W, padx=5, pady=5)
res_range_val = tk.StringVar()
res_range_meas = tk.Entry(dynamic_frame, text=res_range_val, width=20)
res_range_meas.grid(column=1, row=1, sticky=tk.W, padx=5, pady=5)
res_range_meas.grid_columnconfigure(1, weight=1)
#res_range_val.set('9.1e8,1e10,2.44e7,2.5e7,9.1e8,1e10,2.44e7,2.5e7,9.1e8,1e10,2.44e7,2.5e7,9.1e8,1e10,2.44e7,2.5e7,9.1e8,1e10,2.44e7,2.5e7,9.1e8,1e10,2.44e7,2.5e7,9.1e8,1e10,2.44e7,2.5e7,9.1e8,1e10,2.44e7,2.5e7,9.1e8,1e10,2.44e7,2.5e7,9.1e8,1e10,2.44e7,2.5e7')
#res_range_val.set('2.44e7,2.5e7,3.02e7,3.12e7,3.99e7,4.15e7,5.85e7,6.22e7,1.1e8,1.23e8,9.1e8,1e10,2.44e7,2.5e7,3.02e7,3.12e7,3.99e7,4.15e7,5.85e7,6.22e7,1.1e8,1.23e8,9.1e8,1e10,2.44e7,2.5e7,3.02e7,3.12e7,3.99e7,4.15e7,5.85e7,6.22e7,1.1e8,1.23e8,9.1e8,1e10,2.44e7,2.5e7,3.02e7,3.12e7,3.99e7,4.15e7,5.85e7,6.22e7,1.1e8,1.23e8,9.1e8,1e10,2.44e7,2.5e7,3.02e7,3.12e7,3.99e7,4.15e7,5.85e7,6.22e7,1.1e8,1.23e8,9.1e8,1e10,2.44e7,2.5e7,3.02e7,3.12e7,3.99e7,4.15e7,5.85e7,6.22e7,1.1e8,1.23e8,9.1e8,1e10,2.44e7,2.5e7,3.02e7,3.12e7,3.99e7,4.15e7,5.85e7,6.22e7,1.1e8,1.23e8,9.1e8,1e10,2.44e7,2.5e7,3.02e7,3.12e7,3.99e7,4.15e7,5.85e7,6.22e7,1.1e8,1.23e8,9.1e8,1e10,2.44e7,2.5e7,3.02e7,3.12e7,3.99e7,4.15e7,5.85e7,6.22e7,1.1e8,1.23e8,9.1e8,1e10,2.44e7,2.5e7,3.02e7,3.12e7,3.99e7,4.15e7,5.85e7,6.22e7,1.1e8,1.23e8,9.1e8,1e10,2.44e7,2.5e7,3.02e7,3.12e7,3.99e7,4.15e7,5.85e7,6.22e7,1.1e8,1.23e8,9.1e8,1e10,2.44e7,2.5e7,3.02e7,3.12e7,3.99e7,4.15e7,5.85e7,6.22e7,1.1e8,1.23e8,9.1e8,1e10,2.44e7,2.5e7,3.02e7,3.12e7,3.99e7,4.15e7,5.85e7,6.22e7,1.1e8,1.23e8,9.1e8,1e10,2.44e7,2.5e7,3.02e7,3.12e7,3.99e7,4.15e7,5.85e7,6.22e7,1.1e8,1.23e8,9.1e8,1e10,2.44e7,2.5e7,3.02e7,3.12e7,3.99e7,4.15e7,5.85e7,6.22e7,1.1e8,1.23e8,9.1e8,1e10,2.44e7,2.5e7,3.02e7,3.12e7,3.99e7,4.15e7,5.85e7,6.22e7,1.1e8,1.23e8,9.1e8,1e10,2.44e7,2.5e7,3.02e7,3.12e7,3.99e7,4.15e7,5.85e7,6.22e7,1.1e8,1.23e8,9.1e8,1e10,2.44e7,2.5e7,3.02e7,3.12e7,3.99e7,4.15e7,5.85e7,6.22e7,1.1e8,1.23e8,9.1e8,1e10,2.44e7,2.5e7,3.02e7,3.12e7,3.99e7,4.15e7,5.85e7,6.22e7,1.1e8,1.23e8,9.1e8,1e10,2.44e7,2.5e7,3.02e7,3.12e7,3.99e7,4.15e7,5.85e7,6.22e7,1.1e8,1.23e8,9.1e8,1e10')
#res_range_val.set('2.4875e7,2.5125e7,2.6395e7,2.6928e7,2.8275e7,2.8846e7,3.0443e7,3.1059e7,3.2972e7,3.3638e7,3.5959e7,3.6685e7,3.9540e7,4.0339e7,4.3914e7,4.4801e7,4.9376e7,5.0374e7,5.6390e7,5.7529e7,6.5726e7,6.7053e7,7.8766e7,8.0358e7,9.8263e7,1.0024e8,1.3058e8,1.3322e8,1.9459e8,1.9852e8,3.8168e8,3.8939e8,9.9e8,1.01e10,2.4875e7,2.5125e7,2.6395e7,2.6928e7,2.8275e7,2.8846e7,3.0443e7,3.1059e7,3.2972e7,3.3638e7,3.5959e7,3.6685e7,3.9540e7,4.0339e7,4.3914e7,4.4801e7,4.9376e7,5.0374e7,5.6390e7,5.7529e7,6.5726e7,6.7053e7,7.8766e7,8.0358e7,9.8263e7,1.0024e8,1.3058e8,1.3322e8,1.9459e8,1.9852e8,3.8168e8,3.8939e8,9.9e8,1.01e10,2.4875e7,2.5125e7,2.6395e7,2.6928e7,2.8275e7,2.8846e7,3.0443e7,3.1059e7,3.2972e7,3.3638e7,3.5959e7,3.6685e7,3.9540e7,4.0339e7,4.3914e7,4.4801e7,4.9376e7,5.0374e7,5.6390e7,5.7529e7,6.5726e7,6.7053e7,7.8766e7,8.0358e7,9.8263e7,1.0024e8,1.3058e8,1.3322e8,1.9459e8,1.9852e8,3.8168e8,3.8939e8,9.9e8,1.01e10,2.4875e7,2.5125e7,2.6395e7,2.6928e7,2.8275e7,2.8846e7,3.0443e7,3.1059e7,3.2972e7,3.3638e7,3.5959e7,3.6685e7,3.9540e7,4.0339e7,4.3914e7,4.4801e7,4.9376e7,5.0374e7,5.6390e7,5.7529e7,6.5726e7,6.7053e7,7.8766e7,8.0358e7,9.8263e7,1.0024e8,1.3058e8,1.3322e8,1.9459e8,1.9852e8,3.8168e8,3.8939e8,9.9e8,1.01e10,2.4875e7,2.5125e7,2.6395e7,2.6928e7,2.8275e7,2.8846e7,3.0443e7,3.1059e7,3.2972e7,3.3638e7,3.5959e7,3.6685e7,3.9540e7,4.0339e7,4.3914e7,4.4801e7,4.9376e7,5.0374e7,5.6390e7,5.7529e7,6.5726e7,6.7053e7,7.8766e7,8.0358e7,9.8263e7,1.0024e8,1.3058e8,1.3322e8,1.9459e8,1.9852e8,3.8168e8,3.8939e8,9.9e8,1.01e10,2.4875e7,2.5125e7,2.6395e7,2.6928e7,2.8275e7,2.8846e7,3.0443e7,3.1059e7,3.2972e7,3.3638e7,3.5959e7,3.6685e7,3.9540e7,4.0339e7,4.3914e7,4.4801e7,4.9376e7,5.0374e7,5.6390e7,5.7529e7,6.5726e7,6.7053e7,7.8766e7,8.0358e7,9.8263e7,1.0024e8,1.3058e8,1.3322e8,1.9459e8,1.9852e8,3.8168e8,3.8939e8,9.9e8,1.01e10,2.4875e7,2.5125e7,2.6395e7,2.6928e7,2.8275e7,2.8846e7,3.0443e7,3.1059e7,3.2972e7,3.3638e7,3.5959e7,3.6685e7,3.9540e7,4.0339e7,4.3914e7,4.4801e7,4.9376e7,5.0374e7,5.6390e7,5.7529e7,6.5726e7,6.7053e7,7.8766e7,8.0358e7,9.8263e7,1.0024e8,1.3058e8,1.3322e8,1.9459e8,1.9852e8,3.8168e8,3.8939e8,9.9e8,1.01e10,2.4875e7,2.5125e7,2.6395e7,2.6928e7,2.8275e7,2.8846e7,3.0443e7,3.1059e7,3.2972e7,3.3638e7,3.5959e7,3.6685e7,3.9540e7,4.0339e7,4.3914e7,4.4801e7,4.9376e7,5.0374e7,5.6390e7,5.7529e7,6.5726e7,6.7053e7,7.8766e7,8.0358e7,9.8263e7,1.0024e8,1.3058e8,1.3322e8,1.9459e8,1.9852e8,3.8168e8,3.8939e8,9.9e8,1.01e10,2.4875e7,2.5125e7,2.6395e7,2.6928e7,2.8275e7,2.8846e7,3.0443e7,3.1059e7,3.2972e7,3.3638e7,3.5959e7,3.6685e7,3.9540e7,4.0339e7,4.3914e7,4.4801e7,4.9376e7,5.0374e7,5.6390e7,5.7529e7,6.5726e7,6.7053e7,7.8766e7,8.0358e7,9.8263e7,1.0024e8,1.3058e8,1.3322e8,1.9459e8,1.9852e8,3.8168e8,3.8939e8,9.9e8,1.01e10,2.4875e7,2.5125e7,2.6395e7,2.6928e7,2.8275e7,2.8846e7,3.0443e7,3.1059e7,3.2972e7,3.3638e7,3.5959e7,3.6685e7,3.9540e7,4.0339e7,4.3914e7,4.4801e7,4.9376e7,5.0374e7,5.6390e7,5.7529e7,6.5726e7,6.7053e7,7.8766e7,8.0358e7,9.8263e7,1.0024e8,1.3058e8,1.3322e8,1.9459e8,1.9852e8,3.8168e8,3.8939e8,9.9e8,1.01e10,2.4875e7,2.5125e7,2.6395e7,2.6928e7,2.8275e7,2.8846e7,3.0443e7,3.1059e7,3.2972e7,3.3638e7,3.5959e7,3.6685e7,3.9540e7,4.0339e7,4.3914e7,4.4801e7,4.9376e7,5.0374e7,5.6390e7,5.7529e7,6.5726e7,6.7053e7,7.8766e7,8.0358e7,9.8263e7,1.0024e8,1.3058e8,1.3322e8,1.9459e8,1.9852e8,3.8168e8,3.8939e8,9.9e8,1.01e10')
#res_range_val.set('2.375e7,2.625e7,2.713e7,3.000e7,3.164e7,3.497e7,3.794e7,4.193e7,4.738e7,5.236e7,6.307e7,6.970e7,9.429e7,1.042e8,1.867e8,2.063e8,9.5e9,1.05e10')
#res_range_val.set('2.375e7,2.625e7,2.375e7,2.625e7,2.713e7,3.000e7,2.713e7,3.000e7,3.164e7,3.497e7,3.164e7,3.497e7,3.794e7,4.193e7,3.794e7,4.193e7,4.738e7,5.236e7,4.738e7,5.236e7,6.307e7,6.970e7,6.307e7,6.970e7,9.429e7,1.042e8,9.429e7,1.042e8,1.867e8,2.063e8,9.5e9,1.05e10,9.5e9,1.05e10,2.375e7,2.625e7,2.375e7,2.625e7,2.713e7,3.000e7,2.713e7,3.000e7,3.164e7,3.497e7,3.164e7,3.497e7,3.794e7,4.193e7,3.794e7,4.193e7,4.738e7,5.236e7,4.738e7,5.236e7,6.307e7,6.970e7,6.307e7,6.970e7,9.429e7,1.042e8,9.429e7,1.042e8,1.867e8,2.063e8,9.5e9,1.05e10,9.5e9,1.05e10,2.375e7,2.625e7,2.375e7,2.625e7,2.713e7,3.000e7,2.713e7,3.000e7,3.164e7,3.497e7,3.164e7,3.497e7,3.794e7,4.193e7,3.794e7,4.193e7,4.738e7,5.236e7,4.738e7,5.236e7,6.307e7,6.970e7,6.307e7,6.970e7,9.429e7,1.042e8,9.429e7,1.042e8,1.867e8,2.063e8,9.5e9,1.05e10,9.5e9,1.05e10,2.375e7,2.625e7,2.375e7,2.625e7,2.713e7,3.000e7,2.713e7,3.000e7,3.164e7,3.497e7,3.164e7,3.497e7,3.794e7,4.193e7,3.794e7,4.193e7,4.738e7,5.236e7,4.738e7,5.236e7,6.307e7,6.970e7,6.307e7,6.970e7,9.429e7,1.042e8,9.429e7,1.042e8,1.867e8,2.063e8,9.5e9,1.05e10,9.5e9,1.05e10,2.375e7,2.625e7,2.375e7,2.625e7,2.713e7,3.000e7,2.713e7,3.000e7,3.164e7,3.497e7,3.164e7,3.497e7,3.794e7,4.193e7,3.794e7,4.193e7,4.738e7,5.236e7,4.738e7,5.236e7,6.307e7,6.970e7,6.307e7,6.970e7,9.429e7,1.042e8,9.429e7,1.042e8,1.867e8,2.063e8,9.5e9,1.05e10,9.5e9,1.05e10,2.375e7,2.625e7,2.375e7,2.625e7,2.713e7,3.000e7,2.713e7,3.000e7,3.164e7,3.497e7,3.164e7,3.497e7,3.794e7,4.193e7,3.794e7,4.193e7,4.738e7,5.236e7,4.738e7,5.236e7,6.307e7,6.970e7,6.307e7,6.970e7,9.429e7,1.042e8,9.429e7,1.042e8,1.867e8,2.063e8,9.5e9,1.05e10,9.5e9,1.05e10,2.375e7,2.625e7,2.375e7,2.625e7,2.713e7,3.000e7,2.713e7,3.000e7,3.164e7,3.497e7,3.164e7,3.497e7,3.794e7,4.193e7,3.794e7,4.193e7,4.738e7,5.236e7,4.738e7,5.236e7,6.307e7,6.970e7,6.307e7,6.970e7,9.429e7,1.042e8,9.429e7,1.042e8,1.867e8,2.063e8,9.5e9,1.05e10,9.5e9,1.05e10,2.375e7,2.625e7,2.375e7,2.625e7,2.713e7,3.000e7,2.713e7,3.000e7,3.164e7,3.497e7,3.164e7,3.497e7,3.794e7,4.193e7,3.794e7,4.193e7,4.738e7,5.236e7,4.738e7,5.236e7,6.307e7,6.970e7,6.307e7,6.970e7,9.429e7,1.042e8,9.429e7,1.042e8,1.867e8,2.063e8,9.5e9,1.05e10,9.5e9,1.05e10,2.375e7,2.625e7,2.375e7,2.625e7,2.713e7,3.000e7,2.713e7,3.000e7,3.164e7,3.497e7,3.164e7,3.497e7,3.794e7,4.193e7,3.794e7,4.193e7,4.738e7,5.236e7,4.738e7,5.236e7,6.307e7,6.970e7,6.307e7,6.970e7,9.429e7,1.042e8,9.429e7,1.042e8,1.867e8,2.063e8,9.5e9,1.05e10,9.5e9,1.05e10,2.375e7,2.625e7,2.375e7,2.625e7,2.713e7,3.000e7,2.713e7,3.000e7,3.164e7,3.497e7,3.164e7,3.497e7,3.794e7,4.193e7,3.794e7,4.193e7,4.738e7,5.236e7,4.738e7,5.236e7,6.307e7,6.970e7,6.307e7,6.970e7,9.429e7,1.042e8,9.429e7,1.042e8,1.867e8,2.063e8,9.5e9,1.05e10,9.5e9,1.05e10,2.375e7,2.625e7,2.375e7,2.625e7,2.713e7,3.000e7,2.713e7,3.000e7,3.164e7,3.497e7,3.164e7,3.497e7,3.794e7,4.193e7,3.794e7,4.193e7,4.738e7,5.236e7,4.738e7,5.236e7,6.307e7,6.970e7,6.307e7,6.970e7,9.429e7,1.042e8,9.429e7,1.042e8,1.867e8,2.063e8,9.5e9,1.05e10,9.5e9,1.05e10,2.375e7,2.625e7,2.375e7,2.625e7,2.713e7,3.000e7,2.713e7,3.000e7,3.164e7,3.497e7,3.164e7,3.497e7,3.794e7,4.193e7,3.794e7,4.193e7,4.738e7,5.236e7,4.738e7,5.236e7,6.307e7,6.970e7,6.307e7,6.970e7,9.429e7,1.042e8,9.429e7,1.042e8,1.867e8,2.063e8,9.5e9,1.05e10,9.5e9,1.05e10,2.375e7,2.625e7,2.375e7,2.625e7,2.713e7,3.000e7,2.713e7,3.000e7,3.164e7,3.497e7,3.164e7,3.497e7,3.794e7,4.193e7,3.794e7,4.193e7,4.738e7,5.236e7,4.738e7,5.236e7,6.307e7,6.970e7,6.307e7,6.970e7,9.429e7,1.042e8,9.429e7,1.042e8,1.867e8,2.063e8,9.5e9,1.05e10,9.5e9,1.05e10,2.375e7,2.625e7,2.375e7,2.625e7,2.713e7,3.000e7,2.713e7,3.000e7,3.164e7,3.497e7,3.164e7,3.497e7,3.794e7,4.193e7,3.794e7,4.193e7,4.738e7,5.236e7,4.738e7,5.236e7,6.307e7,6.970e7,6.307e7,6.970e7,9.429e7,1.042e8,9.429e7,1.042e8,1.867e8,2.063e8,9.5e9,1.05e10,9.5e9,1.05e10,2.375e7,2.625e7,2.375e7,2.625e7,2.713e7,3.000e7,2.713e7,3.000e7,3.164e7,3.497e7,3.164e7,3.497e7,3.794e7,4.193e7,3.794e7,4.193e7,4.738e7,5.236e7,4.738e7,5.236e7,6.307e7,6.970e7,6.307e7,6.970e7,9.429e7,1.042e8,9.429e7,1.042e8,1.867e8,2.063e8,9.5e9,1.05e10,9.5e9,1.05e10,2.375e7,2.625e7,2.375e7,2.625e7,2.713e7,3.000e7,2.713e7,3.000e7,3.164e7,3.497e7,3.164e7,3.497e7,3.794e7,4.193e7,3.794e7,4.193e7,4.738e7,5.236e7,4.738e7,5.236e7,6.307e7,6.970e7,6.307e7,6.970e7,9.429e7,1.042e8,9.429e7,1.042e8,1.867e8,2.063e8,9.5e9,1.05e10,9.5e9,1.05e10,2.375e7,2.625e7,2.375e7,2.625e7,2.713e7,3.000e7,2.713e7,3.000e7,3.164e7,3.497e7,3.164e7,3.497e7,3.794e7,4.193e7,3.794e7,4.193e7,4.738e7,5.236e7,4.738e7,5.236e7,6.307e7,6.970e7,6.307e7,6.970e7,9.429e7,1.042e8,9.429e7,1.042e8,1.867e8,2.063e8,9.5e9,1.05e10,9.5e9,1.05e10,2.375e7,2.625e7,2.375e7,2.625e7,2.713e7,3.000e7,2.713e7,3.000e7,3.164e7,3.497e7,3.164e7,3.497e7,3.794e7,4.193e7,3.794e7,4.193e7,4.738e7,5.236e7,4.738e7,5.236e7,6.307e7,6.970e7,6.307e7,6.970e7,9.429e7,1.042e8,9.429e7,1.042e8,1.867e8,2.063e8,9.5e9,1.05e10,9.5e9,1.05e10,2.375e7,2.625e7,2.375e7,2.625e7,2.713e7,3.000e7,2.713e7,3.000e7,3.164e7,3.497e7,3.164e7,3.497e7,3.794e7,4.193e7,3.794e7,4.193e7,4.738e7,5.236e7,4.738e7,5.236e7,6.307e7,6.970e7,6.307e7,6.970e7,9.429e7,1.042e8,9.429e7,1.042e8,1.867e8,2.063e8,9.5e9,1.05e10,9.5e9,1.05e10,2.375e7,2.625e7,2.375e7,2.625e7,2.713e7,3.000e7,2.713e7,3.000e7,3.164e7,3.497e7,3.164e7,3.497e7,3.794e7,4.193e7,3.794e7,4.193e7,4.738e7,5.236e7,4.738e7,5.236e7,6.307e7,6.970e7,6.307e7,6.970e7,9.429e7,1.042e8,9.429e7,1.042e8,1.867e8,2.063e8,9.5e9,1.05e10,9.5e9,1.05e10')
res_range_val.set('7.877e7,9.296e7,1.198e8,1.441e8')
# 
# # Typical high resistance range: 5e8,1e10
#Typical low resistance range: 2.38e7,2.5e7

#'0.9e8,4.77e8,4.77e8,1e10,3.8e7,4.53e7,4.53e7,4.98e7,2.2e7,2.38e7,2.38e7,2.5e7,1.5e7,1.61e7,1.61e7,1.67e7,1.18e7,1.22e7,1.22e7,1.25e7,0.99e7,1.02e7,1.02e7,1.04e7'
#4.77e8,1e10,4.53e7,4.98e7,2.38e7,2.5e7,1.61e7,1.67e7,1.22e7,1.25e7,1.02e7,1.04e7
#Dev1-R2-B_C4-1

# Tries Limit
max_tries_label = tk.Label(dynamic_frame, text = 'Max Tries')
max_tries_label.grid(column=0, row=2, sticky=tk.W, padx=5, pady=5)
max_tries_val = tk.IntVar()
max_tries_meas = tk.Entry(dynamic_frame, text=max_tries_val, width=5)
max_tries_meas.grid(column=1, row=2, sticky=tk.W, padx=5, pady=5)
max_tries_meas.grid_columnconfigure(1, weight=1)
max_tries_val.set(200)

# Follow up tests
follow_up = tk.BooleanVar()
follow_up_check = tk.Checkbutton(dynamic_frame, text='Retention Sweep', variable=follow_up)
follow_up.set(True)
follow_up_check.grid(column=0, row=3, sticky=tk.W, padx=5, pady=5)

follow_time_label = tk.Label(dynamic_frame, text='Added Sampling Time (s)')
follow_time_label.grid(column=0, row=4, sticky=tk.W, padx=5, pady=5)
follow_time_val = tk.IntVar()
follow_time_meas = tk.Entry(dynamic_frame, text=follow_time_val, width=5)
follow_time_meas.grid(column=1, row=4, sticky=tk.W, padx=5, pady=5)
follow_time_meas.grid_columnconfigure(1, weight=1)
follow_time_val.set(0)

follow_num_label = tk.Label(dynamic_frame, text='Samples')
follow_num_label.grid(column=0, row=5, sticky=tk.W, padx=5, pady=5)
follow_num_val = tk.IntVar()
follow_num_meas = tk.Entry(dynamic_frame, text=follow_num_val, width=5)
follow_num_meas.grid(column=1, row=5, sticky=tk.W, padx=5, pady=5)
follow_num_meas.grid_columnconfigure(1, weight=1)
follow_num_val.set(24000)

start_delay_label = tk.Label(dynamic_frame, text="Start Delay (s)")
start_delay_label.grid(column=0, row=6, sticky='w', padx=5, pady=5)
start_delay_val = tk.IntVar()
start_delay_meas = tk.Entry(dynamic_frame, text=start_delay_val, width=5)
start_delay_meas.grid(column=1, row=6, sticky='w', padx=5, pady=5)
start_delay_meas.grid_columnconfigure(1, weight=1)
start_delay_val.set(1)

auto_range_on = tk.BooleanVar()
auto_range_check = tk.Checkbutton(dynamic_frame, text="Auto Range", variable=auto_range_on)
auto_range_check.grid(column=0, row=7, sticky='w', padx=5, pady=5)
auto_range_on.set(False)

## Plot Settings
plot_frame = tk.Frame(dashboard_frame, bd=2, width=375)
plot_frame.grid(column=1, row=2, sticky=tk.W+tk.E)
plot_frame.grid_columnconfigure(0, weight=1)
# Button to start sweep
ft = font.Font(size='14', weight='bold')
save_dir_btn = tk.Button(plot_frame, text='Target Resistance Sweep', command=dynamic_scan_call, font=ft)     # THIS STARTS THE SCAN!!
save_dir_btn.grid(column=0, row=0, columnspan=3, sticky='EW', padx=5, pady=5)
#Button for single pulse
ft = font.Font(size='14', weight='bold')
save_dir_btn = tk.Button(plot_frame, text='Start Pulsing', command=single_pulse, font=ft)
save_dir_btn.grid(column=0, row=1, columnspan=3, sticky='EW', padx=5, pady=5)
# Follow loop only?
ft = font.Font(size='14', weight='bold')
# follow_time = follow_time_val.get()
#         follow_num = follow_num_val.get()

save_dir_btn = tk.Button(plot_frame, text='Follow Loop', 
                         command=follow_loop, font=ft)
save_dir_btn.grid(column=0, row=2, columnspan=3, sticky='EW', padx=5, pady=5)
# Progress bar for the scan
progress_var = tk.DoubleVar()
pb = ttk.Progressbar(plot_frame, orient='horizontal', mode='determinate', variable=progress_var)
# Labels for the current sens and current
sens_var = tk.StringVar()
curr_var = tk.StringVar()
sens_label = tk.Label(plot_frame, textvariable=sens_var)
curr_label = tk.Label(plot_frame, textvariable=curr_var)
def set_e_stop():
    global e_stop
    e_stop = True
stop_btn = tk.Button(plot_frame, text='Stop', command=set_e_stop, 
                    bg='red', fg='white', font=ft)

# Keep previous plots
hold_on = tk.BooleanVar()
plot_hold_on_check = tk.Checkbutton(plot_frame, text='Hold On', variable=hold_on)
hold_on.set(True)
plot_hold_on_check.grid(column=0, row=4, padx=5, pady=5)
# Set linthresh for symlog
plot_linthresh_label = tk.Label(plot_frame, text='Symlog linthresh')
plot_linthresh_val = tk.DoubleVar()
plot_linthresh_meas = tk.Entry(plot_frame, text=plot_linthresh_val, width=10)
plot_linthresh_val.set(1e-9)
# Plot Scale
plot_scale_label = tk.Label(plot_frame, text='Y-Scale')
plot_scale_label.grid(column=1, row=4, sticky=tk.W, padx=5, pady=5)
plot_scale_options = ['linear', 'log', 'symlog']
plot_scale_clicked = tk.StringVar()
plot_scale_clicked.set(plot_scale_options[0])
plot_scale_menu = tk.OptionMenu(plot_frame, plot_scale_clicked ,*plot_scale_options)
plot_scale_menu.config(width=6)
plot_scale_menu.grid(column=2, row=4, sticky=tk.E, padx=5, pady=5)

plot_linthresh_val.trace('w', scale_change)
plot_scale_clicked.trace('w', scale_change)

## Figure settings
font = {'size'   : 10}
rc('font', **font)

# Measurement Figures
figure_frame = tk.Frame(root)
figure_frame.grid(column=1, row=0, sticky='e')


fig_meas, ax_meas = plt.subplots(1, 1)
fig_meas.set_size_inches(8, 3.5)
canvas_meas = mplTk.FigureCanvasTkAgg(fig_meas, master=figure_frame)
canvas_meas.get_tk_widget().grid(column=0, row=0, sticky='e')
ax_meas.set_xlabel('Measurement #')
ax_meas.set_ylabel('Current (A)')
fig_meas.tight_layout()


frame_meas_toolbar = tk.Frame(figure_frame)
frame_meas_toolbar.grid(row=1, column=0, sticky=tk.E)
toolbar_meas = mplTk.NavigationToolbar2Tk(canvas_meas, frame_meas_toolbar)
toolbar_meas.update()

# Pulse figure
fig_in, ax_in = plt.subplots(1, 1)
fig_in.set_size_inches(8, 3.5)
canvas_in = mplTk.FigureCanvasTkAgg(fig_in, master=figure_frame)
canvas_in.get_tk_widget().grid(column=0, row=2, sticky=tk.W+tk.E)
ax_in.set_xlabel('Measurement #')
ax_in.set_ylabel('Pulse Voltage (V)')
fig_in.tight_layout()

frame_in_toolbar = tk.Frame(figure_frame)
frame_in_toolbar.grid(row=3, column=0, sticky=tk.E)
toolbar_in = mplTk.NavigationToolbar2Tk(canvas_in, frame_in_toolbar)
toolbar_in.update()

root.mainloop()
