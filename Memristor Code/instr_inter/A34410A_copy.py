import pyvisa as visa
from time import sleep
from time import perf_counter as timer
from warnings import warn

class A34410A:
    integ_times = [.006, .02, 0.06, .2, 1, 2, 10, 100]

    # integration time that stays below the maximum bandwidth of the
    # associated sensitivity in sens_table over at SR570.py 
    integ_times_table = [10, 10, 10, 10, 10, 10, 1, 1, 1, 1, 1, 1, .06, .06, .06,
                         .006, .006, .006, .006, .006, .006, .006, .006, .006,
                          .006, .006, .006, .006]

    def __init__(self, name='GPIB0::22::INSTR'): 
        try:
            rm = visa.ResourceManager()
            if not name in rm.list_resources(name):
                raise visa.VisaIOError(1) #randomly chosen error code
            self.DMM = rm.open_resource(name)  ## Digital Multi Meter Agilent 34410A
            self.DMM.write('*RST') 
            sleep(.1)
            self.DMM.write('CONF:VOLT:DC') 
            sleep(.1)
            self.DMM.write('VOLT:DC:NPLC 10') 
            sleep(.1)
            self.DMM.write('VOLT:DC:IMP:AUTO 1')
            sleep(.1)
            #self.DMM.write('VOLT::DC:ZERO:AUTO ONCE')

        except visa.VisaIOError:
            raise Exception("Not attached to Digital Multimeter")
    
    def set_dc_current_mode(self):
        self.DMM.write('CONF:CURR:DC') 
        sleep(.1)

    def read(self):
        self.DMM.write('READ?')
        return float(self.DMM.read())
    
    def set_integ_time(self,t=1):
        assert(t in self.integ_times)
        self.DMM.write(f'VOLT:DC:NPLC {t}')

    def write(self, command):
        self.DMM.write(command)
        # return float(self.DMM.read())

    def get(self):
        return self.DMM.read()

    def delay_sample(self, sample_count=100, trig_delay=0, sleep_time=1):
        self.DMM.write(f'SAMP:COUN {sample_count}') # take 100 samples
        self.DMM.write('SAMP:SOUR IMM') 
        self.DMM.write('TRIG:SOUR BUS')
        self.DMM.write(f'TRIG:DEL {trig_delay}')
        self.DMM.write('INIT')
        self.DMM.write('*TRG')
        self.DMM.write('FETC?')

        sleep(sleep_time)
        return self.DMM.read()
    
    def timer_sample(self, sample_count=100, trig_delay=0, source_timer=.1):
        self.write(f'SAMP:COUN {sample_count}') # take 100 samples
        sleep(.1)
        self.write('SAMP:SOUR TIM')
        sleep(.1) 
        self.write('TRIG:SOUR BUS')
        sleep(.1)
        self.write('VOLT:DC:RANG:AUTO ONCE')
        sleep(.1)
        self.write(f'TRIG:DEL {trig_delay}')
        sleep(.1)
        self.write(f'SAMP:TIM {source_timer}')
        sleep(.1)
        self.write('SAMP:TIM?')
        print("Sample time set: ", self.get())
        sleep(.1)
        
        self.write('STAT:QUES?')
        stat = int(self.get())
        print("Status:", stat) 
        if (stat & 4 == 4):
            warn("Sample timing violation")
        
        sleep(.1)
        self.write('INIT')
        sleep(.1)
        self.write('*TRG')
        start_time = timer()

        # this command does not clear the condition register upon reading
        self.write('STAT:OPER:COND?') 
        cond = self.get()
        try:
            cond = int(cond)
        except:
            raise Exception("Can't convert into integer.")
        finally:
            print("Condition start:", cond)

        # loop while multimeter is measuring
        # meaning 16 & cond == 16
        while(16 & cond == 16):
            sleep(.1)
            self.write('STAT:OPER:COND?')
            cond = self.get()
            try:
                cond = int(cond)
            except:
                raise Exception("Can't convert into integer.")

        print("Condition end:", cond)
        # sleep(trig_delay + sample_count*source_timer + sleep_time)
        end_time = timer()
        elapsed_time = end_time - start_time

        self.write('FETC?')
        sleep(.1)
        return self.get(), elapsed_time
        
        