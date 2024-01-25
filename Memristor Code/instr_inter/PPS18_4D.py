import pyvisa as visa
from time import sleep

class PPS18_4D:
    def __init__(self, name='GPIB0::12::INSTR') -> None:
        rm = visa.ResourceManager()
        if not name in rm.list_resources(name):
            raise visa.VisaIOError(1) #randomly chosen error code
        self.dc_ps = rm.open_resource(name)
        self.dc_ps.write('OCP 1') # Turns on overcurrent protection
        sleep(.1)
        self.dc_ps.write('VOFF 0') # Sets the offset voltage to 0
        sleep(.1)
        self.dc_ps.write('VSET 0') # Set initial voltage to 0
        sleep(.1)
        self.dc_ps.write('OUT 1') # Enables output on channel 1
        sleep(.1)
    
    def set_voltage(self, v):
        self.dc_ps.write(f'VSET {v}')
        sleep(.25)