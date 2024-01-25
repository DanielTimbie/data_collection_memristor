import pyvisa as visa
from time import sleep

class ESP301:
    def __init__(self,
                 name='GPIB0::1::INSTR'):
        self.rm = visa.ResourceManager()

        self.motor = self.rm.open_resource(name)
        # self.motor.term_chars = '\n'
        # self.motor.baud_rate = 921600
        # self.motor.data_bits = 8
        # self.motor.stop_bits = visa.constants.StopBits.one
        # self.motor.parity = visa.constants.Parity.none


        # try:
        #     if not name in self.rm.list_resources(name):
        #         raise visa.VisaIOError(1) #randomly chosen error code
        #     DGen=self.rm.open_resource(name)

        # except visa.VisaIOError:
        #     raise Exception("Not attached to delay generator")

    def move_x_axis(self,distance):
        self.motor.write("2PR%f"%distance)
        sleep(1)

    def move_y_axis(self,distance):
        self.motor.write("1PR%f"%distance)
        sleep(1)

    def move_z_axis(self,distance):
        self.motor.write("3PR%f"%distance)
        sleep(1)

# esp = ESP301()
# esp.move_x_axis(-0.3)