import serial
import serial.tools.list_ports
import time
from time import sleep


class handlerCore:
    def __init__(self, handshakes=2):
        self.valvestate = 0
        self.pumpstate = 0
        self.titanxstates = [0 for i in range(5)]
        self.handshakes = handshakes
        self.titanx_states = {
            "Probe 1": [0, 2, 1, 12, 0],
            "Probe 2": [0, 3, 1, 12, 0],
            "Probe 3": [0, 4, 1, 12, 0],
            "Probe 4": [0, 5, 1, 12, 0],
            "Probe 5": [0, 6, 1, 12, 0],
            "Probe 6": [0, 7, 1, 12, 0],
            "Probe 7": [0, 8, 1, 12, 0],
            "Probe 8": [0, 9, 1, 12, 0],
            "Probe 9": [0, 10, 1, 12, 0],
            "Probe 10": [0, 11, 1, 12, 0],
            "Probe 11": [0, 12, 1, 12, 0],
            "Probe 12": [0, 1, 2, 12, 0],
            "Probe 13": [0, 1, 3, 12, 0],
            "Probe 14": [0, 1, 4, 12, 0],
            "Probe 15": [0, 1, 5, 12, 0],
            "Probe 16": [0, 1, 6, 12, 0],
            "Probe 17": [0, 1, 7, 12, 0],
            "Probe 18": [0, 1, 8, 12, 0],
            "Probe 19": [0, 1, 9, 12, 0],
            "Probe 20": [0, 1, 10, 12, 0],
            "Probe 21": [0, 1, 11, 12, 0],
            "Probe 22": [0, 1, 12, 12, 0],
            "Probe 23": [0, 1, 1, 1, 0],
            "Probe 24": [0, 1, 1, 2, 0],
            "SSC": [0, 1, 1, 8, 0],
            "PFA(half-MeAc)": [0, 1, 1, 9, 0],
            "EtOH(MeAc)": [0, 1, 1, 10, 0],
            "Image": [0, 1, 1, 11, 0],
            "Cleave": [0, 1, 1, 12, 0],
        }

    def get_heartbeat(self, comport, connect_code="MARLIN", timeout=10.0):
        try:
            ti = time.time()
            no_timeout = True
            s = serial.Serial(comport, 9600, timeout=0.5)
            readcmd = "5\n".encode("ascii")
            while no_timeout:
                s.write(readcmd)
                returnedstr = s.read_until()
                t_elapsed = time.time() - ti
                if len(returnedstr) > 0:
                    no_timeout = False
                elif t_elapsed > timeout:
                    no_timeout = False
            s.close()
            if returnedstr == "MARLIN":
                return True
            else:
                return False
        except (OSError, serial.SerialException):
            return False

    def connect(self, connect_code="MARLIN", timeout=10.0):
        ports = ["COM%s" % (i + 1) for i in range(256)]
        result = []
        for port in ports:
            heartbeat = self.get_heartbeat(
                port, connect_code=connect_code, timeout=timeout
            )
            if heartbeat:
                result.append(port)
        if len(result) == 0:
            raise ValueError("No MARLIN detected.")
        elif len(result) == 1:
            ti = time.time()
            no_timeout = True
            self.serial_handle = serial.Serial(result[0], 9600, timeout=0.5)
            readcmd = "5\n".encode("ascii")
            while no_timeout:
                self.serial_handle.write(readcmd)
                returnedstr = self.serial_handle.read_until()
                t_elapsed = time.time() - ti
                if len(returnedstr) > 0:
                    no_timeout = False
                elif t_elapsed > timeout:
                    no_timeout = False
            if returnedstr == "MARLIN":
                self.serial_handle.timeout = 10.0
                print ("Connected.")
            else:
                raise ValueError("MARLIN connection timeout.")
        else:
            raise ValueError("More than one MARLIN detected.")

    def updatestate(self, valvestate, pumpstate, titanxstates):
        self.valvestate = valvestate
        self.pumpstate = pumpstate
        self.titanxstates = titanxstates

    def sendstate(self, valvestate, pumpstate, titanxstates):

        self.updatestate(valvestate, pumpstate, titanxstates)
        no_handshake = True
        handshake_failed = False
        handshake_attempts = 0

        while no_handshake:
            valvestr = "4" + str(self.valvestate)
            pumpstr = str(self.pumpstate)
            pumpstr = "3" + ("0" * (4 - len(pumpstr)) + pumpstr)

            titanxstrlist = []

            for titannum, titanxstate in enumerate(self.titanxstates):
                if titanxstate != 0:
                    titanxstr = str(titanxstate)
                    titanxstr = (
                        "2" + str(titannum) + ("0" * (2 - len(titanxstr)) + titanxstr)
                    )
                    titanxstrlist.append(titanxstr)

            cmdlist = [""] + [valvestr] + [pumpstr] + titanxstrlist

            for cmd in cmdlist:
                sendstr = cmd + "\n"
                statestr = sendstr.encode("ascii")
                self.serial_handle.write(statestr)
                time.sleep(0.25)

            readcmd = "0\n".encode("ascii")
            self.serial_handle.write(readcmd)
            returnedstr = self.serial_handle.read_until()[:-1]
            self.serial_handle.reset_output_buffer()
            self.serial_handle.reset_input_buffer()
            checkstr = (
                "["
                + ",".join([str(state) for state in self.titanxstates])
                + "];"
                + str(self.valvestate)
                + ";"
                + str(self.pumpstate)
            )

            print returnedstr.strip()

            if returnedstr.strip() == checkstr.strip():
                no_handshake = False
            handshake_attempts += 1
            if handshake_attempts >= self.handshakes:
                raise Exception("Handshake failed.")

    def set_valve_state(self, titanx_state_name, valvestate):
        titanxstates = self.titanx_states[titanx_state_name]
        self.sendstate(valvestate, self.pumpstate, titanxstates)

    def set_pump_state(self, pumpstate):
        self.sendstate(self.valvestate, pumpstate, self.titanxstates)
