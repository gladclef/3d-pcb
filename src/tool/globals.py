from tool.units import *


class _Globals:
    _instance: "_Globals" = None

    def __init__(self):
        print("Setting _Globals._instance")
        self._instance = self

        #####################
        # measurable values #
        #####################
        
        self.BREADBOARD_SPACING = in2mm(0.1)
        """
        The standard pin spacing for breadboards.
        Limits how large traces can be. Defaults to 0.1 inches.
        """
        self.NOZZLE_DIAMETER = 0.4
        """
        Diameter of the 3d printer's nozzle. Limits how small the space between
        traces can be. Defaults to 0.4 millimeters.
        """
        self.LAYER_HEIGHT = 0.1
        """
        How tall we expect the sliced layers to be. Affects how long the
        opening to the traces can be. Defaults to 0.1 millimeters.
        """
        self.BOARD_THICKNESS = 2.0
        """
        Expected thickness of the printed circuit board. Determines the
        location of bottom layer traces. Defaults to 2 millimeters.
        """
        self.VIA_DIAMETER = 0.4
        """
        Size of the drill holes used for trace vias.
        Default is 0.4 millimeters.
        """
        self.THROUGH_HOLE_DIAMETER = 0.4
        """
        Size of the drill holes used for through-hole mounted components.
        Default is 0.4 millimeters.
        """

        ##################
        # tunable values #
        ##################

        self.TRACE_OPENING_CLEARANCE = 0.2
        """
        Additional value added to (or subtraced from) the trace opening, which
        is otherwise based on the WIRE_DIAMETER. Increasing this value will
        make it easier to place wires into the circuit board.
        """
        self.TRACE_DIAMETER_CLEARANCE = 0.1
        """
        Additional value added to (or subtraced from) the trace diameter, which
        is otherwise based on the WIRE_DIAMETER. Increasing this value will
        make it more likely that the wire will comfortably sit in its trace and
        slide through the trace more easily.
        """
        self.VIA_DIAMETER_CLEARANCE = 0.1
        """
        Additional value added to (or subtraced from) the via diameter.
        Increasing this value will make it easier to fit wires through the vias.
        """
        self.THROUGH_HOLE_DIAMETER_CLEARANCE = 0.1
        """
        Additional value added to (or subtraced from) the drill hole diameter.
        Increasing this value will make it easier to fit through-hole components
        through the drill holes.
        """
    
    @classmethod
    def get_instance(cls) -> "_Globals":
        if cls._instance is not None:
            return cls._instance
        return cls()

board_parameters = _Globals.get_instance()