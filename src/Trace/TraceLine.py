from FileIO.Line import Line as FLine
from Geometry.LineSegment import LineSegment
from Geometry.Point import Point


class TraceLine(LineSegment):
    def __init__(self, fline: FLine, xy0: Point, xy1: Point):
        super().__init__(xy0, xy1, fline)

        self._is_end: bool = False
        """ True if this edge only connects to one other edge.
        Mutually exclusive with is_solo and is_inner_edge. """
        self._is_solo: bool = False
        """ True if this edge doesn't connect to any other edges.
        Mutually exclusive with is_end and is_inner_edge. """
        self._is_inner_edge: bool = False
        """ True if this edge connects to another edge on both ends.
        Mutually exclusive with is_end and is_solo. """
        self.joined_ends: list[int] = []
        """ If this contains a 0, then self.xy0 connects to another edge.
        If this contains a 1, then self.xy1 connects to another edge. """
        self.is_branch: bool = False
        """ True if this edge is an end that branches off of another SingleTrace. """

    @property
    def is_end(self):
        return self._is_end
    
    @is_end.setter
    def is_end(self, val: bool):
        if val:
            assert (not self._is_solo) and (not self._is_inner_edge)
        self._is_end = val
    
    @property
    def is_solo(self):
        return self._is_solo
    
    @is_solo.setter
    def is_solo(self, val: bool):
        if val:
            assert (not self._is_end) and (not self._is_inner_edge)
        self._is_solo = val

    @property
    def is_inner_edge(self):
        return self._is_inner_edge
    
    @is_inner_edge.setter
    def is_inner_edge(self, val: bool):
        if val:
            assert (not self._is_solo) and (not self._is_solo)
        self._is_inner_edge = val

    @property
    def xy_points(self) -> tuple[Point, Point]:
        return (self.xy0, self.xy1)

    @property
    def joined_points(self) -> list[Point]:
        """ Return the points that correspond to any indexes in self.joined_ends """
        return [self.xy_points[i] for i in self.joined_ends]

    @property
    def unjoined_ends(self) -> list[int]:
        """ If this contains a 0, then self.xy0 does not connect to another edge.
        If this contains a 1, then self.xy1 does not connect to another edge. """
        return list(filter(lambda idx: idx not in self.joined_ends, [0, 1]))

    @property
    def unjoined_points(self) -> list[Point]:
        """ Return the points that correspond to any indexes in self.unjoined_ends """
        return [self.xy_points[i] for i in self.unjoined_ends]

    @classmethod
    def uncategorized(cls, trace_lines: list["TraceLine"]) -> list["TraceLine"]:
        """ Selects the edges that have not been classified as an end, solo, or inner edge. """
        return list(filter(lambda tl: not tl.is_end and not tl.is_solo and not tl.is_inner_edge, trace_lines))

    @classmethod
    def only_ends(cls, trace_lines: list["TraceLine"]) -> list["TraceLine"]:
        return list(filter(lambda tl: tl.is_end, trace_lines))

    @classmethod
    def only_solos(cls, trace_lines: list["TraceLine"]) -> list["TraceLine"]:
        return list(filter(lambda tl: tl.is_solo, trace_lines))

    @classmethod
    def only_inners(cls, trace_lines: list["TraceLine"]) -> list["TraceLine"]:
        return list(filter(lambda tl: tl.is_inner_edge, trace_lines))

    @classmethod
    def not_solo(cls, trace_lines: list["TraceLine"]) -> list["TraceLine"]:
        return list(filter(lambda tl: not tl.is_solo, trace_lines))
    
    def copy_properties(self, copy_from: "TraceLine"):
        self.is_end, self.is_solo, self.is_inner_edge = False, False, False
        self.is_end = copy_from.is_end
        self.is_solo = copy_from.is_solo
        self.is_inner_edge = copy_from.is_inner_edge
        self.joined_ends = copy_from.unjoined_ends
        self.is_branch = copy_from.is_branch

    def reversed(self) -> "TraceLine":
        ret = TraceLine(self.fline, self.xy1, self.xy0)
        ret.copy_properties(self)
        return ret

    def __repr__(self):
        return f"TraceLine<{self.xy0},{self.xy1},{self.fline}>"