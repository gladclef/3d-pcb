from abc import ABC, abstractmethod

import vtk

from Geometry.Point import Point
from Trace.AbstractVtkPointTracker import AbstractVtkPointTracker as PntInc
from FileIO.Line import Line as FLine
from Geometry.Path import Path
from Trace.PipeShape import PipeShape, DEFAULT_PIPE_SHAPE
from Trace.TraceLine import TraceLine
from tool.globals import board_parameters as g
from tool.units import *

class AbstractTrace(Path, PntInc, ABC):
    """
    Represents a wire trace on a PCB.
    """
    def __init__(self,
                 source_lines: list[FLine],
                 segments: list[TraceLine],
                 shape: PipeShape=None):
        Path.__init__(self, source_lines, segments)

        # set some defaults
        if shape is None:
            shape = DEFAULT_PIPE_SHAPE(g.DEFAULT_WIRE_DIAMETER)

        self.shape = shape
        self.segments: list[TraceLine] = self.segments

    def segments_at_xypnt(self, pnt: Point) -> list[TraceLine]:
        return super().segments_at_xypnt(pnt)

    def insert_xypnt(self, new_xy_pnt: Point, old_segment: TraceLine) -> tuple[TraceLine, TraceLine]:
        prev_segment, next_segment = super().insert_xypnt(new_xy_pnt, old_segment)
        
        new_prev_segment = TraceLine(prev_segment.fline, prev_segment.xy0, prev_segment.xy1)
        new_next_segment = TraceLine(next_segment.fline, next_segment.xy0, next_segment.xy1)
        new_prev_segment.copy_properties(old_segment)
        new_next_segment.copy_properties(old_segment)

        return new_prev_segment, new_next_segment
    
    def append_xypnt(self, new_xy_pnt: Point, at_start: bool, fline: FLine = None) -> TraceLine:
        segment = super().append_xypnt(new_xy_pnt, at_start, fline)

        new_segment = TraceLine(segment.fline, segment.xy0, segment.xy1)
        new_segment.is_end = True

        if 1 in new_segment.joined_ends:
            self.segments[1].is_end = False
            self.segments[1].is_inner_edge = False
        else:
            self.segments[-2].is_end = False
            self.segments[-2].is_inner_edge = False

        return new_segment

    def remove_xypnt(self, xy_pnt: Point) -> tuple[list[TraceLine], list[TraceLine]]:
        old_segments, new_segments = super().remove_xypnt(xy_pnt)

        new_new_segments: list[TraceLine] = []
        for segment in new_segments:
            new_new_segment = TraceLine(segment.fline, segment.xy0, segment.xy1)
            new_new_segment.copy_properties(old_segments[0])
            new_new_segments.append(new_new_segment)
        
        return old_segments, new_segments

    @abstractmethod
    def to_vtk(self) -> vtk.vtkPolyData:
        raise NotImplementedError