import copy
from typing import TYPE_CHECKING

from FileIO.Line import Line as FLine
from Geometry.LineSegment import LineSegment


class Path:
    """
    Represents a 2D path in the XY plane.

    A path consists of back-to-back line segments.
    """

    def __init__(self,
                 source_lines: list[FLine],
                 xy_points: list[tuple[float, float]],
                 segments: list[tuple[int, int] | tuple[int, int, FLine]] | list[LineSegment]
    ):
        # sanity check/normalize input
        assert isinstance(source_lines, list) and isinstance(source_lines[0], FLine)
        new_xy_points = []
        for xy_idx, xy in enumerate(xy_points):
            txy = xy if isinstance(xy, tuple) else tuple(xy)
            try:
                x, y = txy
            except:
                raise ValueError(f"Error in {self.__class__.__name__}(): " + f"expected xy_points to be a list of 2 tuples, but value at xy_points[{xy_idx}]={xy}")
            new_xy_points.append(txy)

        self.source_lines = source_lines

        # build the list of segments the first time, to get the ordered points
        xypntindicies_2_segments = self._build_segments(xy_points, segments)
        new_segments = [s for i, s in xypntindicies_2_segments]
        ordered_xy_points = self._get_xy_points_in_segment_order(new_segments)

        # build the list of segments a second time, now that the points can be ordered
        new_segments_idxs = []
        for s in new_segments:
            idx1 = ordered_xy_points.index(s.xy1)
            idx2 = ordered_xy_points.index(s.xy2)
            new_segments_idxs.append(tuple(sorted([idx1, idx2])))
        self._xy_points = ordered_xy_points
        self._xypntindicies_2_segments = self._build_segments(ordered_xy_points, new_segments_idxs)
    
    @property
    def xy_points(self) -> list[tuple[float, float]]:
        """
        Returns the list of XY points that define this path.
        """
        return self._xy_points

    @property
    def edges(self) -> list[tuple[int, int]]:
        return [tuple([self.xy_points.index(p) for p in self.segment_xypnts(s)]) for s in self.segments]

    @property
    def segments(self):
        """
        Returns a list of LineSegment objects that make up this Path.
        """
        return [s[1] for s in self._xypntindicies_2_segments]
    
    def segments_at_xypnt(self, pnt_or_pntidx: tuple[float, float] | int) -> list[LineSegment]:

        """
        Returns a list of LineSegments that contain the given XY point.

        Args:
            pnt_or_pntidx (tuple[float, float] | int): Either an XY coordinate or its index in self._xy_points.

        Returns:
            list[LineSegment]: The segments that include the specified point.
        """
        try:
            x, y = pnt_or_pntidx
            pnt_idx = self.xy_points.index(pnt_or_pntidx)
        except:
            pnt_idx: int = pnt_or_pntidx
        
        matching_segments = filter(lambda s: pnt_idx in s[0], self._xypntindicies_2_segments)
        return [s[1] for s in matching_segments]
    
    def segment_xypnts(self, segment: LineSegment) -> tuple[tuple[float, float], tuple[float, float]]:
        """
        Returns the XY points of a given LineSegment.

        Args:
            segment (LineSegment): The LineSegment to find the endpoints for.

        Returns:
            tuple[tuple[float, float], tuple[float, float]]: The two XY coordinates that make up the segment.
        """
        for s in self._xypntindicies_2_segments:
            if s[1] == segment:
                return s[0]
    
    @staticmethod
    def _build_segment(xy_points: list[tuple[float, float]], segment: tuple[int, int] | tuple[int, int, FLine], source_line: FLine = None) -> LineSegment:
        """
        Builds a LineSegment from two XY points.

        Params
        ------
        xy_points: list[tuple[float, float]]
            List of x,y pairs.
        segment: tuple[int, int]
            The indices of the start and end point indicies in xy_points.

        Returns
        -------
        segment: LineSegment
            A new LineSegment object.
        """
        x1, y1 = xy_points[segment[0]]
        x2, y2 = xy_points[segment[1]]
        source_line = None if len(segment) < 3 else segment[2]
        return LineSegment((x1, y1), (x2, y2), source_line=source_line)

    @classmethod
    def _build_segments(cls,
                        xy_points: list[tuple[float, float]], 
                        segments: list[tuple[int, int] | tuple[int, int, FLine]] | list[LineSegment]
    ) -> list[tuple[tuple[tuple, tuple], LineSegment]]:
        """
        Builds the internal representation of this Path from a list of segments.

        Params
        ------
        xy_points: list[tuple[float, float]]
            List of x,y pairs.
        segments: (list[tuple[int, int]] | list[LineSegment])
            A list of either LineSegments or xy point pairs.
        
        Returns
        -------
        xypntindicies_2_segments: list[tuple[tuple[tuple, tuple], LineSegment]]
            A list of [xy_point, segment] pairs, one for each point in the given segments.
        """
        xypntindicies_2_segments: list[tuple[tuple[tuple, tuple], LineSegment]] = []

        # Insert new segments and register their existence.
        for s in segments:
            if not isinstance(s, LineSegment):
                s = cls._build_segment(xy_points, s)
            k = (s.xy1, s.xy2)
            xypntindicies_2_segments.append((k, s))
        
        return xypntindicies_2_segments

    @staticmethod
    def _get_xy_points_in_segment_order(segments: list[LineSegment]) -> list[tuple[float, float]]:
        """
        Returns
        -------
        ordered_xy_points: list[tuple[float, float]]
            The x,y point pairs in the given segments, reordered to be in the order found in the segments.
        """
        ordered_xy_points = []

        if len(segments) == 1:
            # special case: only one segment
            ordered_xy_points = [segments[0].xy1, segments[0].xy2]
        
        else:
            prev_segment: LineSegment = None
            for segment_idx, segment in enumerate(segments):
                if prev_segment is not None:
                    if prev_segment.xy1 in [segment.xy1, segment.xy2]:
                        xy_common, xy_other = prev_segment.xy1, prev_segment.xy2
                    else:
                        assert prev_segment.xy2 in [segment.xy1, segment.xy2]
                        xy_common, xy_other = prev_segment.xy2, prev_segment.xy1

                    if segment_idx == 1:
                        # special case, first segment
                        ordered_xy_points.append(xy_other)

                    ordered_xy_points.append(xy_common)

                prev_segment = segment
            
            # special case, last segment
            if segment.xy1 == ordered_xy_points[-1]:
                ordered_xy_points.append(segment.xy2)
            else:
                assert segment.xy2 == ordered_xy_points[-1]
                ordered_xy_points.append(segment.xy1)
        
        return ordered_xy_points