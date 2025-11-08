import copy

from FileIO.Line import Line as FLine
from Geometry.LineSegment import LineSegment
from Geometry.Point import Point


class Path:
    """
    Represents a 2D path in the XY plane.

    A path consists of back-to-back line segments.
    """

    def __init__(self,
                 source_lines: list[FLine],
                 segments: list[LineSegment]
    ):
        # sanity check/normalize input
        assert isinstance(source_lines, list)
        assert all([isinstance(l, FLine) for l in source_lines])
        assert isinstance(segments, list)
        assert all([isinstance(s, LineSegment) for s in segments])

        self.source_lines = source_lines
        self.segments = segments
    
    @property
    def xy_points(self) -> list[Point]:
        """
        Returns the list of XY points that define this path.
        """
        assert all([self.segments[i].xy1 == self.segments[i+1].xy0 for i in range(len(self.segments)-1)])
        return [s.xy0 for s in self.segments] + [self.segments[-1].xy1]

    @property
    def edges(self) -> list[tuple[Point, Point]]:
        return [s.xy_points for s in self._segments]
    
    def segments_at_xypnt(self, pnt: Point) -> list[LineSegment]:

        """
        Returns a list of LineSegments that contain the given XY point.

        Args:
            pnt Point: An XY point found in self.xy_points.

        Returns:
            list[LineSegment]: The segments that include the specified point.
        """
        assert pnt in self.xy_points
        matching_segments = list(filter(lambda s: pnt in s.xy_points, self.segments))
        assert 1 <= len(matching_segments) <= 2
        if len(matching_segments) == 1 or matching_segments[0].xy1 == pnt:
            return matching_segments
        else:
            return [matching_segments[1], matching_segments[0]]
    
    def insert_xypnt(self, new_xy_pnt: Point, old_segment: LineSegment) -> tuple[LineSegment, LineSegment]:
        """ Inserts a new xy point between the two points given by old_segment.

        The old segment will be removed and two new segments made on either side of it.

        Parameters
        ----------
        new_xy_pnt : Point
            The new point to be inserted into this path.
        old_segment : LineSegment
            The old segment to be split into two segments.

        Returns
        -------
        tuple[LineSegment, LineSegment]
            The two new segments added around the new point.
        """
        assert old_segment in self.segments

        # remove the old segment
        old_index = self.segments.index(old_segment)
        self.segments.remove(old_segment)

        # create new segments
        prev_segment = LineSegment(old_segment.xy0, new_xy_pnt, old_segment.source_line)
        next_segment = LineSegment(new_xy_pnt, old_segment.xy1, old_segment.source_line)

        # insert the new segments
        self.segments.insert(old_index+1, prev_segment)
        self.segments.insert(old_index+2, next_segment)

        return prev_segment, next_segment
    
    def append_xypnt(self, new_xy_pnt: Point, at_start: bool, fline: FLine = None) -> LineSegment:
        """ Adds a new xy point at one end or the other of this path.

        Parameters
        ----------
        new_xy_pnt : Point
            The new point to be added to this path.
        at_start : bool
            True if the point should be added to the beggining of the path.
            False if it should be added to the end of the path.
        fline : FLine
            The file line from which the new segment will be created.

        Returns
        -------
        LineSegment
            The new segments added at this point.
        """
        # create the new segment
        xy0 = new_xy_pnt if at_start else self.segments[-1].xy1
        xy1 = self.segments[0].xy0 if at_start else new_xy_pnt
        new_segment = LineSegment(xy0, xy1, fline)

        # insert the new segments
        if at_start:
            self.segments.insert(0, new_segment)
        else:
            self.segments.append(new_segment)

        return new_segment

    def remove_xypnt(self, xy_pnt: Point) -> tuple[list[LineSegment], list[LineSegment]]:
        """ Removes the given point from this instance.

        To do this, we remove the point and generate new segments,
        combining the properties of the adjacent segments.
        TODO combine the properties instead of just using the previous
        segment's properties.

        Parameters
        ----------
        xy_pnt : Point
            The point to be removed.

        Returns
        -------
        old_segments: list[LineSegment]
            The old line segments that were removed as part of removing the point.
            There will be at most two of these.
        new_segments: list[LineSegment]
            The new line segments that were added as part of removing the point.
            There will be at most one of these.
        """
        assert xy_pnt in self.xy_points
        
        # get the previous and next segments
        old_segments = self.segments_at_xypnt(xy_pnt)
        if len(old_segments) == 1:
            if old_segments[0].xy0 == xy_pnt:
                prev_segment, next_segment = old_segments[0], None
            elif old_segments[0].xy1 == xy_pnt:
                prev_segment, next_segment = None, old_segments[0]
            else:
                raise RuntimeError
        else:
            prev_segment, next_segment = old_segments[0], old_segments[1]
        
        # remove the old segments
        for segment in old_segments:
            self.segments.remove(segment)

        # build new segments
        new_segments: list[LineSegment] = []
        if len(old_segments) > 1:
            new_segment = LineSegment(prev_segment.xy0, next_segment.xy1, prev_segment.fline)
            new_segments.append(new_segment)
        else:
            # nothing to do, no new segments to build
            pass

        return old_segments, new_segments

    def change_segment_points(self, segment: LineSegment, new_xy0: Point = None, new_xy1: Point = None) -> LineSegment:
        assert segment in self.segments
        new_pnt = new_xy0 or new_xy1

        # simple case, no changes
        if new_xy0 is None and new_xy1 is None:
            return segment
        
        # both points change
        if new_xy0 is not None and new_xy1 is not None:
            segment = self.change_segment_points(segment, new_xy0, None)
            return self.change_segment_points(segment, None, new_xy1)

        # special case, only one segment to this path
        if len(self.segments) == 1:
            prev_seg, next_seg = self.insert_xypnt(new_pnt)
            if new_xy0 is not None:
                self.remove_xypnt(prev_seg.xy0)
                return next_seg
            if new_xy1 is not None:
                self.remove_xypnt(next_seg.xy1)
                return prev_seg
        
        # special case, this segment is at one of the ends
        if (segment == self.segments[0] and new_xy0 is not None) or (segment == self.segments[-1] and new_xy1 is not None):
            at_start = segment == self.segments[0]
            fline = segment.source_line
            new_segment = self.append_xypnt(new_pnt, at_start, fline)
            self.remove_xypnt(segment.xy1 if at_start else segment.xy0)
            return new_segment
        
        # standard_case, this segment is in the middle
        prev_seg, next_seg = self.insert_xypnt(new_pnt, segment)
        if new_xy0 is not None:
            self.remove_xypnt(segment.xy0)
        else:
            self.remove_xypnt(segment.xy1)