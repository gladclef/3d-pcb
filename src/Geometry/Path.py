import copy

from FileIO.Line import Line as FLine
from Geometry.Line import Line
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
    
    def _check_structure_is_valid(self):
        assert all([self.segments[i].xy1 == self.segments[i+1].xy0 for i in range(len(self.segments)-1)]), "Segments not in order!"

    @property
    def xy_points(self) -> list[Point]:
        """
        Returns the list of XY points that define this path.
        """
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
    
    def get_previous_segment(self, segment: LineSegment) -> LineSegment | None:
        ret = None if segment == self.segments[0] else self.segments[self.segments.index(segment)-1]
        assert (ret is None) or (ret.xy1 == segment.xy0)
        return ret
    
    def get_next_segment(self, segment: LineSegment) -> LineSegment | None:
        ret = None if segment == self.segments[-1] else self.segments[self.segments.index(segment)+1]
        assert (ret is None) or (ret.xy0 == segment.xy1)
        return ret

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

        self._check_structure_is_valid()

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

        self._check_structure_is_valid()

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

        self._check_structure_is_valid()

        return old_segments, new_segments
    
    def remove_segment(
            self,
            segment: LineSegment,
            modify_adjoining_segments = True,
            max_distance_from_center = 2.0
        ) -> tuple[LineSegment|None, LineSegment|None]:
        """ Removes the given segment from this path.

        If modify_adjoining_segments is False, then the removal is done with a
        simple deletion of the segment.

        If modify_adjoining_segments is True, then the segments before and
        after the given segment are modified to extend to the intersection
        points of their respective lines.

        Example of before removal, removal with modify_adjoining_segments == False,
        and removal with modify_adjoining_segments == True:

                Before        No Modification     Modification
            A------------B    A                  A-------------AD
                         |      \                              /
                         C        \                           /
                        /           \                        / 
                       /              \                     /  
                      D                 D                  D   

        Parameters
        ----------
        segment : LineSegment
            The line segment to be removed.
        modify_adjoining_segments : bool, optional
            True to adjust the previous and next segments, by default True
        max_distance_from_center : float, optional
            The maximum distance from the center of segment that the
            intersection point of the previous and next segments can be, in the
            case that modify_adjoining_segments is True.
        
        Returns
        -------
        tuple[LineSegment|None, LineSegment|None]:
            The previous and next line segments to the given segment.
        """
        assert segment in self.segments

        # get the previous and next segments
        prev_segment = self.get_previous_segment(segment)
        next_segment = self.get_next_segment(segment)

        # simple case, modify_adjoining_segments == False
        if modify_adjoining_segments is False:
            for pnt in segment.xy_points:
                self.remove_xypnt(pnt)
            self._check_structure_is_valid()
            return prev_segment, next_segment
        
        # simple case, segment is at the beggining or end
        if segment == self.segments[0]:
            assert prev_segment is None
            self.remove_xypnt(segment.xy1)
            self._check_structure_is_valid()
            return prev_segment, next_segment
        elif segment == self.segments[-1]:
            assert next_segment is None
            self.remove_xypnt(segment.xy0)
            self._check_structure_is_valid()
            return prev_segment, next_segment

        # start by finding the previous and next segments
        segment_idx = self.segments.index(segment)
        prev_segment = self.segments[segment_idx-1]
        next_segment = self.segments[segment_idx+1]
        assert prev_segment.xy1 == segment.xy0
        assert next_segment.xy0 == segment.xy1

        # check for lines that are parallel
        center_pnt = segment.distance_along_line(segment.length, segment.xy0)
        if prev_segment.is_parallel_to(next_segment):
            self.segments.remove(segment)
            new_pnt = prev_segment.distance_along_line(max_distance_from_center, center_pnt)
            prev_segment.xy1 = new_pnt
            next_segment.xy0 = new_pnt
            self._check_structure_is_valid()
            return prev_segment, next_segment

        # find the segments intersection point
        prev_line = Line.from_two_points(prev_segment.xy0, prev_segment.xy1)
        next_line = Line.from_two_points(next_segment.xy0, next_segment.xy1)
        intersection_pnt = prev_line.intersection(next_line)
        if intersection_pnt.distance(center_pnt) <= max_distance_from_center:
            self.segments.remove(segment)
            prev_segment.xy1 = intersection_pnt
            next_segment.xy0 = intersection_pnt
            self._check_structure_is_valid()
            return prev_segment, next_segment
        else:
            center_to_intersection = Line.from_two_points(center_pnt, intersection_pnt)
            closest_pnt = center_to_intersection.distance_along_line(max_distance_from_center, center_pnt)
            self.segments.remove(segment)
            prev_segment.xy1 = closest_pnt
            next_segment.xy0 = closest_pnt
            self._check_structure_is_valid()
            return prev_segment, next_segment

    def change_segment_points(
            self,
            segment: LineSegment,
            new_xy0: Point = None,
            new_xy1: Point = None
        ) -> tuple[LineSegment|None, LineSegment|None]:
        """ Moves the points for the given segment and related components to
        the new locations.

        Parameters
        ----------
        segment : LineSegment
            The segment to be moved.
        new_xy0 : Point, optional
            The new point to move the start of the segment to, or None. By default None
        new_xy1 : Point, optional
            The new point to move the end of the segment to, or None. By default None

        Returns
        -------
        tuple[LineSegment|None, LineSegment|None]
            The previous and next segments that are next to the given segment.
        """
        assert segment in self.segments

        # simple case, no changes
        if new_xy0 is None and new_xy1 is None:
            return
        
        # get the previous and next segments
        prev_segment = None if segment == self.segments[0] else self.segments[self.segments.index(segment)-1]
        next_segment = None if segment == self.segments[-1] else self.segments[self.segments.index(segment)+1]
        
        # move the segment points
        if new_xy0 is not None:
            segment.xy0 = new_xy0
            if prev_segment is not None:
                prev_segment.xy1 = new_xy0
        else:
            segment.xy1 = new_xy1
            if next_segment is not None:
                next_segment.xy0 = new_xy1
                
        self._check_structure_is_valid()

        return prev_segment, next_segment