import copy
from typing import TYPE_CHECKING

from FileIO.Line import Line as FLine
from Geometry.LineSegment import LineSegment

if TYPE_CHECKING:
    from Trace.SingleTrace import _TraceLine


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
        # import here to avoid an import cycle
        from Trace.SingleTrace import _TraceLine

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
            idx1, idx2 = tuple(sorted([idx1, idx2]))
            new_segments_idxs.append(_TraceLine(s.source_line, idx1, idx2))
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
    
    def segment_xypnts_indicies(self, segment: LineSegment) -> tuple[int, int]:
        for s in self._xypntindicies_2_segments:
            if s[1] == segment:
                return s[0]
        assert False, f"Failed to locate segment {segment} in this path!"

    def segment_xypnts(self, segment: LineSegment) -> tuple[tuple[float, float], tuple[float, float]]:
        """
        Returns the XY points of a given LineSegment.

        Args:
            segment (LineSegment): The LineSegment to find the endpoints for.

        Returns:
            tuple[tuple[float, float], tuple[float, float]]: The two XY coordinates that make up the segment.
        """
        xy1_idx, xy2_idx = self.segment_xypnts_indicies(segment)
        xy1, xy2 = self.xy_points[xy1_idx], self.xy_points[xy2_idx]
        assert xy1 == segment.xy1
        assert xy2 == segment.xy2
        return xy1, xy2
    
    def insert_xypnt(self, new_xy_pnt: tuple[float, float], old_segment: LineSegment) -> tuple[LineSegment, LineSegment]:
        """ Inserts a new xy point between the two points given by old_segment.

        The old segment will be removed and two new segments made on either side of it.

        Parameters
        ----------
        new_xy_pnt : tuple[float, float]
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
        for old_xypnts2seg_idx, (k, s) in enumerate(self._xypntindicies_2_segments):
            if s == old_segment:
                old_xypnts2seg = (k, s)
        assert old_xypnts2seg
        self._xypntindicies_2_segments.remove(old_xypnts2seg)

        # insert the new point
        prev_xy_pnt = old_segment.xy1
        prev_xy_pnt_idx = self._xy_points.index(prev_xy_pnt)
        self._xy_points.insert(prev_xy_pnt_idx+1, new_xy_pnt)

        # create new segments
        prev_segment = LineSegment(old_segment.xy1, new_xy_pnt, old_segment.source_line)
        next_segment = LineSegment(new_xy_pnt, old_segment.xy2, old_segment.source_line)
        xy1_idx, xy12_idx, xy2_idx = self.xy_points.index(prev_segment.xy1), self.xy_points.index(prev_segment.xy2), self.xy_points.index(next_segment.xy2)
        assert xy1_idx >= 0
        assert xy12_idx >= 0
        assert xy2_idx >= 0

        # insert the new segments
        self._xypntindicies_2_segments.insert(old_xypnts2seg_idx+1, ((xy12_idx, xy2_idx), next_segment))
        self._xypntindicies_2_segments.insert(old_xypnts2seg_idx+1, ((xy1_idx, xy12_idx), prev_segment))

        return prev_segment, next_segment

    def remove_xypnt(self, xy_pnt: tuple[float, float]) -> tuple[list[LineSegment], list[LineSegment]]:
        """ Removes the given point from this instance.

        To do this, we remove the point and generate new segments,
        combining the properties of the adjacent segments.
        TODO combine the properties instead of just using the previous
        segment's properties.

        Parameters
        ----------
        xy_pnt : tuple[float, float]
            The point to be removed.

        Returns
        -------
        old_segments: list[LineSegment]
            The old line segments that were removed as part of removing the point.
            There will be two of these in the case of a non-forking trace.
        new_segments: list[LineSegment]
            The new line segments that were added as part of removing the point.
            There will be one of these in the case of a non-forking trace.
        """
        assert xy_pnt in self.xy_points
        
        # get the previous and next segments
        old_segments = self.segments_at_xypnt(xy_pnt)
        prev_segments: list[LineSegment] = []
        next_segments: list[LineSegment] = []
        for segment in old_segments:
            if segment.xy1 == xy_pnt:
                next_segments.append(segment)
            else:
                prev_segments.append(segment)
        
        # Get the indicies of the next segments.
        # We'll insert our new segments at these indicies.
        prev_xy2seg_indicies: list[int] = []
        next_xy2seg_indicies: list[int] = []
        for segment in prev_segments:
            prev_xypnts_idxs = self.segment_xypnts_indicies(segment)
            prev_xy2seg_indicies.append(self._xypntindicies_2_segments.index((prev_xypnts_idxs, segment)))
        for segment in next_segments:
            next_xypnts_idxs = self.segment_xypnts_indicies(segment)
            next_xy2seg_indicies.append(self._xypntindicies_2_segments.index((next_xypnts_idxs, segment)))
        assert all([idx >= 0 for idx in prev_xy2seg_indicies])
        assert all([idx >= 0 for idx in next_xy2seg_indicies])

        new_xy2seg_indicies = copy.copy(next_xy2seg_indicies)
        for i, idx in enumerate(new_xy2seg_indicies):
            for prev_idx in prev_xy2seg_indicies:
                if prev_idx < idx:
                    new_xy2seg_indicies[i] -= 1
        
        # remove the old point and segments
        pnt_idx = self.xy_points.index(xy_pnt)
        assert pnt_idx >= 0
        _xypntindicies_2_segments_length_old = len(self._xypntindicies_2_segments)
        self._xy_points.remove(xy_pnt)
        for s in copy.copy(self._xypntindicies_2_segments):
            if s[1] in old_segments:
                self._xypntindicies_2_segments.remove(s)
        _xypntindicies_2_segments_length_new = len(self._xypntindicies_2_segments)
        assert _xypntindicies_2_segments_length_old - _xypntindicies_2_segments_length_new == len(old_segments)

        # build new segments
        new_segments: list[LineSegment] = []
        prev_segment = prev_segments[0]
        for new_idx, next_segment in zip(reversed(new_xy2seg_indicies), reversed(next_segments)):
            xy1, xy2 = prev_segment.xy1, next_segment.xy2
            new_segment = LineSegment(xy1, xy2, next_segment.source_line)
            new_segments.append(new_segment)
            xy1_idx, xy2_idx = self.xy_points.index(xy1), self._xy_points.index(xy2)
            assert xy1_idx >= 0
            assert xy2_idx >= 0
            self._xypntindicies_2_segments.insert(new_idx, ((xy1_idx, xy2_idx), new_segment))

        return old_segments, new_segments

    @staticmethod
    def _build_segment(xy_points: list[tuple[float, float]], segment: "_TraceLine") -> LineSegment:
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
        x1, y1 = xy_points[segment.xy1_idx]
        x2, y2 = xy_points[segment.xy2_idx]
        source_line = segment.fline
        return LineSegment((x1, y1), (x2, y2), source_line=source_line)

    @classmethod
    def _build_segments(cls,
                        xy_points: list[tuple[float, float]], 
                        segments: list[tuple[int, int] | tuple[int, int, FLine]] | list[LineSegment]
    ) -> list[tuple[tuple[int, int], LineSegment]]:
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
        xypntindicies_2_segments: list[tuple[tuple[int, int], LineSegment]] = []

        # Insert new segments and register their existence.
        for s in segments:
            if not isinstance(s, LineSegment):
                s = cls._build_segment(xy_points, s)
            k = (xy_points.index(s.xy1), xy_points.index(s.xy2))
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