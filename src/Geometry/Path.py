from typing import TYPE_CHECKING

from Geometry.LineSegment import LineSegment


class Path:
    """
    Represents a 2D path in the XY plane.

    A path consists of back-to-back line segments.
    """
    def __init__(self, xy_points: list[tuple[float, float]], segments: list[tuple[int, int]] | list["LineSegment"]):
        # normalize input
        xy_points = [tuple(ab) for ab in xy_points]

        self._xy_points = xy_points
        self._xypntindicies_segments: list[tuple[tuple[tuple, tuple], LineSegment]] = []

        self._build_segments(segments)
    
    @property
    def xy_points(self) -> list[tuple[float, float]]:
        return self._xy_points

    @property
    def segments(self):
        return [s[1] for s in self._xypntindicies_segments]
    
    def segments_at_xypnt(self, pnt_or_pntidx: tuple[float, float] | int) -> list[LineSegment]:
        try:
            x, y = pnt_or_pntidx
            pnt_idx = self.xy_points.index(pnt_or_pntidx)
        except:
            pnt_idx: int = pnt_or_pntidx
        
        matching_segments = filter(lambda s: pnt_idx in s[0], self._xypntindicies_segments)
        return [s[1] for s in matching_segments]
    
    def segment_xypnts(self, segment: LineSegment) -> tuple[tuple[float, float], tuple[float, float]]:
        for s in self._xypntindicies_segments:
            if s[1] == segment:
                return s[0]
    
    def _build_segment(self, segment: tuple[int, int]) -> "LineSegment":
        x1, y1 = self._xy_points[segment[0]]
        x2, y2 = self._xy_points[segment[1]]
        return LineSegment((x1, y1), (x2, y2))

    def _build_segments(self, segments: list[tuple[int, int]] | list["LineSegment"]):
        """ Inserts the given segments into self._segments. Any segments not in this list are removed. """
        # normalize input
        segments: list[LineSegment] = [(s if isinstance(s, LineSegment) else self._build_segment(s)) for s in segments]

        # Get the point indicies that exist currently.
        # Any that aren't found during insertion will be removed.
        to_keep: list[tuple] = []
        to_remove: int = 0

        # Insert new segments and register their existance.
        xypntindicies_to_index = {s[0]: i for i, s in enumerate(self._xypntindicies_segments)}
        for s in segments:
            k = (s.xy1, s.xy2)
            if k not in xypntindicies_to_index:
                self._xypntindicies_segments.append((k, s))
                to_keep.append(self._xypntindicies_segments[-1])
            else:
                to_remove += 1

        # Remove keys that are not in the given list of segments.
        if to_remove > 0:
            new_xypntindicies_segments = []
            for xys in to_keep:
                new_xypntindicies_segments.append(xys)
            self._xypntindicies_segments = new_xypntindicies_segments