from abc import ABC

class AbstractPath(ABC):
    """
    Represents a 2D path in the XY plane.

    A path consists of back-to-back line segments.
    """
    def __init__(self, xy_points: list[tuple[float, float]], segments: list[tuple[int, int]] | list["Segment"]):
        # import here to avoid cyclic imports
        from Segment import Segment

        # normalize input
        xy_points = [tuple(ab) for ab in xy_points]

        self._xy_points = xy_points
        self._xypntindicies_segments: list[tuple[tuple[tuple, tuple], Segment]] = []
        
        self._build_segments(segments)
    
    @property
    def xy_points(self) -> list[tuple[float, float]]:
        return self._xy_points

    @property
    def segments(self):
        return [s[1] for s in self._xypntindicies_segments]

    def _build_segments(self, segments: list[tuple[int, int]] | list["Segment"]):
        """ Inserts the given segments into self._segments. Any segments not in this list are removed. """
        # import here to avoid cyclic imports
        from Segment import Segment

        # normalize input
        segments: list[Segment] = [(s if isinstance(s, Segment) else Segment(self, s)) for s in segments]

        # Get the point indicies that exist currently.
        # Any that aren't found during insertion will be removed.
        to_keep: list[tuple] = []
        to_remove: int = 0

        # Insert new segments and register their existance.
        xypntindicies_to_index = {s[0]: i for i, s in enumerate(self._xypntindicies_segments)}
        for s in segments:
            k = s.xy_points
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