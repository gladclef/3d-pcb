import numpy as np
import vtk

from Path import Path
import geometry as geo
from Segment import Segment

class PathCorner:
    def __init__(self, parent: Path, segments: tuple[Segment, Segment], bend_radius: float):
        """
        Parameters
        ----------
        parent : Path
            The parent instance that contains this instance.
        segments : tuple[Segment, Segment]
            The two segments that this corner is found on. They should
            share one of the end points.
        bend_radius : float
            How large of a radius the bend should be.
        """
        # verify that the segments share an end point, and that the shared point is in the middle
        if segments[0].xy1 == segments[1].xy0:
            if segments[1].xy1 == segments[0].xy0:
                raise ValueError("The segments must share only one common end point.")
        elif segments[1].xy1 == segments[0].xy0:
            segments = [segments[1], segments[0]]
        else:
            raise ValueError("The segments must share a common end point.")
        
        # check for a valid bend radius
        if bend_radius < parent.shape.radius:
            raise RuntimeError("Error in PathCorner.min_corner_increment_angle: " + f"bend_radius {bend_radius} is less than the pipe shape radius {parent.shape.radius}!")

        self.parent = parent
        self.segments = segments
        self._bend_radius = bend_radius
        self.debug = False

    @property
    def bend_radius(self):
        return self._bend_radius
    
    @bend_radius.setter
    def bend_radius(self, val: float):
        self._bend_radius = val

    def get_arc_seg_intersections(self):
        """
        Find the location where the arc intersects both line segments.

        To do this we:
            1. assume both lines are infinitely long
            2. get the points from the lines that is a distance (the bend radius) away from the line
            3. choose the point for each line that is in the same direction as the other line
            4. adjust both lines towards each other to go through the closer points
            5. find the intersection of the adjusted lines
        """
        debug = True
        a0, a1 = self.segments[0].xy0, self.segments[0].xy1
        b0, b1 = self.segments[1].xy0, self.segments[1].xy1

        # check that the middle point is shared
        if a1 != b0:
            raise ValueError("Error in PathCorner.get_arc_seg_intersections(): " + f"the midpoint between the segments ({a0},{a1}) and ({b0},{b1}) isn't shared! " + f"{a1} != {b0}")

        # don't search for the arc segment intersection of two segments with the same(ish) slope
        angle_a = geo.line_angle((a0, a1))
        angle_b = geo.line_angle((b0, b1))
        angle_diff = geo.normalize_angle(angle_b - angle_a)
        if self.debug:
            print(f"{angle_diff=}")
        if angle_diff < 1/1e6 or angle_diff > (2*np.pi - 1/1e6):
            tangent = geo.get_tangent_line((a0, a1), a1)
            return geo.distance_along_line(tangent, self.bend_radius, a1)

        # 2. get the bend radius distance points away from each line
        def distance_along_tangent(line: tuple[tuple[float, float], tuple[float, float]], from_point: tuple[float, float]) -> tuple[tuple[float, float], tuple[float, float]]:
            tangent = geo.get_tangent_line(line, from_point)
            if self.debug:
                line_slope_intercept = geo.line_two_points_to_slope_intercept(line)
                tangent_slope_intercept = geo.line_two_points_to_slope_intercept(tangent)
                print(f"{line_slope_intercept=}\n{tangent_slope_intercept=}")
            ret0 = geo.distance_along_line((tangent[0], tangent[1]), self.bend_radius, from_point)
            ret1 = geo.distance_along_line((tangent[1], tangent[0]), self.bend_radius, from_point)
            return ret0, ret1
        pnts_dist_adj_a = distance_along_tangent((a0, a1), a1)
        if self.debug:
            print(f"{pnts_dist_adj_a=}")
        pnts_dist_adj_b = distance_along_tangent((b0, b1), a1)
        if self.debug:
            print(f"{pnts_dist_adj_b=}")

        # 3. choose the point for each line that is in the same direction as the other line
        def get_point_on_same_side(line: tuple[tuple[float, float], tuple[float, float]],
                                   ref_point: tuple[float, float],
                                   two_points: tuple[tuple[float, float], tuple[float, float]]) -> tuple[float, float]:
            if geo.is_point_on_right(line, two_points[0]):
                if geo.is_point_on_right(line, ref_point):
                    return two_points[0]
                else:
                    return two_points[1]
            else: # point 0 is on the left
                if not geo.is_point_on_right(line, ref_point):
                    return two_points[0]
                else:
                    return two_points[1]
        pnt_dist_adj_a = get_point_on_same_side((a0, a1), b1, pnts_dist_adj_a)
        pnt_dist_adj_b = get_point_on_same_side((b1, b0), a0, pnts_dist_adj_b)
        if self.debug:
            print(f"{pnt_dist_adj_a=}\n{pnt_dist_adj_b=}")

        # 4. adjust both lines towards each other to go through the closer points
        def get_line_through_point(parallel_line: tuple[tuple[float, float], tuple[float, float]], point: tuple[float, float]) -> tuple[tuple[float, float], tuple[float, float]]:
            (x1, y1), (x2, y2) = parallel_line
            slope, _ = geo.line_two_points_to_slope_intercept(parallel_line)

            # check for infinity or zero
            if abs(slope) >= geo.INF_THRESH:
                if y2 > y1:
                    return point, (point[0], point[1] + 10)
                else:
                    return point, (point[0], point[1] - 10)
            elif abs(slope) <= geo.ZERO_THRESH:
                if x2 > x1:
                    return point, (point[0] + 10, point[1])
                else:
                    return point, (point[0] - 10, point[1])
            
            # get the new line
            return point, (point[0]+10, slope*10 + point[1])
        adj_line_a = get_line_through_point((a0, a1), pnt_dist_adj_a)
        adj_line_b = get_line_through_point((b0, b1), pnt_dist_adj_b)
        if self.debug:
            print(f"{adj_line_a=}\n{adj_line_b=}")

        # 5. find the intersection of the adjusted lines
        adj_lines_intersection = geo.lines_intersection(adj_line_a, adj_line_b)
        if self.debug:
            print(f"{adj_lines_intersection=}")

        return adj_lines_intersection

    def get_arc_properties(self) -> tuple[float, float, float, tuple[float, float]]:
        """ Get common properties about the arc that this path corner follows.

        Returns
        -------
        angle_diff: float
            The difference between the angle of the first and second segments.
            This value plus the angle of the first segment will equal the angle
            of second segment (angle_seg1 + angle_diff == angle_seg2).
        mid_angle: float
            The angle that is halfway between the angle of the first and second segments.
        arc_length: float
            The total length of the arc. Note that this is the length for a perfect
            arc. The actual length will be different due to the path corner being
            composed of multiple increments.
        arc_center: tuple[float, float]
            The x,y point at the center of the circle that the arc is a part of.
        """

        # Get some segment properties
        angle_a = geo.line_angle(self.segments[0].xy_points)
        angle_b = geo.line_angle(self.segments[1].xy_points)
        if self.debug:
            print(f"{angle_a=}\n{angle_b=}")

        # Compute the arc angles.
        # This is more easily done by adjusting the angles such that angle_a = 0.
        angle_b = geo.normalize_angle(angle_b - angle_a)
        angle_diff = angle_b if angle_b < np.pi else angle_b-2*np.pi
        mid_angle = geo.normalize_angle(angle_diff/2 + angle_a)
        if self.debug:
            print(f"{angle_b=}\n{angle_diff=}\n{angle_diff/2+angle_a=}\n{mid_angle=}")

        # Compute the length of the arc necessary based on
        # the given bend radius and the two angles a and b.
        circle_proportion = abs(angle_diff) / (2*np.pi)
        arc_length = circle_proportion * 2*np.pi * self.bend_radius

        # Get the center for where to place the arc
        arc_center = self.get_arc_seg_intersections()

        return angle_diff, mid_angle, arc_length, arc_center

    def get_center_line_points(self):
        """
        Get the points that define the path for the path corner.

        These are the points along the "center line" of the path corner,
        as pairs of (xy, θ), where:
            - xy is the x,y location of the point
            - θ is the direction that point is facing (0-2π)
        
        The first point in the path corner corresponds to the end of the
        first segment, and the last point corresponds to the start of the
        second segment. There will therefore always be at least two points.

        Returns
        -------
        center_line_points: list[tuple[tuple[float, float], float]]
            Pairs of (xy, θ), one per point on the center line.
        """
        corner_increment_angle = np.deg2rad(5)
        angle_diff, mid_angle, arc_length, arc_center = self.get_arc_properties()

        # Determine the number of increments to be used.
        # Note that the start of the first increment and the end of the last
        # increment is where the segment ends will be.
        nincrements = int(abs(angle_diff) / (corner_increment_angle*2 - np.deg2rad(1))) - 1
        nincrements = max(nincrements, 1)

        # Find the tangent from the arc center to the path corner (aka the end of the first segment).
        seg_tangent = geo.get_tangent_line(self.segments[0].xy_points, self.segments[0].xy1)
        tan_tst_pnt = geo.distance_along_line(seg_tangent, 1.0, self.segments[0].xy1)
        if geo.is_point_on_right(self.segments[0].xy_points, arc_center):
            if geo.is_point_on_right(self.segments[0].xy_points, tan_tst_pnt):
                seg_tangent = (seg_tangent[1], seg_tangent[0])
                tan_tst_pnt2 = geo.distance_along_line(seg_tangent, 1.0, self.segments[0].xy1)
                assert not geo.is_point_on_right(self.segments[0].xy_points, tan_tst_pnt2)
        else:
            if not geo.is_point_on_right(self.segments[0].xy_points, tan_tst_pnt):
                seg_tangent = (seg_tangent[1], seg_tangent[0])
                tan_tst_pnt2 = geo.distance_along_line(seg_tangent, 1.0, self.segments[0].xy1)
                assert geo.is_point_on_right(self.segments[0].xy_points, tan_tst_pnt2)

        # Get the angle for the start of the first increment.
        seg_tan_angle = geo.line_angle(seg_tangent)
        first_tan_angle = seg_tan_angle

        # Build out the list of increments.
        first_angle = self.segments[0].angle
        inc_angle = angle_diff / nincrements
        ret: list[tuple[tuple[float, float], float]] = []
        for i in range(nincrements+1):
            tan_angle = first_tan_angle + i*inc_angle
            pnt_xy = np.cos(tan_angle)*self.bend_radius, np.sin(tan_angle)*self.bend_radius
            pnt_xy = pnt_xy[0] + arc_center[0], pnt_xy[1] + arc_center[1]
            pnt_angle = first_angle + i*inc_angle
            ret.append((pnt_xy, pnt_angle))
        
        return ret

    def draw_corner(self):
        """ For debugging: draw the corner center line to be followed during shape generation. """
        import matplotlib.pyplot as plt
        
        angle_diff, mid_angle, arc_length, arc_center = self.get_arc_properties()
        center_line_points = self.get_center_line_points()

        # create the plot
        fig, ax = plt.subplots(figsize=(10,10))

        # draw the segments
        dx0, dy0 = self.segments[0].x1 - self.segments[0].x0, self.segments[0].y1 - self.segments[0].y0
        ax.arrow(*self.segments[0].xy0, dx0, dy0, color="tab:blue")
        dx1, dy1 = self.segments[1].x1 - self.segments[1].x0, self.segments[1].y1 - self.segments[1].y0
        ax.arrow(*self.segments[1].xy0, dx1, dy1, color="tab:blue")

        # draw the center point
        ax.add_patch(plt.Circle(arc_center, self.bend_radius/10, color="tab:purple"))

        # draw the arc
        for pnt_xy, pnt_ang in center_line_points:
            ax.arrow(*pnt_xy, np.cos(pnt_ang)*self.bend_radius/2, np.sin(pnt_ang)*self.bend_radius/2, color="black")
        for i in range(len(center_line_points)-1):
            pnt0_xy, pnt0_ang = center_line_points[i]
            pnt1_xy, pnt1_ang = center_line_points[i+1]
            ax.plot([pnt0_xy[0], pnt1_xy[0]], [pnt0_xy[1], pnt1_xy[1]], color="tab:orange")
        
        # show the plot
        ax.set_aspect('equal')
        ax.set_xlim(-3.5, 3.5)
        ax.set_ylim(-3.5, 3.5)
        plt.show(block=False)


if __name__ == "__main__":
    import json

    import matplotlib.pyplot as plt

    from SingleTrace import SingleTrace
    from PipeShape import PipeBasicBox
    from units import *

    vidframes_out = None #open("output/renders/input.txt", "w")
    unittests_out = open("tests/test_PathCorners_regression_vals.json", "w")

    if unittests_out is not None:
        unittests_out.write("[\n")

    for i, deg100 in enumerate(range(14500, -14500, -25)):
        deg = deg100 / 100
        rad = np.deg2rad(deg)

        xy_points = [
            [-3.5, 0],
            [0, 0],
            [0+np.cos(rad)*3.5, 0+np.sin(rad)*3.5]
        ]
        segments = [
            [0, 1],
            [1, 2]
        ]
        shape = PipeBasicBox(awg2mm(26))
        trace = SingleTrace(xy_points, segments, shape)
        corner = PathCorner(trace, trace.segments, 1)

        if vidframes_out is not None:
            corner.draw_corner()
            plt.savefig(f"output/renders/{i}_{deg}.png")
            plt.close('all')
            vidframes_out.write(f"file '{i}_{deg}.png'\n")
            vidframes_out.write(f"duration {1/30}\n")
            # ffmpeg command: something like "ffmpeg -f concat -i 'input.txt' -c:v -r 20 out.mp4"
        
        if unittests_out is not None:
            angle_diff, mid_angle, arc_length, arc_center = corner.get_arc_properties()
            center_line_points = corner.get_center_line_points()
            regression_val = {
                "deg100": deg100,
                "deg": deg,
                "rad": rad,
                "angle_diff": angle_diff,
                "mid_angle": mid_angle,
                "arc_length": arc_length,
                "arc_center": arc_center,
                "center_line_points": center_line_points
            }
            regression_str = "" if i == 0 else ",\n"
            regression_str += json.dumps(regression_val).replace("\n", " ")
            unittests_out.write(regression_str)
    
    if vidframes_out is not None:
        vidframes_out.close()
    if unittests_out is not None:
        unittests_out.write("\n]")
        unittests_out.close()