class Point:
    def __init__(self, x: float, y: float, delta: float=1e-4):
        self.x = x
        self.y = y
        self.delta = delta

    def almost_equal(self, other: "Point", delta: float=None) -> bool:
        diff = self - other
        xdiff, ydiff = abs(diff.x), abs(diff.y)
        if delta is None:
            delta = self.delta
        if xdiff < delta and ydiff < delta:
            return True
        return False
    
    def __add__(self, other: "Point") -> "Point":
        return Point(self.x + other.x, self.y + other.y, max(self.delta, other.delta))
    
    def __sub__(self, other: "Point") -> "Point":
        return Point(self.x - other.x, self.y - other.y, max(self.delta, other.delta))
    
    def __mul__(self, other: "Point") -> "Point":
        return Point(self.x * other.x, self.y * other.y, max(self.delta, other.delta))
    
    def __truediv__(self, other: "Point") -> "Point":
        return Point(self.x / other.x, self.y / other.y, max(self.delta, other.delta))
    
    def __abs__(self) -> "Point":
        return Point(abs(self.x), abs(self.y), self.delta)
    
    def __hash__(self):
        return hash(self.x + self.y + self.delta)
    
    def __eq__(self, other: "Point"):
        if not hasattr(other, "x") or not hasattr(other, "y") or not hasattr(other, "delta"):
            return False
        if self.x == other.x and self.y == other.y and self.delta == other.delta:
            return True
        return False
    
    def __repr__(self) -> str:
        nplaces = len(f"{self.delta:f}".split(".")[1].rstrip("0"))
        fstr = f":0.0{nplaces}f"
        return ("<{"+fstr+"},{"+fstr+"}>").format(self.x, self.y)
    
class Pointz(Point):
    def __init__(self, x: float, z: float):
        super().__init__(x, z)
    
    @property
    def z(self) -> float:
        return self.y
    
    @z.setter
    def z(self, value: float):
        self.y = value


if __name__ == "__main__":
    a = Point(0, 1)
    b = Point(2, 3)

    print(a - b)