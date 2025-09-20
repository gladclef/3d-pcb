import dataclasses


@dataclasses.dataclass
class Line:
    v: str
    lineno: int
    sourcefile: str

    @staticmethod
    def from_file(file: str):
        ret: list[Line] = []

        with open(file, 'r') as fin:
            lineno = 0
            for l in fin:
                ret.append(Line(l, lineno, file))
                lineno += 1
        
        return ret