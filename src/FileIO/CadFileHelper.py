import re


class CadFileHelper:
    """ Assistant class to help parse gencad files, such as those exported by KiCad. """

    def __init__(self, start_tag: str | re.Pattern, end_tag: str | re.Pattern):
        self.start_tag = start_tag
        self.end_tag = end_tag

    def _matches(self, line: str, tag: str | re.Pattern):
        sline = line.strip()

        if isinstance(tag, str):
            return sline == tag
        
        else:
            match = tag.match(sline)
            return match is not None

    def start_matches(self, line: str):
        return self._matches(line, self.start_tag)

    def end_matches(self, line: str):
        return self._matches(line, self.end_tag)

    def get_next_region(self, lines: list[str]) -> tuple[list[str], list[str], list[str]]:
        """ Find the next region in the given lines as defined by the start and end tags.

        Parameters
        ----------
        lines : list[str]
            The lines in the file to search through.

        Returns
        -------
        pre_lines: list[str]
            The lines that come before the region.
        region_lines: list[str]
            The lines for the first found region, including the start and end tag lines.
        post_lines: list[str]
            The lines that come after the region.
        """
        in_region = False
        region_start = -1
        region_end = -1

        for line_idx, line in enumerate(lines):
            start_matches = False

            if self.start_matches(line):
                if not in_region:
                    start_matches = True
                    in_region = True
                    region_start = line_idx
                else:
                    if not self.end_matches(line):
                        raise RuntimeError
            
            if not start_matches and self.end_matches(line):
                if in_region:
                    in_region = False
                    region_end = line_idx
                    break
                else:
                    raise RuntimeError
                
        if region_start >= 0:
            post_lines = [] if region_end == len(lines)-1 else lines[region_end+1:]
            return lines[:region_start], lines[region_start:region_end+1], post_lines
        else:
            return lines, [], []