import re

from FileIO.Line import Line as FLine

class CadFileHelper:
    """
    Assistant class to help parse gencad files, such as those exported by KiCad.
    """

    def __init__(self, start_tag: str | re.Pattern, end_tag: str | re.Pattern, allow_multiple_starts=False, ignore_false_endings=False):
        """
        Initialize the `CadFileHelper` instance.

        Parameters
        ----------
        start_tag : str or re.Pattern
            The starting tag or pattern to identify the beginning of a region.
        end_tag : str or re.Pattern
            The ending tag or pattern to identify the end of a region.
        allow_multiple_starts : bool, optional
            True to ignore multiple starting match lines without an ending match line. Default is False.
        ignore_false_endings : bool, optional
            Whether to ignore lines matching `end_tag` when not in a region. Default is False.
        """
        self.start_tag = start_tag
        """ The starting tag or pattern for identifying regions. """
        self.end_tag = end_tag
        """ The ending tag or pattern for identifying regions. """
        self.allow_multiple_starts = allow_multiple_starts
        """ True to ignore multiple starting match lines without an ending match line. """
        self.ignore_false_endings = ignore_false_endings
        """ Whether to ignore lines matching `end_tag` when not in a region. """

    def _matches(self, line: FLine, tag: str | re.Pattern) -> bool:
        """
        Check if the given line matches the specified tag or pattern.

        Parameters
        ----------
        line : FLine
            The line to check.
        tag : str or re.Pattern
            The tag or pattern to match against.

        Returns
        -------
        bool
            True if the line matches the tag, False otherwise.

        """
        sline = line.v.strip()

        if isinstance(tag, str):
            return sline == tag
        
        else:
            match = tag.match(sline)
            return match is not None

    def start_matches(self, line: FLine) -> bool:
        """
        Check if the given line matches the start tag.

        Parameters
        ----------
        line : FLine
            The line to check.

        Returns
        -------
        bool
            True if the line matches the start tag, False otherwise.

        """
        return self._matches(line, self.start_tag)

    def end_matches(self, line: FLine) -> bool:
        """
        Check if the given line matches the end tag.

        Parameters
        ----------
        line : FLine
            The line to check.

        Returns
        -------
        bool
            True if the line matches the end tag, False otherwise.

        """
        return self._matches(line, self.end_tag)

    def get_next_region(self, lines: list[FLine]) -> tuple[list[FLine], list[FLine], list[FLine]]:
        """
        Find the next region in the given lines as defined by the start and end tags.

        Parameters
        ----------
        lines : list[FLine]
            The lines in the file to search through.

        Returns
        -------
        pre_lines: list[FLine]
            The lines that come before the region.
        region_lines: list[FLine]
            The lines for the first found region, including the start and end tag lines.
        post_lines: list[FLine]
            The lines that come after the region.

        Raises
        ------
        RuntimeError
            If a second start tag is encountered within a region or if an end tag is found outside of a region (and `ignore_false_endings` is False).

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
                        if not self.allow_multiple_starts:
                            raise RuntimeError("Error in CadFileHelper.get_next_region(): expected end to region but instead found second start!")
            
            if not start_matches and self.end_matches(line):
                if in_region:
                    in_region = False
                    region_end = line_idx
                    break
                else:
                    if not self.ignore_false_endings:
                        raise RuntimeError("Error in CadFileHelper.get_next_region(): found end region but we're not in a region!")
                
        if region_start >= 0:
            if region_end == -1:
                region_end = len(lines)-1
            post_lines = [] if region_end == len(lines)-1 else lines[region_end+1:]
            return lines[:region_start], lines[region_start:region_end+1], post_lines
        else:
            return lines, [], []