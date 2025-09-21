# 3d-pcb
Converts KiCad gencad ".cad" files into 3d printable circuit boards. I want to be able to design a PCB, print it the same day, and iterate.

## Examples

Small example [hello_light](/examples/hello_light/):
<img src="/examples/hello_light/images/6_fabbed.JPG" width="600" />

## Usage

#### Installation

1. [Download and install python](https://www.python.org/downloads/)
2. [Download and install git](https://git-scm.com/downloads)
3. Create a new directory to download this project into. For example `%USERPROFILE%\AppData\Local\3D-PCB` on Windows or `~/prg/3D-PCB` on Linux/Mac.
4. Open a terminal in the newly created directory. On Windows start git bash. On Linux/Mac use the shortcut `Win+Space`/`Cmd+Space` and type "terminal". Then use the "cd" command to move to the newly created directory, for example "cd ~/AppData/Local/3D-PCB".
5. Download the project by running "git clone https://github.com/gladclef/3d-pcb.git" in the terminal.

#### 3D-PCB Creation

The steps to produce a 3D-PCB are:

1. Design schematic and layout
2. Run "python src/Board.py" from the downloaded git directory.
3. Create a new Blender project
    a. Import the traces, vias, and components
    b. Size a board to match
    c. Boolean-subtract the traces, vias, and components
    d. Export the ready-to-print board
4. Print the 3D-PCB
5. Fill the traces with 24 AWG tinned copper wire
6. Solder in the components

## Limitations

Below is a list of known limitations to layouts based on th current state of the project:

- Limited to 2 layers
- KiCad "GenCad" export files supported
- No forking traces (traces must go from one through-hole to another, not from a trace to a through-hole)
- No sharp bends (90 degree bends are OK)
- Trace connections to through-holes should provide ~0.5mm of space North of the trace.
- Vias [aren't currently supported](https://github.com/gladclef/3d-pcb/issues/1)
- Through-hole components only
