# Fluid Solver

This project implements a lumped element thermal/fluid simulation with the intent of being used in real time simulation. To ensure real-time execution several strategies are used:
- Dynamic memory allocation is using during construction, but avoided during model execution
- Thermo/fluid equations are in an explicit form
- Equations are calculated in a discrete time manner and iteration is bounded

## Solver phases
Simulation is split into 3 main phases.
- Construction
  - Thermo/fluid elements are added to the network, connected, and given parameters and initial conditions
- Simplification/Validation
  - Depending on settings some elements are lumped together to reduce the network complexity
  - Elements and connections are validated to ensure correctness
  - Any final run-tim buffers are created
- Execution
    - There are three modes of execution.
      - Free Run - Where the model will run forever until program termination
      - Fixed Length - The Model will run until a specified simulation time elapses
      - Step Mode - The user calls a step function and the model executes a single time step
    - The solution of a network is divided into 2 phases
     - First, all the edges execute and calculate flows given external inputs, internal state, and the connected nodes
     - Then all the nodes execute. Nodes calculate a new internal state based on external inputs, previous internal state, and the flows of connected edges

## Model structure
A `FluidNetwork` is composed of several fluid elements. These elements form an undirected graph.
Elements fall into these major categories:
- `ThermalNode` - A thermal mass with termperature state.
- `ThermalEdge` - An edge either end connected to either a `ThermalNode`, `FluidNode`, or a `ThermalBoundary`. `ThermalEdge`'s can model conduction, convection, or radiation.
- `ThermalBoundary` - A temperature boudary condition driven by an external input.
- `FluidNode` - A fixed control volume of fluid. With an internal state for both gas and liquid.
- `FluidEdge` - An edge that can model various components. Such as valves, pipes, orifices, external flow sources, lumped edges etc. It can be connected to a `FluidNode`, `FluidBoundary`, or `FluidJunction`
- `FluidBoundary` - A fluid boudary condition driven by an external input. Temperature, Pressure, liquid fraction, etc.
- `FluidJunction` - A placeholder used for building more complex networks. During the simplification phase these will be replaced with either
  - A `FluidNode` given a volume based on nearby edges.
  - A  lumped edge formed by connected edges. For example two valves in series may be lumped together into one edge based on the state of both valves.

## Usage
TODO