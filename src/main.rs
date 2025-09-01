use fluid_solver::Edge;
use fluid_solver::FluidNetwork;
use fluid_solver::Node;
use fluid_solver::elements::EdgeType;
use fluid_solver::elements::{ThermalBoundary, ThermalEdgeImpl, ThermalNode};

fn main() {
    // Create a simple network: ThermalBoundary -> ThermalNode -> ThermalEdge
    let mut network = FluidNetwork::new();

    // Add boundary
    let boundary = network.add_node(
        Node::ThermalBoundary(ThermalBoundary {
            name: None,
            temperature: 300.0,
        }),
        None,
    );

    // Add thermal node
    let thermal_node = network.add_node(
        Node::Thermal(ThermalNode {
            name: None,
            thermal_mass: 10.0,
            initial_temp: 280.0,
        }),
        Some("ThermNode"),
    );

    // Add edge between them
    let edge_impl = ThermalEdgeImpl {
        name: None,
        conductance: 5.0,
        emissivity: None,
    };
    network.add_edge(
        Edge {
            name: None,
            behavior: Box::new(edge_impl),
            edge_type: EdgeType::Thermal,
        },
        boundary,
        thermal_node,
        None,
    );

    println!(
        "Network built with {} nodes and {} edges",
        network.nodes.len(),
        network.edges.len()
    );
    println!("{:?}", network);

    // Test lookup by name
    if let Some(node_idx) = network.find_node_by_name("ThermNode") {
        println!("Found node 'ThermNode' at index {}", node_idx.index());
    } else {
        println!("Node 'ThermNode' not found");
    }

    if let Some(edge_idx) = network.find_edge_by_name("Edge_2") {
        println!("Found edge 'Edge_2' at index {}", edge_idx.index());
    } else {
        println!("Edge lookup failed");
    }
}
