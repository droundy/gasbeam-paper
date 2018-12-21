use clapme::ClapMe;

/// Inputs
#[derive(Debug, ClapMe)]
pub struct Inputs {
    /// The desired flow rate
    flow_rate: f64,
}

fn main() {
    let input = Inputs::from_args();
}
