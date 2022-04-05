
mod simtools;
use simtools::{simsolver};
use simsolver::{*};
use simmodel::{*};

fn main() {
    
    let mut model = SpaceStateModel::new(2, 1, 2);

    model.set_mat_a(&[0.0, 1.0, -1.0, -1.0]).unwrap();
    model.set_mat_b(&[0.0, 1.0]).unwrap();
    model.init_state(&[0.0, 0.0]).unwrap();

    let mut solver = Simulator::<SpaceStateModel>::new(20.0, 0.1, SolverType::Euler, model);
    solver.run_sim();

    
    
}
