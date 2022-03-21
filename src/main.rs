
mod SimTools;
use SimTools::{SimSolver};
use SimSolver::{*};
use SimModel::{*};


fn main() {
    let mut solver = Simulator::<SpaceStateModel>::new(10.0, 0.01, SolverType::Euler);
    let mut model = SpaceStateModel::new(2, 1, 2);

    solver.run_sim();
    model.init_state(&[0.0, 0.0]);

    println!("Hello, world!");
}
