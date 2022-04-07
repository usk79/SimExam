extern crate nalgebra as na;
use na::{U2, U3, Dynamic, ArrayStorage, VecStorage, Matrix, OMatrix, DMatrix};

mod simtools;
use simtools::{simsolver};
use simsolver::{*};
use simmodel::{*};

struct NewModel {
    model: SpaceStateModel,
}

impl NewModel {
    pub fn new() -> Self {
        let mut model = SpaceStateModel::new(2, 1, 2);
        model.set_mat_a(&[0.0, 1.0, -1.0, -1.0]).unwrap();
        model.set_mat_b(&[0.0, 1.0]).unwrap();
        model.init_state(&[0.0, 0.0]).unwrap();

        NewModel {
            model : model
        }
    }
}

impl Model for NewModel {
    fn slopefunc(&self, x: &DMatrix<f64>) -> DMatrix<f64> {
        self.model.slopefunc(x)
    }

    fn get_state_info(&self) -> Vec<String> {
        self.model.get_state_info()
    }

    fn set_state(&mut self, newstate: DMatrix<f64>) {
        self.model.set_state(newstate);
    }

    fn get_state(&self) -> &DMatrix<f64> {
        &self.model.get_state()
    }

    fn calc_nextstate(&mut self, delta_t : f64, solvertype: &SolverType) { 
        self.model.set_u(&[1.0]).unwrap();

        self.model.calc_nextstate(delta_t, solvertype);
    }
}

struct RLCCircuit {
    model: SpaceStateModel,
}

impl RLCCircuit {
    pub fn new(r: f64, l: f64, c: f64) -> Self {
        let mut model = SpaceStateModel::new(2, 1, 2);
        model.set_mat_a(&[-r / l, -1.0 / (l * c), 1.0, 0.0]).unwrap();
        model.set_mat_b(&[1.0 / l, 0.0]).unwrap();
        model.init_state(&[0.0, 0.0]).unwrap();

        Self {
            model : model
        }
    }
}

impl Model for RLCCircuit {
    fn slopefunc(&self, x: &DMatrix<f64>) -> DMatrix<f64> {
        self.model.slopefunc(x)
    }

    fn get_state_info(&self) -> Vec<String> {
        vec!["i".to_string(), "q".to_string()]
    }

    fn set_state(&mut self, newstate: DMatrix<f64>) {
        self.model.set_state(newstate);
    }

    fn get_state(&self) -> &DMatrix<f64> {
        &self.model.get_state()
    }

    fn calc_nextstate(&mut self, delta_t : f64, solvertype: &SolverType) { 
        self.model.set_u(&[1.0]).unwrap();

        self.model.calc_nextstate(delta_t, solvertype);
    }
}

fn main() {
    
    let model = NewModel::new();
    let rlc = RLCCircuit::new(10.0, 0.3e-3, 0.1e-6);    // http://programmer-life.net/tfunc_040.html　のパラメータ
    let rlc2 = RLCCircuit::new(10.0, 100e-3, 100e-6);   

    let mut solver = Simulator::<NewModel>::new(10.0, 0.01, SolverType::RungeKutta, model);
    solver.run_sim();
    solver.export_sim("./test.csv");

    let mut rlcsim = Simulator::<RLCCircuit>::new(0.001, 0.0000001, SolverType::RungeKutta, rlc);
    rlcsim.run_sim();
    rlcsim.export_sim("./rlc.csv");

    let mut rlcsim2 = Simulator::<RLCCircuit>::new(0.2, 0.00001, SolverType::RungeKutta, rlc2);
    rlcsim2.run_sim();
    rlcsim2.export_sim("./rlc2.csv");
    

}
