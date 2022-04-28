extern crate nalgebra as na;
use na::{U2, U3, Dynamic, ArrayStorage, VecStorage, Matrix, OMatrix, DMatrix};

use std::fmt;

mod simtools;
use simtools::{simsolver};
use simsolver::{*};
use simmodel::{*};

struct NewModel {
    model: SpaceStateModel,
}

impl NewModel {
    pub fn new() -> Self {
        let mut model = SpaceStateModel::new(2, 1, 1);
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

    fn get_signals_info(&self) -> Vec<String> {
        self.model.get_signals_info()
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

    fn get_allsignals(&self) -> Vec<f64> { 
        self.model.get_allsignals()
    }
}

impl fmt::Display for NewModel {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        self.model.fmt(f)
    }
}


struct RLCCircuit {
    model: SpaceStateModel,
}

impl RLCCircuit {
    pub fn new(r: f64, l: f64, c: f64, init_i: f64, init_q: f64) -> Self {
        let mut model = SpaceStateModel::new(2, 1, 2);
        model.set_mat_a(&[-r / l, -1.0 / (l * c), 1.0, 0.0]).unwrap();
        model.set_mat_b(&[1.0 / l, 0.0]).unwrap();
        model.set_mat_c(&[r, 0.0, 0.0, 1.0 / c]).unwrap();
        model.init_state(&[init_i, init_q]).unwrap();

        Self {
            model : model
        }
    }
}

impl Model for RLCCircuit {
    fn slopefunc(&self, x: &DMatrix<f64>) -> DMatrix<f64> {
        self.model.slopefunc(x)
    }

    fn get_signals_info(&self) -> Vec<String> {
        vec!["v", "i", "q", "Vr", "Vc"]
            .iter().map(|x| x.to_string()).collect::<Vec<String>>()
    }

    fn set_state(&mut self, newstate: DMatrix<f64>) {
        self.model.set_state(newstate);
    }

    fn get_state(&self) -> &DMatrix<f64> {
        &self.model.get_state()
    }

    fn calc_nextstate(&mut self, delta_t : f64, solvertype: &SolverType) {
        let x = self.get_state();
        // let f = DMatrix::from_vec(1, 2, vec![1.0, -9.7e3]); // pole = [-50, -60]
        let f = DMatrix::from_vec(1, 2, vec![-6.00e+00, -9.95e+03]); // pole = [-20 + 10j, -20 - 10j]
        let u = -f * x;
        
        self.model.set_u(&u.iter().map(|x| *x).collect::<Vec<f64>>()).unwrap();
        //self.model.set_u(&[0.0]).unwrap();

        self.model.calc_nextstate(delta_t, solvertype);
    }

    fn get_allsignals(&self) -> Vec<f64> { 
        self.model.get_allsignals()
    }
}

impl fmt::Display for RLCCircuit {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        self.model.fmt(f)
    }
}

fn main() {
    /*
    let model = NewModel::new();
    let rlc = RLCCircuit::new(10.0, 0.3e-3, 0.1e-6, 0.0, 0.0);    // http://programmer-life.net/tfunc_040.html　のパラメータ
    let rlc2 = RLCCircuit::new(10.0, 100e-3, 100e-6, 0.1, 0.0002);   
    println!("{}", model);
    println!("{}", rlc2);

    let mut solver = Simulator::<NewModel>::new(10.0, 0.01, SolverType::RungeKutta, model);
    solver.run_sim();
    solver.export_sim("./test.csv");

    let mut rlcsim = Simulator::<RLCCircuit>::new(0.001, 0.0000001, SolverType::RungeKutta, rlc);
    rlcsim.run_sim();
    rlcsim.export_sim("./rlc.csv");

    let mut rlcsim2 = Simulator::<RLCCircuit>::new(1.0, 0.00001, SolverType::RungeKutta, rlc2);
    rlcsim2.run_sim();
    rlcsim2.export_sim("./rlc2.csv");

    rlcsim2.timeplot("rlc4", (1000, 300));*/

    let mut model = SpaceStateModel::from_tf(&[3.0, 1.0, 1.0, 5.0, 4.0], &[2.0, 2.0, 3.0, 4.0, 5.0]).unwrap();
    println!("{}", model);
    model.set_u(&[1.0]);
    let mut tfsim = Simulator::<SpaceStateModel>::new(1.0, 0.001, SolverType::RungeKutta, model);
    tfsim.run_sim();
    tfsim.export_sim("./tfsim.csv");
    tfsim.timeplot("tfsim", (1000, 300));
    
}
